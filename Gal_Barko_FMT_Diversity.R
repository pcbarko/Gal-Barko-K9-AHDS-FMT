
#Initialize session ===

setwd("~/Box/Gal_Barko_FMT_Final/Manuscript folder/Data and R Code")

library(phyloseq)
library(tidyverse)
library(ggplot2)
library(ape)
library(dplyr)
library(gridExtra)
library(grid)
library(lme4)
library(lmerTest)
library(emmeans)
library(vegan)
library(microbiome)
library(ranacapa)
library(ggpubr)
library(psych)
library(rstatix)

#Load and Inspect Data ====

##The phyloseq object containing the raw, unfiltered data is loaded for analysis of alpha diversity. 
physeq <- readRDS("FMT_physeq_raw.RDS") 

## The filtered and agglomerated relative abundnace phyloseq object for analysis of beta deversity/ordination.
physeq_glom_ra <- readRDS("FMT_physeq_glom_ra.RDS") 

#Alpha diversity ====

alpha_diversity <- estimate_richness(physeq)

meta <- as.data.frame(sample_data(physeq))

rownames(alpha_diversity) <- rownames(sample_data(physeq))

shannon <-estimate_richness(physeq, measures = "Shannon")

shannon <- cbind(shannon, meta)

#Assess normality of shannon diversity indices

hist(alpha_diversity$Shannon, main="Shannon diversity", xlab="", breaks = 15)

shapiro.test(alpha_diversity$Shannon)

##The histogram looks a bit skewed but the Shapiro-Wilk test is not violated. 
##Consider the distribution of Shannon diversity indices to be approximately normal.

#Compare Shannon diversity between dogs with AHDS (pre-treatment) and Healthy donors

#subset by admission time
shannon_admit <- shannon[shannon$Collection_Time == "Admission", ]

t.test(shannon_admit$Shannon ~ shannon_admit$Diagnosis)
##there is a significant difference between healthy dogs and dogs with AHDS before treatmen

#Linear model to assess differences in alpha diversity due to treatment group and time

model_shannon <- lmer(formula = Shannon ~ Treatment_Group*Collection_Time + (1|DogID),
              data = shannon)

#residuals plot
qqnorm(residuals(model_shannon))

#pairwise comparisons
emm_shannon <- emmeans(model_shannon, pairwise ~ Collection_Time | Treatment_Group)

summary(emm_shannon, adjust = "fdr")

emm_shannon2 <- emmeans(model_shannon, pairwise ~ Treatment_Group | Collection_Time)

summary(emm_shannon2, adjust = "fdr")


#Plots of Shannon diversity ====

#Baseline admission samples AHDS vs Healtyh Donor

gplot <- ggboxplot(shannon_admit, x = "Diagnosis", y = "Shannon", 
                    fill = "Diagnosis", 
                    add = "jitter") + 
  xlab("") +
  ylab("Shannon Diversity Index") + 
  theme(legend.position="bottom") +
  stat_compare_means(method = "t.test",
                     comparisons = list(c("Donor", "AHDS")), 
                     label = "p.signif")

gplot

#Post-FMT faceted by Collection_Time
gplot2 <- ggboxplot(shannon, x = "Treatment_Group", y = "Shannon", 
                   fill = "Treatment_Group", 
                   add = "jitter", 
                   facet.by = "Collection_Time", 
                   palette = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  labs(fill = "Treatment Group", x= "", y = "Shannon Diversity Index") +
  theme(legend.position="bottom")

gplot2


emdf <- as.data.frame(summary(emm_shannon2, adjust = "fdr"))
contrasts <- strsplit(as.character(emdf$contrasts.contrast), split = "-", fixed = T)

group1 <- trimws(sapply(contrasts, `[`, 1))
group2 <- trimws(sapply(contrasts, `[`, 2))

Collection_Time <- emdf$emmeans.Collection_Time
pval <- emdf$contrasts.p.value
y.position <- c(4.1, 4.3, 4.5, 4.1, 4.3, 4.5, 4.1, 4.3, 4.5)

emstats <- data.frame(Collection_Time, group1, group2, pval, y.position)
emstats$signif <- ifelse(emstats$pval < 0.05, "*", "ns")

gplot2 <- gplot2 + stat_pvalue_manual(emstats, label = "signif") 

gplot2

#Post-FMT faceted by Treatment_Group
gplot3 <- ggboxplot(shannon, x = "Collection_Time", y = "Shannon", 
                    fill = "Collection_Time", 
                    add = "jitter", 
                    facet.by = "Treatment_Group") +
  labs(fill = "Collection Time", x= "", y = "Shannon Diversity Index") +
  theme(legend.position="bottom")

emdf2 <- as.data.frame(summary(emm_shannon, adjust = "fdr"))
contrasts2 <- strsplit(as.character(emdf2$contrasts.contrast), split = "-", fixed = T)

group1 <- trimws(sapply(contrasts2, `[`, 1))
group2 <- trimws(sapply(contrasts2, `[`, 2))

Treatment_Group <- emdf2$emmeans.Treatment_Group
pval2 <- round(emdf2$contrasts.p.value, 3)
y.position <- c(4.1, 4.3, 4.5, 4.1, 4.3, 4.5, 4.1, 4.3, 4.5)

emstats2 <- data.frame(Treatment_Group, group1, group2, pval2, y.position)
emstats2$signif <- ifelse(emstats2$pval2 < 0.05, "*", "ns")

gplot3 <- gplot3 + stat_pvalue_manual(emstats2, label = "signif")

gplot3

##Output

tiff("shannon.tiff", height = 50, width = 30, units = "cm", res = 600)
ggarrange(ggarrange(gplot, ncol = 2), 
          gplot3, 
          gplot2,
          nrow = 3,
          labels = c("A", "B", "C")
) 
dev.off()

#Beta Disversity and Ordination ====

##AHDS vs Donor at Admission

physeq.ord <- ordinate(subset_samples(physeq_glom_ra, 
                                      Collection_Time == "Admission"),
                                      "PCoA", "bray")

physeq.ord2 <- ordinate(physeq_glom_ra, "PCoA", "bray")

plot_ordination(physeq_glom_ra, physeq.ord, type="scree", 
                title="Bray-Curtis - Scree Plot") +
  coord_cartesian(xlim = c(0, 20))

pcoa1 <- plot_ordination(
  subset_samples(physeq_glom_ra, Collection_Time == "Admission"), 
                physeq.ord, 
                type="samples", 
                color="Diagnosis",
                axes = c(1,2)) + 
  geom_point(size = 5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position="bottom")

pcoa1

##is the observed separation by diagnosis statistically significant at baseline (admission)?
###use PERMANOVA 

otu_admit <- otu_table(subset_samples(physeq_glom_ra, Collection_Time == "Admission"))

otu_admit <- as.data.frame(t(otu_admit))

veg.dist_admit <- vegdist(otu_admit, distance = "bray")

admit <- data.frame(sample_data(subset_samples(physeq_glom_ra, Collection_Time == "Admission")))

adonis(veg.dist_admit ~ Diagnosis, data = admit, permutations = 1000)
##pval = 0.003996

pcoa1 <- pcoa1 + annotate(geom = "text", x = 0.23, y = -0.3, 
                          label = "Diagnosis P = 0.003") + 
  theme(plot.margin=unit(c(5,5,0,10),"mm")) + 
  geom_vline(xintercept = 0, 
             linetype="dotted", 
             color = "grey", 
             size=1) + 
  geom_hline(yintercept = 0, 
             linetype="dotted", 
             color = "grey", 
             size=1)

pcoa1

##Treatment Groups vs Collection Time

pcoa2 <- plot_ordination(physeq_glom_ra, 
                physeq.ord2, 
                type="samples", 
                color="Treatment_Group",
                shape = "Collection_Time",
                axes = c(1,2)) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_point(size = 5) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position="bottom") + labs(color = "Treatment Group", 
                                         shape = "Collection Time") + 
  geom_vline(xintercept = 0, 
             linetype="dotted", 
             color = "grey", 
             size=1) + 
  geom_hline(yintercept = 0, 
             linetype="dotted", 
             color = "grey", 
             size=1)

##are the observed clusters by treatment group and collection statistically significant?

otu <- otu_table(subset_samples(physeq_glom_ra, Diagnosis == "AHDS"))

otu <- as.data.frame(t(otu))

meta <- data.frame(sample_data(subset_samples(physeq_glom_ra, Diagnosis == "AHDS")))

veg.dist <- vegdist(otu, distance = "bray")

adonis(veg.dist ~ Treatment_Group * Collection_Time + DogID, data = meta, permutations = 1000)
###treatment_group p = 0.008991
###Collection_Time p = 0.048941

pcoa2 <- pcoa2 + annotate(geom = "text", x = 0.23, y = -0.47, 
                          label = "Treatment Group P = 0.009 \n Collection Time P = 0.049") + 
  theme(plot.margin=unit(c(5,5,0,10),"mm"))

pcoa2

##Output
tiff("PCoA.tiff", height = 20, width = 42, units = "cm", res = 600)
ggarrange(pcoa1,
          pcoa2,
          ncol = 2,
          labels = c("A", "B")
) 
dev.off()

#Shared OTUs and Divergence in Donor-Recipient Pairs ====
###There is probably a more elegent way to do this, but this works

##generate statistics from subsets of phyloseq object
physeq2 <- subset_samples(physeq, Treatment_Group == "Donor" | Treatment_Group == "FMT")

##Admission Samples
physeq2_admit <- subset_samples(physeq2, Collection_Time == "Admission")

sample_data(physeq2_admit)

physeq2_admit1 <- subset_samples(physeq2_admit, DogID == "1" | DogID == "2")
shared_1v2_admit <-  ntaxa(filter_taxa(physeq2_admit1, function(x) sum(x >= 1) == (2), TRUE))
dist_1v2_admit <- distance(physeq2_admit1, method = "jsd", type = "samples")


physeq2_admit2 <- subset_samples(physeq2_admit, DogID == "7" | DogID == "6")
shared_7v6_admit <-  ntaxa(filter_taxa(physeq2_admit2, function(x) sum(x >= 1) == (2), TRUE))
dist_7v6_admit <- distance(physeq2_admit2, method = "jsd", type = "samples")


physeq2_admit3 <- subset_samples(physeq2_admit, DogID == "9" | DogID == "8")
shared_9v8_admit <-  ntaxa(filter_taxa(physeq2_admit3, function(x) sum(x >= 1) == (2), TRUE))
dist_9v8_admit <- distance(physeq2_admit3, method = "jsd", type = "samples")


physeq2_admit4 <- subset_samples(physeq2_admit, DogID == "11" | DogID == "12")
shared_11v12_admit <-  ntaxa(filter_taxa(physeq2_admit4, function(x) sum(x >= 1) == (2), TRUE))
dist_11v12_admit <- distance(physeq2_admit4, method = "jsd", type = "samples")


##Discharge Samples
physeq2_DC <- subset_samples(physeq2, Collection_Time == "Discharge")

physeq2_DC1 <- subset_samples(physeq2_DC, DogID == "1" | DogID == "2")
shared_1v2_DC <-  ntaxa(filter_taxa(physeq2_DC1, function(x) sum(x >= 1) == (2), TRUE))
dist_1v2_DC <- distance(physeq2_DC1, method = "jsd", type = "samples")


physeq2_DC2 <- subset_samples(physeq2_DC, DogID == "7" | DogID == "6")
shared_7v6_DC <-  ntaxa(filter_taxa(physeq2_DC2, function(x) sum(x >= 1) == (2), TRUE))
dist_7v6_DC <- distance(physeq2_DC2, method = "jsd", type = "samples")


physeq2_DC3 <- subset_samples(physeq2_DC, DogID == "9" | DogID == "8")
shared_9v8_DC <-  ntaxa(filter_taxa(physeq2_DC3, function(x) sum(x >= 1) == (2), TRUE))
dist_9v8_DC <- distance(physeq2_DC3, method = "jsd", type = "samples")


physeq2_DC4 <- subset_samples(physeq2_DC, DogID == "11" | DogID == "12")
shared_11v12_DC <-  ntaxa(filter_taxa(physeq2_DC4, function(x) sum(x >= 1) == (2), TRUE))
dist_11v12_DC <- distance(physeq2_DC4, method = "jsd", type = "samples")


##30 Day Samples
physeq2_30 <- subset_samples(physeq2, Collection_Time == "30days")

physeq2_30_1 <- subset_samples(physeq2_30, DogID == "1" | DogID == "2")
shared_1v2_30 <-  ntaxa(filter_taxa(physeq2_30_1, function(x) sum(x >= 1) == (2), TRUE))
dist_1v2_30 <- distance(physeq2_30_1, method = "jsd", type = "samples")


physeq2_30_2 <- subset_samples(physeq2_30, DogID == "7" | DogID == "6")
shared_7v6_30 <-  ntaxa(filter_taxa(physeq2_30_2, function(x) sum(x >= 1) == (2), TRUE))
dist_7v6_30 <- distance(physeq2_30_2, method = "jsd", type = "samples")


physeq2_30_3 <- subset_samples(physeq2_30, DogID == "9" | DogID == "8")
shared_9v8_30 <-  ntaxa(filter_taxa(physeq2_30_3, function(x) sum(x >= 1) == (2), TRUE))
dist_9v8_30 <- distance(physeq2_30_3, method = "jsd", type = "samples")


physeq2_30_4 <- subset_samples(physeq2_30, DogID == "11" | DogID == "12")
shared_11v12_30 <-  ntaxa(filter_taxa(physeq2_30_4, function(x) sum(x >= 1) == (2), TRUE))
dist_11v12_30 <- distance(physeq2_30_4, method = "jsd", type = "samples")


##Make df

pair1 <- data.frame("Pair" = rep("Pair 1", 3), 
                    "Timepoint" = c("Admission", "Discharge", "30 Days"), 
                    nShared = c(shared_1v2_admit, shared_1v2_DC, shared_1v2_30),
                    jsd = c(dist_1v2_admit, dist_1v2_DC, dist_1v2_30))

pair2 <- data.frame("Pair" = rep("Pair 2", 3), 
                    "Timepoint" = c("Admission", "Discharge", "30 Days"), 
                    nShared = c(shared_7v6_admit, shared_7v6_DC, shared_7v6_30),
                    jsd = c(dist_7v6_admit, dist_7v6_DC, dist_7v6_30))

pair3 <- data.frame("Pair" = rep("Pair 3", 3), 
                    "Timepoint" = c("Admission", "Discharge", "30 Days"), 
                    nShared = c(shared_9v8_admit, shared_9v8_DC, shared_9v8_30),
                    jsd = c(dist_9v8_admit, dist_9v8_DC, dist_9v8_30))

pair4 <- data.frame("Pair" = rep("Pair 4", 3), 
                    "Timepoint" = c("Admission", "Discharge", "30 Days"), 
                    nShared = c(shared_11v12_admit, shared_11v12_DC, shared_11v12_30),
                    jsd = c(dist_11v12_admit, dist_11v12_DC, dist_11v12_30))

shared_OTU <- rbind(pair1, pair2, pair3, pair4)

shared_OTU$Timepoint <- factor(shared_OTU$Timepoint, levels = c("Admission", "Discharge", "30 Days"))

##Plots of Shared OTUs and Divergence

p1 <- ggplot(data = shared_OTU, aes(x = shared_OTU$Timepoint, 
                                    y = shared_OTU$nShared, 
                                    color = shared_OTU$Pair, 
                                    group = shared_OTU$Pair)) +
  geom_point() +
  geom_line()+
  xlab("") +
  ylab("Number of Shared OTUs") +
  labs(color="Donor-Recipent Pair") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(plot.margin=unit(c(5,5,0,10),"mm"))


p2 <- ggplot(data = shared_OTU, aes(x = shared_OTU$Timepoint, 
                                    y = shared_OTU$jsd, 
                                    color = shared_OTU$Pair, 
                                    group = shared_OTU$Pair)) +
  geom_point() +
  geom_line()+
  xlab("") +
  ylab("Jensen-Shannon Divergence") +
  labs(color="Donor-Recipent Pair") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(plot.margin=unit(c(5,5,0,10),"mm"))


physeq.ord2 <- ordinate(physeq2, "PCoA", "jsd")

plot_ordination(physeq2, physeq.ord2, type="scree", 
                title="Bray-Curtis - Scree Plot") +
  coord_cartesian(xlim = c(0, 15))

p3 <- plot_ordination(physeq2, 
                      physeq.ord2, 
                      type="samples", 
                      color="Treatment_Group",
                      shape = "Collection_Time",
                      axes = c(1,2)) + 
  geom_point(size = 5) +
  theme_bw() + 
  ggtitle("Principal Coordinate Analysis: \nJensen-Shannon Divergence (Axes 1/2)") +
  theme(plot.title = element_text(hjust = 0.5))



##Statistical Comparisons

shared_OTU %>%
  group_by(Timepoint) %>%
  get_summary_stats(nShared)

shared_OTU %>%
  group_by(Timepoint) %>%
  shapiro_test(nShared)

ggqqplot(shared_OTU, "nShared", facet.by = "Timepoint")

res.aov <- anova_test(data = shared_OTU[, 1:3], dv = nShared, wid = Pair, within = Timepoint)
get_anova_table(res.aov)

pwc <- shared_OTU %>%
  pairwise_t_test(
    nShared ~ Timepoint, paired = TRUE,
    p.adjust.method = "BH"
  )
pwc

shared_OTU %>%
  group_by(Timepoint) %>%
  get_summary_stats(jsd)

shared_OTU %>%
  group_by(Timepoint) %>%
  shapiro_test(jsd)

ggqqplot(shared_OTU, "jsd", facet.by = "Timepoint")

res.aov2 <- anova_test(data = shared_OTU[, c(1,2,4)], dv = jsd, wid = Pair, within = Timepoint)
get_anova_table(res.aov2)

pwc2 <- shared_OTU %>%
  pairwise_t_test(
    jsd ~ Timepoint, paired = TRUE,
    p.adjust.method = "BH"
  )
pwc2

##Output

pwc$y.position <- c(150, 160, 170)
p1 <- p1 + stat_pvalue_manual(pwc, label = "p.adj.signif")

pwc2$y.position <- c(0.72, 0.77, 0.82)
p2 <- p2 + stat_pvalue_manual(pwc2, label = "p.adj.signif")


tiff("jsd_shared_plot.tiff", height = 15, width = 30, units = "cm", res = 600)
ggarrange(p1, p2, 
          ncol = 2,
          common.legend = TRUE, legend = "bottom",
          labels = c("A", "B")
) 
dev.off()



