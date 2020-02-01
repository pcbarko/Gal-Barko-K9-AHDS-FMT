
#Initialize session ===

library(phyloseq)
library(tidyverse)
library(ggplot2)
library(ape)
library(dplyr)
library(gridExtra)
library(kableExtra)
library(grid)
library(lme4)
library(lmerTest)
library(emmeans)
library(vegan)
library(microbiome)
library(ranacapa)
library(ggpubr)
library(psych)

source("Functions_Microbiome.R")

#Load and Inspect Data

##The phyloseq object containing the raw, unfiltered data is loaded for analysis of alpha diversity. 
physeq <- readRDS("FMT_physeq_raw.RDS") 

physeq
summarize_phyloseq(physeq)

## The filtered and agglomerated relative abundnace phyloseq object for analysis of beta deversity/ordination.
physeq_glom_ra <- readRDS("FMT_physeq_glom_ra.RDS") 

physeq_glom_ra
summarize_phyloseq(physeq_glom_ra)

#Analysis of microbiota alpha diversity ====

alpha_diversity <- estimate_richness(physeq)

meta <- as.data.frame(sample_data(physeq))

rownames(alpha_diversity) <- rownames(sample_data(physeq))

shannon <-estimate_richness(physeq, measures = "Shannon")

shannon <- cbind(shannon, meta)

#Assess normality of shannon diversity indices

hist(alpha_diversity$Shannon, main="Shannon diversity", xlab="", breaks = 15)

shapiro.test(alpha_diversity$Shannon)

#The histogram looks a bit skewed but the Shapiro-Wilk test is not violated. 
#We will consider the distribution of Shannon diversity indices to be approximately normal.

#Compare Shannon diversity between dogs with AHDS (pre-treatment) and Healthy donors

#subset by admission time
shannon_admit <- shannon[shannon$Collection_Time == "Admission", ]

t.test(shannon_admit$Shannon ~ shannon_admit$Diagnosis)
#there is a significant difference between healthy dogs and dogs with AHDS before treatmen

#Linear model to assess differences in alpha diversity due to treatment group and time
model_shannon <- lmer(formula = Shannon ~ Treatment_Group*Collection_Time + (1|DogID),
              data = shannon)

summary(model_shannon)

#residuals plot
qqnorm(residuals(model_shannon))

#pairwise comparisons
emm_shannon <- emmeans(model_shannon, pairwise ~ Collection_Time | Treatment_Group)

summary(emm_shannon, adjust = "fdr")

emm_shannon2 <- emmeans(model_shannon, pairwise ~ Treatment_Group | Collection_Time)

summary(emm_shannon2, adjust = "fdr")

#Summary statistics for Shannon diversity ====

#baseline data: AHDS vs Healthy 
admit <- shannon[shannon$Collection_Time == "Admission", ] 

stats <- describeBy(admit$Shannon, group = admit$Diagnosis)

descStats <- rbind(stats$Donor, stats$AHDS)

rownames(descStats) <- c("Healthy", "AHDS")

descStats <- as.data.frame(descStats)

descStats[, -c(1, 6:7, 11:13)]

#group-wise at admission
stats2 <- describeBy(admit$Shannon, group = admit$Treatment_Group)

descStats2 <- rbind(stats2$Donor, stats2$Saline, stats2$FMT)

rownames(descStats2) <- c("Donor", "Saline", "FMT")

descStats2 <- as.data.frame(descStats2)

#group-wise at discharge
dc <- shannon[shannon$Collection_Time == "Discharge", ]

stats3 <- describeBy(dc$Shannon, group = dc$Treatment_Group)

descStats3 <- rbind(stats3$Donor, stats3$Saline, stats3$FMT)

rownames(descStats3) <- c("Donor", "Saline", "FMT")

descStats3 <- as.data.frame(descStats3)

#group-wise at 30 days
d30 <- shannon[shannon$Collection_Time == "30days", ]

stats4 <- describeBy(d30$Shannon, group = d30$Treatment_Group)

descStats4 <- rbind(stats4$Donor, stats4$Saline, stats2$FMT)

rownames(descStats4) <- c("Donor", "Saline", "FMT")

descStats4 <- as.data.frame(descStats4)

#Plots of Shannon diversity ====

gplot <- ggplot(shannon_admit, aes(y = Shannon, x = Diagnosis, fill = Diagnosis))

gplot + geom_boxplot(outlier.colour="red", outlier.shape=8,
                     outlier.size=4) +
  ggtitle("Shannon Diversity at Admission (Pre-Treatment)") +
  xlab("Diagnosis") +
  ylab("Shannon Diversity Index") +
  stat_summary(fun.y=mean, geom="point", shape=4, size=4) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5))


gplot2 <- ggplot(shannon, aes(y = Shannon, x = Treatment_Group, fill = Treatment_Group))

gplot2 + geom_boxplot(outlier.colour="red", outlier.shape=8,
                      outlier.size=3) +
  ggtitle("Shannon Diversity by Treatment Group and Collection Time") +
  xlab("Treatment Group") +
  ylab("Shannon Diversity Index") +
  facet_wrap(~Collection_Time) +
  stat_summary(fun.y=mean, geom="point", shape=4, size=3) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5))


#Beta Disversity and Ordination ====

physeq.ord <- ordinate(physeq_glom_ra, "PCoA", "bray")

plot_ordination(physeq_glom_ra, physeq.ord, type="scree", 
                title="Bray-Curtis - Scree Plot") +
  coord_cartesian(xlim = c(0, 20))

plot_ordination(physeq_glom_ra, 
                physeq.ord, 
                type="samples", 
                color="Diagnosis",
                shape = "Collection_Time",
                axes = c(1,2)) + 
  geom_point(size = 5) +
  geom_text(aes(label=as.character(physeq_glom_ra@sam_data$Treatment_Group)), hjust = 1.3, size=3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Principal Coordinate Analysis: \nBray-Curtis Dissimilarity (Axes 1/2)") +
  theme(plot.title = element_text(hjust = 0.5))

plot_ordination(physeq_glom_ra, 
                physeq.ord, 
                type="samples", 
                color="Diagnosis",
                shape = "Collection_Time",
                axes = c(1,3)) + 
  ggtitle("PCoA of Bray-Curtis Dissimilarity (Axes 1/3)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_point(size = 5) +
  geom_text(aes(label=as.character(physeq_glom_ra@sam_data$Treatment_Group)), hjust = 1.3, size=3) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Principal Coordinate Analysis: \nBray-Curtis Dissimilarity (Axes 1/3)") +
  theme(plot.title = element_text(hjust = 0.5))

#is the observed separation by diagnosis statistically significant?

otu <- as.data.frame(t(physeq_glom_ra@otu_table))

otu$Diagnosis <- meta$Diagnosis

otu$Colletion_Time <- meta$Collection_Time

otu$Treatment_Group <- meta$Treatment_Group

anosim(otu[, 1:87], otu$Diagnosis, permutations = 999, distance = "bray", strata = NULL,
       parallel = getOption("mc.cores"))

#Is there significant clustering by collection time?
anosim(otu[, 1:87], otu$Colletion_Time, permutations = 999, distance = "bray", strata = NULL,
       parallel = getOption("mc.cores"))

#Is there significant clustering by treatment group
anosim(otu[, 1:87], otu$Treatment_Group, permutations = 999, distance = "bray", strata = NULL,
       parallel = getOption("mc.cores"))







