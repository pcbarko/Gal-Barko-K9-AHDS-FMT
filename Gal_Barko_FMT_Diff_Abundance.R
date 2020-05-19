
#Initialize session ====
library(DEFormats)
library(edgeR)
library(Glimma)
library(gplots)
library(statmod)
library(dplyr)
library(scales)
library(ggplot2)
library(gridExtra)
library(grid)
library(ggpubr)


source("summarizeFit.R")

#load data ====
physeq_glom <- readRDS("FMT_physeq_glom.RDS") 

#create genus_species annotation for plotting
taxa <- read.csv("FMT_taxa_table.csv", row.names = 1)
taxa$gs <- paste(taxa$Genus, taxa$Species, sep = "_")

physeq_glom@sam_data$Group <- as.factor(paste(physeq_glom@sam_data$Treatment_Group, 
                                              physeq_glom@sam_data$Collection_Time, 
                                          sep = "_"))

#use DESeq2 to load
de <- phyloseq_to_deseq2(physeq_glom, ~ Diagnosis)

dge <- as.DGEList(de)

dge <- dge[,order(dge$samples$Treatment_Group, dge$samples$DogID,dge$samples$Collection_Time)]

dge$samples$ind.n <- rep(1:4, each = 3)

dge <- calcNormFactors(dge)

dge$samples$Donor_ID[dge$samples$Treatment_Group == "Donor"] <- 
  dge$samples$DogID[dge$samples$Treatment_Group == "Donor"]


#get logCPM values

logCPM <- cpm(dge, log = TRUE, prior.count = 0.5)

plotDensities(logCPM, legend = F)

#Assess strength of random effects: DogID and DonorID ====

#assess dog strength
design1 <- model.matrix(~dge$samples$Collection_Time)

dog.cor <- duplicateCorrelation(logCPM, design = design1, block = dge$samples$DogID )
dog.cor$cor
#0.2493754

#assess donor strength
temp <- dge$samples$Treatment_Group != "Saline" & 
  dge$samples$Collection_Time != "Admission"

design2 <- model.matrix(~dge$samples$Collection_Time[temp, drop = TRUE])

donor.cor <- duplicateCorrelation(logCPM[,temp], design = design2, 
                                  block = dge$samples$Donor_ID[temp] )
donor.cor$cor
#0.04940949 - low so can ignore


#model as in limmaUsersGuide() section 9.7 ====

designGroup <- model.matrix(~ 0 + Group, data = dge$samples)
colnames(designGroup) <- levels(dge$samples$Group)

corfit <- duplicateCorrelation(logCPM, design = designGroup, block = dge$samples$DogID)
corfit$cor
#0.2281621

fit <- lmFit(logCPM, design = designGroup, block = dge$samples$DogID, correlation = corfit$cor)

summary(fit)

colnames(designGroup)

#contrasts of interest:

cont.matrix <- makeContrasts(ad_AHDSvsDonor = (FMT_Admission+Saline_Admission)/2 - Donor_Admission,
                             dis_FMTvsSaline = FMT_Discharge - Saline_Discharge,
                             dis_FMTvsDonor = FMT_Discharge - Donor_Discharge,
                             dis_SalinevsDonor = Saline_Discharge - Donor_Discharge,
                             d30_FMTvsSaline = FMT_30days - Saline_30days,
                             d30_FMTvsDonor = FMT_30days - Donor_30days,
                             d30_SalinevsDonor = Saline_30days - Donor_30days,
                             levels = designGroup)

fit2 <- eBayes(contrasts.fit(fit, cont.matrix))
codes.res <- decideTests(fit2, p.value = 0.2)

summary(codes.res)

topTable(fit2, coef = 3)

layout(matrix(1:8, 2, 4))
for(i in 1:ncol(cont.matrix))
  hist(fit2$p.value[,i], 20, xlim = c(0,1), main = colnames(cont.matrix)[i])
dev.off()

layout(matrix(1:8, 2, 4))
for(i in 1:ncol(cont.matrix))
  hist(p.adjust(fit2$p.value[,i], "fdr"), 20, xlim = c(0,1), main = colnames(cont.matrix)[i])
dev.off()

sum(rowSums(codes.res !=0)>0)

topTable(fit2, coef = 6)

#Heatmap ====

col.pan <- colorpanel(100, "blue", "white", "red")

heat.h1 <- t(scale(t(logCPM[rowSums(codes.res !=0)>0,])))

rownames(heat.h1) <- taxa$gs[match(rownames(heat.h1), rownames(taxa))]

colnames(heat.h1) <- paste(dge$samples$Group, dge$samples$DogID, sep = "_")

cluster.x.h1 <- hclust(dist(heat.h1))
cluster.y.h1 <- hclust(dist(t(heat.h1)))

temp <- factor(dge$samples$Treatment_Group, levels = c("Saline","FMT","Donor"))
ind <- order(dge$samples$Collection_Time, temp)

heat <- as.data.frame(t(heat.h1[, c(1,4,7,10,13,16,19,22,25,28,31,34,2,5,8,11,14,17,20,23,26,29,33,35,3,6,9,12,15,18,21,24,27,30,33,36)]))

heat$'Treatment Group' <- factor(rep(c(rep("Donor", 4), rep("Saline", 4), rep("FMT", 4)), 3), 
                                 levels = c("Donor", "Saline", "FMT"))
levels(heat$`Treatment Group`)

heat$'Collection Time' <- factor(c(rep("Admission", 12), rep("Discharge", 12), rep("30 Days", 12)), 
                                 levels = c("Admission", "Discharge", "30 Days"))
levels(heat$`Collection Time`)

heat_col <- heat[, c(30, 29)]

heatmap <- pheatmap(t(heat[, -c(29,30)]),
         color = colorRampPalette(c("navy", "white", "firebrick3"))(10),
         cluster_cols = F, 
         cutree_cols = 3, 
         scale = "row",
         gaps_col = c(4, 8, 12, 16, 20, 24, 28, 32, 36), 
         annotation_col = heat_col,
         show_colnames = F
         )

res.out <- summarizeFit(fit2, calcFC = F)

#plot results ====
res.out$taxa <- taxa$gs[match(rownames(res.out), rownames(taxa))]

res <- res.out[, c(23, 1:22)]


res_admit <- res[, c(1, 3:5)]
res_admit <- res_admit[res_admit$rawP.ad_AHDSvsDonor < 0.05, ]
res_admit$Comparison <- as.factor("AHDS vs. Donor")
res_admit$Time <- "Admission"
names(res_admit) <- c("Taxa", "log2FC", "P-Value", "FDR", "Comparison", "Time")

kable <- kable(res_admit[, 1:4]) %>%
  kable_styling("striped", full_width = F) 


res_DC_FMTvSal <- res[, c(1, 6:8)]
res_DC_FMTvSal <- res_DC_FMTvSal[res_DC_FMTvSal$rawP.dis_FMTvsSaline < 0.05, ]
res_DC_FMTvSal$Comp <- as.factor("FMT vs. Saline")
res_DC_FMTvSal$Time <- "Discharge"
names(res_DC_FMTvSal) <- c("Taxa", "log2FC", "P-Value", "FDR", "Comparison", "Time")

res_DC_FMTvDonor <- res[, c(1, 9:11)]
res_DC_FMTvDonor <- res_DC_FMTvDonor[res_DC_FMTvDonor$rawP.dis_FMTvsDonor < 0.05, ]
res_DC_FMTvDonor$Comp <- as.factor("FMT vs. Donor")
res_DC_FMTvDonor$Time <- "Discharge"
names(res_DC_FMTvDonor) <- c("Taxa", "log2FC", "P-Value", "FDR", "Comparison", "Time")

res_DC_SalvDonor <- res[, c(1, 12:14)]
res_DC_SalvDonor <- res_DC_SalvDonor[res_DC_SalvDonor$rawP.dis_SalinevsDonor < 0.05, ]
res_DC_SalvDonor$Comp <- as.factor("Saline vs. Donor")
res_DC_SalvDonor$Time <- "Discharge"
names(res_DC_SalvDonor) <- c("Taxa", "log2FC", "P-Value", "FDR", "Comparison", "Time")

res_30_FMTvSal <- res[, c(1, 15:17)]
res_30_FMTvSal <-  res_30_FMTvSal[res_30_FMTvSal$rawP.d30_FMTvsSaline <0.05, ]
res_30_FMTvSal$Comp <-  as.factor("FMT vs. Saline")
res_30_FMTvSal$Time <- "30 Days"
names(res_30_FMTvSal) <- c("Taxa", "log2FC", "P-Value", "FDR", "Comparison", "Time")

res_30_FMTvDonor <- res[, c(1, 18:20)]
res_30_FMTvDonor <-  res_30_FMTvDonor[res_30_FMTvDonor$rawP.d30_FMTvsDonor <0.05, ]
res_30_FMTvDonor$Comp <-  as.factor("FMT vs. Donor")
res_30_FMTvDonor$Time <- "30 Days"
names(res_30_FMTvDonor) <- c("Taxa", "log2FC", "P-Value", "FDR", "Comparison", "Time")

res_30_salvDonor <- res[, c(1, 21:23)]
res_30_salvDonor <-  res_30_salvDonor[res_30_salvDonor$rawP.d30_SalinevsDonor <0.05, ]
res_30_salvDonor$Comp <-  as.factor("Saline vs. Donor")
res_30_salvDonor$Time <- "30 Days"
names(res_30_salvDonor) <- c("Taxa", "log2FC", "P-Value", "FDR", "Comparison", "Time")


df <- dplyr::bind_rows(res_admit, res_30_FMTvDonor, res_30_FMTvSal, 
                       res_30_salvDonor, res_DC_FMTvDonor, 
                       res_DC_FMTvSal, res_DC_SalvDonor)


df$Time <- as.factor(df$Time)

df$Comparison <- as.factor(df$Comparison)

df$Time <- factor(df$Time, levels = c("Admission", "Discharge", "30 Days"))

levels(df$Comparison)

df$log2FC <- ifelse(df$FDR > 0.2, NA, df$log2FC)

df$Taxa <- ifelse(is.na(df$log2FC), NA, df$Taxa)

gplot1 <- ggplot(data = df, aes(x= factor(Comparison, 
                                levels = c("AHDS vs. Donor", "FMT vs. Saline", "FMT vs. Donor", "Saline vs. Donor")),
                      y = Taxa)) + 
  geom_point(aes(size = abs(log2FC), color = ifelse(df$log2FC > 0, "Increased", "Decreased"))) +
  scale_color_manual(values=c(" blue", "red"), na.translate = F) +
  scale_size(breaks = c(2, 4, 6, 8, 10, 12, 14)) +
  scale_y_discrete(na.translate = FALSE) +
  labs(y="OTU (Genus_species)", 
       x="Statistical Contrast", 
       size="Log2-Fold Change",
       color="") +
  facet_wrap(df$Time, scales = "free_x") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))

png("heatmap.png", height = 30, width = 30, units = "cm", res = 600)
heatmap
dev.off()

##convert png to grob
img <- readPNG("heatmap.png")

grob <- rasterGrob(img)

tiff("Glimma.tiff", height = 30, width = 50, units = "cm", res = 600)
ggarrange(gplot1, grob, nrow = 1, ncol = 2, labels = c("A", "B"))
dev.off()

#Export results ====

write.csv(res, file = "FMT_limma_res.csv")

#additional plots ====

#build list of interesting OTUs

df$OTU <- rownames(res)[match(df$Taxa, res$taxa)]

unique(df$Taxa)

class(res$taxa)

res$taxa <- as.factor(res$taxa)

unique(res$taxa)

taxlist <- c("OTU_614", "OTU_1114", "OTU_1095", "OTU_1765", "OTU_1711", "OTU_1240", "OTU_2304", "OTU_925", "OTU_294", "OTU_1299")

rownames(res[c(36, 39, 44, 49, 70, 73, 77, 79, 68), ])

#extract counts for interesting OTUs, prepare data

otu <- as.data.frame(t(logCPM))
sam <- as.data.frame(sample_data(physeq_glom))

otu <- otu[, taxlist]

head(otu)

names <- taxa$gs[match(names(otu), rownames(taxa))]

names(otu) <- names

otu$Treatment_Group <- sam$Treatment_Group[match(rownames(otu), rownames(sam))]

otu$Collection_Time <- sam$Collection_Time[match(rownames(otu), rownames(sam))]

otu$Diagnosis <- sam$Diagnosis[match(rownames(otu), rownames(sam))]

#summary statistics function

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

names(otu)

#generate summary stats and plot mean abundance

dat_summary1 <- data_summary(data = otu, 
                            varname = "[Eubacterium]_biforme", #name of variable column
                            groupnames = c("Treatment_Group", "Collection_Time")) #names of grouping vars


p1 <- ggplot(dat_summary1, aes(x= Collection_Time, 
                        y = dat_summary1$`[Eubacterium]_biforme`, 
                        group = Treatment_Group,
                        color = Treatment_Group)) +
  geom_errorbar(aes(ymin=dat_summary1$`[Eubacterium]_biforme`-sd, ymax=dat_summary1$`[Eubacterium]_biforme`+sd), width=.1, position=position_dodge(0.05)) +
  geom_line() + 
  geom_point() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("") +
  labs(color = "Treatment Group") +
  ggtitle("Eubacterium biforme") +
  theme(plot.title = element_text(hjust = 0.5, face = "italic"))


dat_summary2 <- data_summary(data = otu, 
                            varname = "Butyricicoccus_pullicaecorum", #name of variable column
                            groupnames = c("Treatment_Group", "Collection_Time")) #names of grouping vars


p2 <- ggplot(dat_summary2, aes(x= Collection_Time, 
                        y = Butyricicoccus_pullicaecorum, 
                        group = Treatment_Group,
                        color = Treatment_Group)) +
  geom_errorbar(aes(ymin=Butyricicoccus_pullicaecorum-sd, ymax=Butyricicoccus_pullicaecorum+sd), width=.1, position=position_dodge(0.05)) +
  geom_line() + 
  geom_point() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("") +
  labs(color = "Treatment Group") +
  ggtitle("Butyricicoccus pullicaecorum") +
  theme(plot.title = element_text(hjust = 0.5, face = "italic"))


dat_summary3 <- data_summary(data = otu, 
                            varname = "Catenibacterium_NA", #name of variable column
                            groupnames = c("Treatment_Group", "Collection_Time")) #names of grouping vars


p3 <- ggplot(dat_summary3, aes(x= Collection_Time, 
                        y = Catenibacterium_NA, 
                        group = Treatment_Group,
                        color = Treatment_Group)) +
  geom_errorbar(aes(ymin=Catenibacterium_NA-sd, ymax=Catenibacterium_NA+sd), width=.1, position=position_dodge(0.05)) +
  geom_line() + 
  geom_point() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("") +
  labs(color = "Treatment Group") +
  ggtitle("Catenibacterium spp.") +
  theme(plot.title = element_text(hjust = 0.5, face = "italic"))



dat_summary4 <- data_summary(data = otu, 
                             varname = "Roseburia_inulinivorans", #name of variable column
                             groupnames = c("Treatment_Group", "Collection_Time")) #names of grouping vars


p4 <- ggplot(dat_summary4, aes(x= Collection_Time, 
                               y = Roseburia_inulinivorans, 
                               group = Treatment_Group,
                               color = Treatment_Group)) +
  geom_errorbar(aes(ymin=Roseburia_inulinivorans-sd, ymax=Roseburia_inulinivorans+sd), width=.1, position=position_dodge(0.05)) +
  geom_line() + 
  geom_point() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("") +
  labs(color = "Treatment Group") +
  ggtitle("Roseburia inulinivorans") +
  theme(plot.title = element_text(hjust = 0.5, face = "italic"))


dat_summary5 <- data_summary(data = otu, 
                             varname = "Prevotella_copri", #name of variable column
                             groupnames = c("Treatment_Group", "Collection_Time")) #names of grouping vars


p5 <- ggplot(dat_summary5, aes(x= Collection_Time, 
                               y = Prevotella_copri, 
                               group = Treatment_Group,
                               color = Treatment_Group)) +
  geom_errorbar(aes(ymin=Prevotella_copri-sd, ymax=Prevotella_copri+sd), width=.1, position=position_dodge(0.05)) +
  geom_line() + 
  geom_point() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("") +
  labs(color = "Treatment Group") +
  ggtitle("Prevotella copri") +
  theme(plot.title = element_text(hjust = 0.5, face = "italic"))

dat_summary6 <- data_summary(data = otu, 
                             varname = "Faecalibacterium_NA", #name of variable column
                             groupnames = c("Treatment_Group", "Collection_Time")) #names of grouping vars


p6 <- ggplot(dat_summary6, aes(x= Collection_Time, 
                               y = Faecalibacterium_NA, 
                               group = Treatment_Group,
                               color = Treatment_Group)) +
  geom_errorbar(aes(ymin=Faecalibacterium_NA-sd, ymax=Faecalibacterium_NA+sd), width=.1, position=position_dodge(0.05)) +
  geom_line() + 
  geom_point() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("") +
  labs(color = "Treatment Group") +
  ggtitle("Faecalibacterium prausnitzii") +
  theme(plot.title = element_text(hjust = 0.5, face = "italic"))


dat_summary7 <- data_summary(data = otu, 
                             varname = "Clostridium_hiranonis", #name of variable column
                             groupnames = c("Treatment_Group", "Collection_Time")) #names of grouping vars


p7 <- ggplot(dat_summary7, aes(x= Collection_Time, 
                               y = Clostridium_hiranonis, 
                               group = Treatment_Group,
                               color = Treatment_Group)) +
  geom_errorbar(aes(ymin=Clostridium_hiranonis-sd, ymax=Clostridium_hiranonis+sd), width=.1, position=position_dodge(0.05)) +
  geom_line() + 
  geom_point() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("") +
  labs(color = "Treatment Group") +
  ggtitle("Clostridium hiranonis") +
  theme(plot.title = element_text(hjust = 0.5, face = "italic"))

dat_summary8 <- data_summary(data = otu, 
                             varname = "Turicibacter_NA", #name of variable column
                             groupnames = c("Treatment_Group", "Collection_Time")) #names of grouping vars


p8 <- ggplot(dat_summary8, aes(x= Collection_Time, 
                               y = Turicibacter_NA, 
                               group = Treatment_Group,
                               color = Treatment_Group)) +
  geom_errorbar(aes(ymin=Turicibacter_NA-sd, ymax=Turicibacter_NA+sd), width=.1, position=position_dodge(0.05)) +
  geom_line() + 
  geom_point() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("") +
  labs(color = "Treatment Group") +
  ggtitle("Turicibacter spp.") +
  theme(plot.title = element_text(hjust = 0.5, face = "italic"))



multiplot <- ggarrange(p1, p2, p3, p4, p5, p6, p8, p7, 
          ncol = 4,
          nrow = 2,
          common.legend = TRUE, legend = "bottom",
          labels = LETTERS[1:10]
) 

tiff("diff_taxa.tiff", height = 15, width = 40, units = "cm", res = 600)
annotate_figure(multiplot, left = text_grob("Mean Abundance (logCPM)", color = "black", rot = 90))
dev.off()

#output data to make tables

tab <- dplyr::bind_rows(res_admit, res_30_FMTvDonor, res_30_FMTvSal, 
                       res_30_salvDonor, res_DC_FMTvDonor, 
                       res_DC_FMTvSal, res_DC_SalvDonor)


write_csv(tab, "res_Glimma_FMT.csv")






              