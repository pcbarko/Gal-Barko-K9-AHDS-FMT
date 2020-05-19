library(mixOmics)
library(WGCNA)
library(ggplot2)
library(microbiome)
library(kableExtra)
library(lme4)
library(lmerTest)
library(emmeans)
library(png)
library(grid)
library(gridExtra)
library(ggpubr)
library(pheatmap)

#be sure to set working directory

#source("Functions_Microbiome.R")

#Load Data ====

physeq <- readRDS("FMT_physeq_glom.RDS") 

summarize_phyloseq(physeq)

#extract OTU table as matrix; need OTUs as columns
OTU <- data.frame(t((otu_table(physeq))))

#get full taxa table for species annotation
taxa <- read.csv("FMT_taxa_table.csv", row.names = 1)
taxa <- taxa[names(OTU), ]

#sample metadata
sample <- as.data.frame(as(sample_data(physeq), "matrix"))

sample <- sample[, -c(4, 6, 7)]

sample$Collection_Time  <- factor(sample$Collection_Time, 
                                  levels(sample$Collection_Time)[c(2, 3, 1)])

table(rownames(sample) == rownames(OTU))

table(rownames(taxa) == colnames(OTU))

#Normalization ====

##need to add a pseudocount offset for CLR transformation
OTU <- OTU +1
sum(which(OTU == 0))

##perform TSS normalization

TSS.divide = function(x){
  x/sum(x)
}

OTU_TSS <- t(apply(OTU, 1, TSS.divide))

OTU_final <- logratio.transfo(OTU_TSS, "CLR")

#Preprocessing for WGCNA ====

##sample network based on euclidian distance

A <- adjacency(t(OTU_final), type = "distance")

##whole network connectivity
k <- as.numeric(apply(A, 2, sum)) -1

##standardized connectivity
Z.k = scale(k)

##designate samples as outliers if Z.k value is below threshold
thresholdZ.k <- -5
outlierColor <- ifelse(Z.k < thresholdZ.k, "red", "black")

##generate cluster tree

sampleTree <- hclust(as.dist(1-A), method = "average")

datColors <- data.frame(outlierC = outlierColor)

##plot the samples dendrogram and the colors underneath

plotDendroAndColors(sampleTree, groupLabels = names(datColors),
                    colors = datColors)

#Choose a set of soft thresholding powers  ====
powers = c(1:20)  

##choose power based on SFT criterion (for unsigned network)
sft <- pickSoftThreshold(OTU_final, powerVector = powers)

##Plot the results: NOTE I was unable to get get the exact code from the tutorial working if I ran them line by line, but they will work if knitted to word doc

par(mfrow = c(1, 2)) 

##SFT index as a function of different powers 
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
     xlab = "soft threshold power", 
     ylab = "SFT, signed R^2",
     type = "n", 
     main = paste("Scale Independence"))

text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
     labels = powers, 
     col = "red") 

abline(h = 0.8, col = "red") 

##Mean connectivity as a function of different powers 
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], 
     xlab = "Soft Threshold (power)", 
     ylab = "Mean Connectivity", 
     type = "n", 
     main = paste("Mean Connectivity")
)

text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")


#Network Construction ====

cor <- WGCNA::cor
mergingThresh = 0.2

net <-  blockwiseModules(OTU_final, 
                       corType = "bicor", 
                       maxBlockSize = 5000, 
                       networkType = "signed hybrid", 
                       power = 4, 
                       minModuleSize = 5, 
                       mergeCutHeight = mergingThresh, 
                       deepSplit = 2, 
                       numericLabels = TRUE, 
                       saveTOMs = TRUE, 
                       pamRespectsDendro = FALSE, 
                       saveTOMFileBase = "FMT_TOM") 

moduleLabelsAutomatic <-  net$colors

##Convert labels to colors for plotting 
moduleColorsAutomatic <-  labels2colors(moduleLabelsAutomatic) 

##A data frame with module eigengenes can be obtained as follows 
MEsAutomatic <-  net$MEs

##Choose a module assignment 
moduleColors <-  moduleColorsAutomatic 

##Define numbers of genes and samples 
nOTU <-  ncol(OTU_final) 
nSamples <-  nrow(OTU_final) 

##Recalculate MEs with color labels 

MEs0 <-  moduleEigengenes(OTU_final, moduleColors)$eigengenes

MEs <-  orderMEs(MEs0)

blocknumber <- 1 

datColors <-  data.frame(moduleColorsAutomatic)[net$blockGenes[[blocknumber]],]

##calculate module membership (module eigengene based connectivity, kME):

datKME <- signedKME(OTU_final, MEs)

##Plot the dendrogram and the module colors underneath

plotDendroAndColors(net$dendrograms[[blocknumber]], 
                    colors = datColors, 
                    groupLabels = c("OTU \nCo-Abundance \nModules"),  
                    dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, 
                    guideHang = 0.05)

#Comparing Module Eigenvalues Among Treatment Groups and Collection Times ====

##extract module eigenvalues (MEs) and create df with sample data and MEs

MEs <- cbind(MEs, sample)

##plot module eigenvalues ====

p1 <- ggplot(MEs, aes(y = MEs$MEturquoise, 
                x = Collection_Time, 
                fill = Treatment_Group)) + 
  geom_boxplot() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("") +
  ylab("Module Eigenvalue") +
  ggtitle("Turquoise Module")

p2 <- ggplot(MEs, aes(y = MEs$MEblue, 
                x = Collection_Time, 
                fill = Treatment_Group)) + 
  geom_boxplot() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("") +
  ylab("Module Eigenvalue") +
  ggtitle("Blue Module")

p3 <- ggplot(MEs, aes(y = MEs$MEyellow, 
                x = Collection_Time, 
                fill = Treatment_Group)) + 
  geom_boxplot() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("") +
  ylab("Module Eigenvalue") +
  ggtitle("Yellow Module")

p4 <- ggplot(MEs, aes(y = MEs$MEbrown, 
                x = Collection_Time, 
                fill = Treatment_Group)) + 
  geom_boxplot() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("") +
  ylab("Module Eigenvalue") +
  ggtitle("Brown Module")

ggarrange(p1, p2, p3, p4, ncol =2, nrow =2, 
          common.legend = T, 
          labels = c("A", "B", "C"))

#Identify taxa contained in each module ====

table(colnames(OTU_final) == rownames(taxa))

module_taxa_df <- data.frame(net$colors, moduleColors)

module_taxa_df$genus <- taxa$Genus[match(rownames(taxa), rownames(module_taxa_df))]
module_taxa_df$species <- as.character(taxa$Species[match(rownames(taxa), rownames(module_taxa_df))])
module_taxa_df$species[is.na(module_taxa_df$species)] <- "spp."

module_taxa_df$Taxa <- paste(module_taxa_df$genus, module_taxa_df$species, sep = " ")

module_taxa_df <- module_taxa_df[, -c(1,3,4)]

##remove gray module
module_taxa_df <- module_taxa_df[ !module_taxa_df$moduleColors == "grey" , ]

module_taxa_df <- module_taxa_df[order(module_taxa_df$moduleColors), ]


#Statistical Analysis of Module Abundance by Treatment Group and Collection Time ====

###The brown module appears to increase in abundnace in FMT recipients
###The brown module also contains taxa previously identified as increased after FMT
##Test for significant changes in brown module abundance

##first, AHDS vs healthy at admission (Pre-FMT)

t.test(MEs$MEbrown ~ MEs$Diagnosis)
#significant P<0.0001


model_brown <- lmer(formula = MEbrown ~ Treatment_Group*Collection_Time + (1|DogID),
                    data = MEs)

summary(model_brown)

brown_posthoc <- emmeans(model_brown, pairwise ~ Collection_Time | Treatment_Group)

brown_posthoc$contrasts
###significant differences in FMT recipients at discharge and 30 days post FMT

brown_posthoc2 <- emmeans(model_brown, pairwise ~ Treatment_Group | Collection_Time )

brown_posthoc2$contrasts
###significant differences in FMT vs saline-recipients at discharge and 30 days


##The turquoise module seems to decrease after FMT and contains many potential pathobionts

t.test(MEs$MEturquoise ~ MEs$Diagnosis)
#not significant


model_turq <- lmer(formula = MEturquoise ~ Treatment_Group*Collection_Time + (1|DogID),
                    data = MEs)

summary(model_turq)

turq_posthoc <- emmeans(model_turq, pairwise ~ Collection_Time | Treatment_Group)

turq_posthoc$contrasts
###significant differences in FMT recipients at discharge and 30 days post FMT

turq_posthoc2 <- emmeans(model_turq, pairwise ~ Treatment_Group | Collection_Time )

turq_posthoc2$contrasts


#Graphics Output ====

##output module dendrogram as png
png("WGCNA.png", height = 20, width = 30, units = "cm", res = 600)
plotDendroAndColors(net$dendrograms[[blocknumber]], 
                    colors = datColors, 
                    groupLabels = c("OTU \nCo-Abundance \nModules"),  
                    dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, 
                    guideHang = 0.05)
dev.off()

##convert png to grob
img <- readPNG("WGCNA.png")

grob <- rasterGrob(img)

##stats to annotate boxplot
brown_stats <- as.data.frame(summary(brown_posthoc))

contrasts <- strsplit(as.character(brown_stats$contrasts.contrast), split = "-", fixed = T)

group1 <- trimws(sapply(contrasts, `[`, 1))
group2 <- trimws(sapply(contrasts, `[`, 2))

Treatment_Group <- brown_stats$emmeans.Treatment_Group
pval <- round(brown_stats$contrasts.p.value, 3)
y.position <- c(0.3, 0.35, 0.4, 0.3, 0.35, 0.4, 0.3, 0.35, 0.4)

emstats <- data.frame(Treatment_Group, group1, group2, pval, y.position)
emstats$signif <- ifelse(emstats$pval < 0.05, "*", "ns")


brown_stats2 <- as.data.frame(summary(brown_posthoc2, adjust = "fdr"))
contrasts2 <- strsplit(as.character(brown_stats2$contrasts.contrast), split = "-", fixed = T)

group1 <- trimws(sapply(contrasts2, `[`, 1))
group2 <- trimws(sapply(contrasts2, `[`, 2))

Collection_Time <- brown_stats2$emmeans.Collection_Time
pval2 <- round(brown_stats2$contrasts.p.value, 3)
y.position <- c(0.3, 0.35, 0.4, 0.3, 0.35, 0.4, 0.3, 0.35, 0.4)

emstats2 <- data.frame(Collection_Time, group1, group2, pval2, y.position)
emstats2$signif <- ifelse(emstats2$pval2 < 0.05, "*", "ns")

##boxplots of brown module with clean formmatting and annotation
gplot_brown <- ggboxplot(MEs, x = "Collection_Time", y = "MEbrown", 
                         fill = "Collection_Time", 
                         add = "jitter",
                         facet.by = "Treatment_Group") +
  labs(fill = "Collection Time", x= "", y = "Brown Module Eigenvalues") +
  theme(legend.position="bottom") + 
  stat_pvalue_manual(emstats, label = "signif")

gplot_brown2 <- ggboxplot(MEs, x = "Treatment_Group", y = "MEbrown", 
          fill = "Treatment_Group", 
          add = "jitter", 
          facet.by = "Collection_Time", 
          palette = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  labs(fill = "Treatment Group", x= "", y = "Brown Module Eigenvalues") +
  theme(legend.position="bottom") + 
  stat_pvalue_manual(emstats2, label = "signif")

##brown module table

brown <- module_taxa_df[module_taxa_df$moduleColors == "brown", ]

brown <- brown[, c(2, 4)]

brown$kMEBrown <- round(brown$kMEBrown, 2)

brown <- brown[order(brown$kMEBrown, decreasing = T), ]

brown <- brown[brown$kMEBrown > 0.5, ]

names(brown) <- c("Taxa", "kME Brown")

brown[4, 1] <- "Eubacterium biforme"
brown[8, 1] <- "Prevotella spp."


tbl <- ggtexttable(brown, rows = NULL, 
                   theme = ttheme(
  tbody.style = tbody_style(color = "black", 
                            face = "italic", 
                            size = 20), 
  colnames.style = colnames_style(face = "bold", 
                                  size = 24)))

##final output

##plots
tiff("WGCNA2.tiff", height = 30, width = 50, units = "cm", res = 600)
ggarrange(grob, tbl, gplot_brown, gplot_brown2, 
          ncol = 2, nrow = 2, 
          labels = c("A", "B", "C", "D"))
dev.off()

tiff("WGCNA2.tiff", height = 60, width = 50, units = "cm", res = 600)
ggarrange(grob, tbl, gplot_brown, gplot_brown2, 
          ncol = 2, nrow = 2, 
          labels = c("A", "B", "C", "D"))
dev.off()

##results table

module_taxa_df$kMEBlue <- datKME$kMEblue[match(rownames(module_taxa_df), rownames(datKME))]

module_taxa_df$kMEBrown <- datKME$kMEbrown[match(rownames(module_taxa_df), rownames(datKME))]

module_taxa_df$kMETurquoise <- datKME$kMEturquoise[match(rownames(module_taxa_df), rownames(datKME))]

module_taxa_df$kMEYellow <- datKME$kMEyellow[match(rownames(module_taxa_df), rownames(datKME))]

write.csv(module_taxa_df, file = "FMT_Network_Modules.csv")

#Supplemental results ====

OTU_final <- as.data.frame(as(OTU_final, "matrix"))

plotNetworkHeatmap(
  OTU_final, 
  names(OTU_final), 
  weights = NULL,
  useTOM = TRUE, 
  power = 4, 
  networkType = "signed", 
  main = "Heatmap of the Co-Abundance Network")

METree <- hclust( as.dist( 1 - cor( net$MEs )), method="average")

plot(METree,main="Clustering of modules merged at 0.2",xlab="",sub="")

abline(h=0.2, col= "red")

##Heatmap of eigenOTUs

MEs$group <- paste(MEs$Treatment_Group, MEs$Collection_Time, MEs$DogID, sep = "_")

whichmodule = "brown" 
datModule = OTU_final[, moduleColorsAutomatic == whichmodule] 
datModule$treatment <- MEs$Treatment_Group[match(rownames(datModule), rownames(MEs))]
datModule$colletion <- MEs$Collection_Time[match(rownames(datModule), rownames(MEs))]

treatOrder <- levels(datModule$treatment)
collectOrder  <-  levels(datModule$colletion)

datModule <- datModule[order(match(datModule$treatment, treatOrder), 
                match(datModule$colletion, collectOrder)), ]


datModule <- datModule[order(match(datModule$colletion, collectOrder),
                             match(datModule$treatment, treatOrder)), ]

rownames(datModule) <- datModule$group

pheatmap(t(datModule[, 1:15]),
         cluster_cols = T, 
         cutree_cols = 2, 
         main = "Heatmap of Brown Module EigenOTUs")


