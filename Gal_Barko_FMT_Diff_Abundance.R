
#Initialize session ====
library(DEFormats)
library(edgeR)
library(Glimma)
library(gplots)
library(statmod)
library(dplyr)
library(scales)
library(ggplot2)


source("summarizeFit.R")

#load data ====
physeq_glom <- readRDS("FMT_physeq_glom.RDS") 

taxa <- read.csv("FMT_taxa_table.csv")
taxa$gs <- paste(taxa$Genus, taxa$Species, sep = "_")

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

glMDSPlot(logCPM, top = nrow(logCPM), labels = rownames(dge$samples),
          groups = dge$samples[,c("group","DogID","Collection_Time","Donor_ID","Treatment_Group","Group")], folder = "Interactive_plots",
          html = "MDSclustering_logCPM+TMM")

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
codes.res <- decideTests(fit2)

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

rownames(heat.h1) <- taxa$gs[match(rownames(heat.h1), taxa$X)]

colnames(heat.h1) <- paste(dge$samples$Group, dge$samples$DogID, sep = "_")

cluster.x.h1 <- hclust(dist(heat.h1))
cluster.y.h1 <- hclust(dist(t(heat.h1)))

temp <- factor(dge$samples$Treatment_Group, levels = c("Saline","FMT","Donor"))
ind <- order(dge$samples$Collection_Time, temp)

heatmap.2(heat.h1[,ind], col = col.pan, Rowv=as.dendrogram(cluster.x.h1),
          Colv=FALSE,
          scale = "none", labRow = rownames(heat.h1), trace = "none", key = T, keysize = 1,
          density.info = "none", margins = c(10,15), symbreaks= FALSE, key.par = list(mar = c(4,1,4,4)),
          key.xlab = "SD from mean (Z-scale)",
          colsep = c(4,8,12,16,20,24,28,32), sepcolor = "black",
          ColSideColors = rep(c("black","green", "purple"), each = 12),
          cexRow = 1)

par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       legend = c("Admission", "Discharge", "30 Days"), # category labels
       col = c("black", "green", "purple"),  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)


res.out <- summarizeFit(fit2)

#plot results ====
res.out$taxa <- taxa$gs[match(rownames(res.out), taxa$X)]

res <- res.out[, c(23, 1:22)]


res_admit <- res[, c(1, 3:5)]
res_admit <- res_admit[res_admit$rawP.ad_AHDSvsDonor < 0.05, ]
res_admit$Comparison <- as.factor("AHDS vs. Donor")
res_admit$Time <- "Admission"
names(res_admit) <- c("Taxa", "log2FC", "P-Value", "FDR", "Comparison", "Time")


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

dev.off()
ggplot(data = df, aes(x= factor(Comparison, 
                                levels = c("AHDS vs. Donor", "FMT vs. Saline", "FMT vs. Donor", "Saline vs. Donor")),
                      y = Taxa)) + 
  geom_point(aes(size = abs(log2FC), color = ifelse(df$log2FC > 0, "Increased", "Decreased"))) +
  scale_color_manual(values=c(" blue", "red"), na.translate = F) +
  scale_size(breaks = c(50, 100, 500, 1000, 2000, 5000)) +
  scale_y_discrete(na.translate = FALSE) +
  labs(y="OTU (Genus_species)", 
       x="Contrast", 
       size="Log2-Fold Change",
       color="") +
  facet_wrap(df$Time, scales = "free_x") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))

#Export results ====

write.csv(res, file = "FMT_limma_res.csv")
