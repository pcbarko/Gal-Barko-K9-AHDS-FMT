
library(dplyr)

res <- read.csv("Sig_OTUs_2020-01-31.csv")

taxa <- read.csv("FMT_taxa_table.csv")

taxa$gs <- paste(taxa$Genus, taxa$Species, sep = "_")

res$taxa <- taxa$gs[match(res$OTU, taxa$X)]

res <- res[, c(32, 1:31)]


res_admit <- res[, c(1, 5:7)]
res_admit <- res_admit[res_admit$rawP.ad_AHDSvsDonor < 0.05, ]
res_admit$Comparison <- "AHDS vs. Donor"
res_admit$Time <- "Admission"
names(res_admit) <- c("Taxa", "log2FC", "P-Value", "FDR", "Comparison", "Time")


res_DC_FMTvSal <- res[, c(1, 8:10)]
res_DC_FMTvSal <- res_DC_FMTvSal[res_DC_FMTvSal$rawP.dis_FMTvsSaline < 0.05, ]
res_DC_FMTvSal$Comp <- "FMT vs. Saline"
res_DC_FMTvSal$Time <- "Discharge"
names(res_DC_FMTvSal) <- c("Taxa", "log2FC", "P-Value", "FDR", "Comparison", "Time")

res_DC_FMTvDonor <- res[, c(1, 11:13)]
res_DC_FMTvDonor <- res_DC_FMTvDonor[res_DC_FMTvDonor$rawP.dis_FMTvsDonor < 0.05, ]
res_DC_FMTvDonor$Comp <- "FMT vs. Donor"
res_DC_FMTvDonor$Time <- "Discharge"
names(res_DC_FMTvDonor) <- c("Taxa", "log2FC", "P-Value", "FDR", "Comparison", "Time")

res_DC_SalvDonor <- res[, c(1, 14:16)]
res_DC_SalvDonor <- res_DC_SalvDonor[res_DC_SalvDonor$rawP.dis_SalinevsDonor < 0.05, ]
res_DC_SalvDonor$Comp <- "Saline vs. Donor"
res_DC_SalvDonor$Time <- "Discharge"
names(res_DC_SalvDonor) <- c("Taxa", "log2FC", "P-Value", "FDR", "Comparison", "Time")

res_30_FMTvSal <- res[, c(1, 17:19)]
res_30_FMTvSal <-  res_30_FMTvSal[res_30_FMTvSal$rawP.d30_FMTvsSaline <0.05, ]
res_30_FMTvSal$Comp <-  "FMT vs. Saline"
res_30_FMTvSal$Time <- "30 Days"
names(res_30_FMTvSal) <- c("Taxa", "log2FC", "P-Value", "FDR", "Comparison", "Time")

res_30_FMTvDonor <- res[, c(1, 20:22)]
res_30_FMTvDonor <-  res_30_FMTvDonor[res_30_FMTvDonor$rawP.d30_FMTvsDonor <0.05, ]
res_30_FMTvDonor$Comp <-  "FMT vs. Donor"
res_30_FMTvDonor$Time <- "30 Days"
names(res_30_FMTvDonor) <- c("Taxa", "log2FC", "P-Value", "FDR", "Comparison", "Time")

res_30_salvDonor <- res[, c(1, 23:25)]
res_30_salvDonor <-  res_30_salvDonor[res_30_salvDonor$rawP.d30_SalinevsDonor <0.05, ]
res_30_salvDonor$Comp <-  "Saline vs. Donor"
res_30_salvDonor$Time <- "30 Days"
names(res_30_salvDonor) <- c("Taxa", "log2FC", "P-Value", "FDR", "Comparison", "Time")


list <- c(res_30_FMTvDonor$Taxa, res_30_FMTvSal$Taxa, 
             res_30_salvDonor$Taxa, res_DC_FMTvDonor$Taxa, 
             res_DC_FMTvSal$Taxa, res_DC_SalvDonor$Taxa)

unique(list)


df <- dplyr::bind_rows(res_admit, res_30_FMTvDonor, res_30_FMTvSal, 
                       res_30_salvDonor, res_DC_FMTvDonor, 
                       res_DC_FMTvSal, res_DC_SalvDonor)

df$Time <- as.factor(df$Time)

df$Comparison <- as.factor(df$Comparison)

df$Time <- factor(df$Time, levels = c("Admission", "Discharge", "30 Days"))

levels(df$Comparison)

#df$Comparison <- factor(df$Comparison, levels = c("FMT vs. Saline", "FMT vs. Donor", "Saline vs. Donor"))

#df <- df[!df$FDR > 0.1, ]

p1 <- ggplot(data = df, aes(x= factor(Comparison, 
                                      levels = c("AHDS vs. Donor", "FMT vs. Saline", "FMT vs. Donor", "Saline vs. Donor")
                                      ), y = Taxa)) + 
  geom_point(aes(size = abs(log2FC), color = ifelse(df$log2FC > 0, "Increased", "Decreased"))) +
  scale_color_manual(values=c(" blue", "red")) +
  scale_size(breaks = c(50, 100, 500, 1000, 2000, 5000)) +
  labs(y="OTU (Genus_Species)", 
       x="Contrast", 
       size="Log2-Fold Change",
       color="Direction of Difference") +
  facet_wrap(df$Time, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5)) 

p1 

names(res_admit) <- c("Taxa", "Log2-Fold Change", "P-Value", "FDR")

p2 <- ggplot(data = res_admit, aes(x= res_admit$Taxa, y = res_admit$`Log2-Fold Change`)) + 
  geom_point(aes(size = abs(res_admit$`Log2-Fold Change`), 
                 color = ifelse(res_admit$`Log2-Fold Change`> 0, "Increased", "Decreased"))) +
  scale_color_manual(values=c("blue", "red")) +
  scale_size(breaks = c(100, 500, 1000, 5000, 8000)) +
  labs(y="Log2-Fold Change", 
       x="OTU (Genus_species)", 
       size="Log2-Fold Change",
       color="Direction of Difference") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))

p2

library(ggpubr)

ggarrange(p1, p2, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1, 
          common.legend = T)

ggarrange(p2, p1, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1, 
          common.legend = T)








