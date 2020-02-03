

library(phyloseq)
library(microbiome)
library(mixOmics)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(kableExtra)
library(igraph)

source("Functions_Microbiome.R")

#Load data ====

physeq <- readRDS("FMT_physeq_glom.RDS") 

physeq
summarize_phyloseq(physeq)

#do this to obtain genus and species for results
tax <- as(tax_table(physeq), "matrix")

tax <- as.data.frame(tax)

tax_species <- read.csv("FMT_taxa_table.csv", row.names = 1)

tax_species$g_s <- paste(tax_species$Genus, tax_species$Species, sep = "_")

tax$Species <- tax_species$g_s[match(rownames(tax), rownames(tax_species))]

tax <- as(tax, "matrix")

tax <- tax_table(tax)

tax_table(physeq) <- tax

physeq

#AHDS vs. Healthy controls ====

#Subset by Admission collection time
admit <- subset_samples(physeq, Collection_Time == "Admission")

admit <- prune_taxa(taxa_sums(admit) > 0, admit)

admit
summarize_phyloseq(admit)

#extract OTU table as matrix; need OTUs as columns
OTU_admit <- as.matrix(t((otu_table(admit))))

##need to add offset for CLR transformation
OTU_admit <- OTU_admit +1
sum(which(OTU_admit == 0))


#perform TSS normalization
OTU_admit.TSS = t(apply(OTU_admit, 1, TSS.divide))

#assign outcome variable
dx <- admit@sam_data$Diagnosis

#PCA
pca = pca(OTU_admit.TSS, ncomp = 10, logratio = 'CLR')

plot(pca)

plotIndiv(pca, 
          comp = c(1,2), # the components to plot
          pch = 16, 
          ind.names = F, 
          group = dx, 
          col.per.group = color.mixo(1:2),
          legend = TRUE,
          title = 'AHDS and Healthy Donors: PCA comp 1 - 2')

## Supervised Analysis with PLSDA

plsda = plsda(X = OTU_admit.TSS, dx, ncomp = nlevels(dx), logratio = 'CLR')

#assess model performance
set.seed(123) #for reproducibility

perf.plsda = perf(plsda, validation = 'Mfold', folds = 4,
                  progressBar = FALSE, nrepeat = 100)

plot(perf.plsda, overlay = 'measure', sd = TRUE)

plotIndiv(plsda, comp = c(1,2), ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'AHDS and Healthy Donors: PLSDA comp 1 - 2')

# Tune sPLSDA

set.seed(33)  # for reproducible results

tune.splsda = tune.splsda(OTU_admit.TSS, dx, ncomp = 2, 
                          logratio = 'CLR',
                          test.keepX = c(seq(5,150, 5)), validation = 'Mfold', 
                          folds = 4, dist = 'max.dist', nrepeat = 100,
                          progressBar = FALSE)


plot(tune.splsda)

select.keepX <- c(10, 15)

select.keepX

#Final sPLSDA model

splsda = splsda(OTU_admit.TSS, dx, ncomp = 2, logratio = 'CLR', keepX = select.keepX) 

plotIndiv(splsda, comp = c(1,2),
          ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE,
          title = 'AHDS vs. Healthy Donors: \nsPLSDA comp 1 - 2')

# Evaluate final sPLSDA model

set.seed(34)  # for reproducible results for this code

perf.splsda = perf(splsda, validation = 'Mfold', folds = 4, 
                   progressBar = FALSE, nrepeat = 100, dist = 'max.dist')

perf.splsda$error.rate

plot(perf.splsda)

#obtain results

selected.OTU.comp1 <- selectVar(splsda, comp = 1)$name

name.var <- tax_table(admit)[, "Species"]

plotLoadings(splsda, 
             comp = 1,
             method = 'mean', 
             contrib = 'max', 
             ndisplay = 10, 
             name.var = name.var, 
             size.title = 1, 
             size.name = 0.5, 
             size.legend = 0.5,
             title = "Donor vs AHDS at Admisstion: \nOTUs Selected on sPLSDA Comp1")

plotLoadings(splsda, 
             comp = 2,
             method = 'mean', 
             contrib = 'max', 
             ndisplay = 10, 
             name.var = name.var, 
             size.title = 1, 
             size.name = 0.5, 
             size.legend = 0.5,
             title = "Donor vs AHDS at Admission: \nOTUs Selected on sPLSDA Comp2")


#Analyve FMT group over time ====

## Subset to Obtain Samples Collected from Dogs Treated with FMT

FMT <- subset_samples(physeq, Treatment_Group == "FMT")

FMT <- prune_taxa(taxa_sums(FMT) > 0, FMT)

FMT
summarize_phyloseq(FMT)

#extract OTU table as matrix; need OTUs as columns
OTU_FMT <- as.matrix(t((otu_table(FMT))))

##need to add offset for CLR transformation
OTU_FMT <- OTU_FMT +1
sum(which(OTU_FMT == 0))


#perform TSS normalization
OTU_FMT.TSS = t(apply(OTU_FMT, 1, TSS.divide))

## Assign "Collection_Time" as Outcome Variable for sPLSDA
collection_time <- FMT@sam_data$Collection_Time

length(collection_time)
str(collection_time)


## Assign DogID for Repeated Measures:
DogID <- FMT@sam_data$DogID

length(DogID)
str(DogID)

## Unsupervised Analysis: PCA
pca2 = pca(OTU_FMT.TSS, ncomp = 10, logratio = 'CLR', multilevel = DogID)

plot(pca2)

plotIndiv(pca2, 
          comp = c(1,2), # the components to plot
          pch = 16, 
          ind.names = F, 
          group = collection_time, 
          col.per.group = color.mixo(1:3),
          legend = TRUE,
          title = 'FMT: PCA comp 1 - 2')

## Supervised Analysis with PLSDA
plsda2 = plsda(X = OTU_FMT.TSS, 
               collection_time, 
               ncomp = nlevels(collection_time), 
               logratio = 'CLR', 
               multilevel = DogID)

set.seed(123)
perf.plsda2 = perf(plsda2, validation = 'Mfold', folds = 4,
                   progressBar = FALSE, nrepeat = 100)

plot(perf.plsda2, overlay = 'measure', sd = TRUE)

plotIndiv(plsda2, comp = c(1,2), ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'FMT: PLSDA comp 1 - 2')

## Tune sPLSDA
set.seed(33)  # for reproducible results

tune.splsda2 = tune.splsda(OTU_FMT.TSS, collection_time, ncomp = 3, 
                           logratio = 'CLR', multilevel = DogID,
                           test.keepX = c(seq(5,150, 5)), validation = 'Mfold', 
                           folds = 4, dist = 'max.dist', nrepeat = 100,
                           progressBar = FALSE)


plot(tune.splsda2)

select.keepX2 <- c(10, 30, 10)

#final model
splsda2 = splsda(OTU_FMT.TSS, 
                 collection_time, 
                 ncomp = 3, 
                 logratio = 'CLR',  
                 multilevel = DogID, 
                 keepX = select.keepX2) 

plotIndiv(splsda2, comp = c(1,2),
          ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE,
          title = 'FMT: sPLSDA comp 1 - 2')

plotIndiv(splsda2, comp = c(1,3),
          ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE,
          title = 'FMT: sPLSDA comp 1 - 3')

## Evaluate sPLSDA
set.seed(34)  # for reproducible results for this code

perf.splsda2 = perf(splsda2, validation = 'Mfold', folds = 4, 
                    progressBar = FALSE, nrepeat = 100, dist = 'max.dist')

perf.splsda2$error.rate

plot(perf.splsda2)

#obtain results
selected.OTU.comp1_FMT = selectVar(splsda2, comp = 1)$name

name.var2 <- tax_table(FMT)[, "Species"]

splsda2$Y <-  factor(splsda2$Y, 
       levels= c("Admission", "Discharge", "30days"))

plotLoadings(splsda2, 
             comp = 1,
             method = 'mean', 
             contrib = 'max', 
             ndisplay = 10, 
             name.var = name.var2, 
             size.title = 1, 
             size.name = 0.5, 
             size.legend = 0.5,
             title = "FMT vs Time: \nOTUs Selected on sPLSDA Comp1")

plotLoadings(splsda2, 
             comp = 2,
             method = 'mean', 
             contrib = 'max', 
             ndisplay = 10, 
             name.var = name.var2, 
             size.title = 1, 
             size.name = 0.5, 
             size.legend = 0.5,
             title = "FMT vs Time: \nOTUs Selected on sPLSDA Comp2")


#Compare FMT vs Saline at Discharge ====

# Subset Phyloseq Object to Obtain Samples Collected at the Time of Discharge

DC <- subset_samples(physeq, Collection_Time == "Discharge")

DC <- subset_samples(DC, Treatment_Group == "FMT" | Treatment_Group == "Saline")

#extract OTU table as matrix; need OTUs as columns
OTU_DC <- as.matrix(t((otu_table(DC))))

##need to add offset for CLR transformation
OTU_DC <- OTU_DC +1
sum(which(OTU_DC == 0))


#perform TSS normalization
OTU_DC.TSS = t(apply(OTU_DC, 1, TSS.divide))

tx2 <- DC@sam_data$Treatment_Group

length(tx2)
str(tx2)

## Unsupervised Analysis: PCA
pca5 = pca(OTU_DC.TSS, ncomp = 8, logratio = 'CLR')

plot(pca5)

## Supervised Analysis with PLSDA

plsda5 = plsda(X = OTU_DC.TSS, tx2, ncomp = nlevels(tx2), logratio = 'CLR')

set.seed(123)
perf.plsda5 = perf(plsda5, validation = 'Mfold', folds = 4,
                   progressBar = FALSE, nrepeat = 100)

plot(perf.plsda5, overlay = 'measure', sd = TRUE)

plotIndiv(plsda5, comp = c(1,2), ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'FMT vs Saline at Discharge: PLSDA comp 1 - 2')

## Tune sPLSDA

set.seed(33)  # for reproducible results

tune.splsda5 = tune.splsda(OTU_DC.TSS, tx2, ncomp = 2, 
                           logratio = 'CLR',
                           test.keepX = c(seq(5,150, 5)), validation = 'Mfold', 
                           folds = 4, dist = 'max.dist', nrepeat = 100,
                           progressBar = FALSE)


plot(tune.splsda5)

select.keepX5 = tune.splsda5$choice.keepX[1:2]

## final sPLSDA model

splsda5 = splsda(OTU_DC.TSS, tx2, ncomp = 2, logratio = 'CLR', keepX = select.keepX5) 

plotIndiv(splsda5, comp = c(1,2),
          ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE,
          title = 'FMT vs Saline at Discharge: \nsPLSDA comp 1 - 2')

## Evaluate sPLSDA

set.seed(34)  # for reproducible results for this code

perf.splsda5 = perf(splsda5, validation = 'Mfold', folds = 4, 
                    progressBar = FALSE, nrepeat = 100, dist = 'max.dist')

perf.splsda5$error.rate

plot(perf.splsda5)

#obtain results
name.var5 <- tax_table(DC)[, "Species"]

plotLoadings(splsda5, 
             comp = 1,
             method = 'mean', 
             contrib = 'max', 
             ndisplay = 10, 
             name.var = name.var5, 
             size.title = 1, 
             size.name = 0.5, 
             size.legend = 0.5, 
             title = "FMT vs Saline at Discharge: \nOTUs Selected on sPLSDA Comp1")

plotLoadings(splsda5, 
             comp = 2,
             method = 'mean', 
             contrib = 'max', 
             ndisplay = 10, 
             name.var = name.var5, 
             size.title = 1, 
             size.name = 0.5, 
             size.legend = 0.5, 
             title = "FMT vs Saline at Discharge: \nOTUs Selected on sPLSDA Comp2")

#Compare FMT vs Saline at 30 days

## Subset Phyloseq Object to Obtain Samples Collected at the Time of Discharge

d30 <- subset_samples(physeq, Collection_Time == "30days")

d30 <- subset_samples(d30, Treatment_Group == "FMT" | Treatment_Group == "Saline")

#extract OTU table as matrix; need OTUs as columns
OTU_d30 <- as.matrix(t((otu_table(d30))))

##need to add offset for CLR transformation
OTU_d30 <- OTU_d30 +1
sum(which(OTU_d30 == 0))


#perform TSS normalization
OTU_d30.TSS = t(apply(OTU_d30, 1, TSS.divide))

tx3 <- d30@sam_data$Treatment_Group

length(tx3)
str(tx3)

## Unsupervised Analysis: PCA
pca6 = pca(OTU_d30.TSS, ncomp = 8, logratio = 'CLR')

plot(pca6)

plotIndiv(pca6, 
          comp = c(1,2), # the components to plot
          pch = 16, 
          ind.names = F, 
          group = tx3, 
          col.per.group = color.mixo(1:2),
          legend = TRUE,
          title = 'FMT vs Saline at 30 Days: PCA comp 1 - 2')

## Supervised Analysis with PLSDA

plsda6 = plsda(X = OTU_d30.TSS, tx3, ncomp = nlevels(tx3), logratio = 'CLR')

set.seed(123)
perf.plsda6 = perf(plsda6, validation = 'Mfold', folds = 4,
                   progressBar = FALSE, nrepeat = 100)

plot(perf.plsda6, overlay = 'measure', sd = TRUE)

plotIndiv(plsda6, comp = c(1,2), ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'FMT vs Saline at 30 Days: PLSDA comp 1 - 2')



## Tune sPLSDA

set.seed(33)  # for reproducible results

tune.splsda6 = tune.splsda(OTU_d30.TSS, tx3, ncomp = 2, 
                           logratio = 'CLR',
                           test.keepX = c(seq(5,150, 5)), validation = 'Mfold', 
                           folds = 4, dist = 'max.dist', nrepeat = 10,
                           progressBar = FALSE)


plot(tune.splsda6)

select.keepX6 = tune.splsda$choice.keepX[1:2]

## Final sPLSDA model
splsda6 = splsda(OTU_d30.TSS, tx3, ncomp = 2, logratio = 'CLR', keepX = select.keepX6) 

plotIndiv(splsda6, comp = c(1,2),
          ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE,
          title = 'FMT vs Saline at 30 Days: sPLSDA comp 1 - 2')

## Evaluate sPLSDA
set.seed(34)  # for reproducible results for this code

perf.splsda6 = perf(splsda6, validation = 'Mfold', folds = 4, 
                    progressBar = FALSE, nrepeat = 100, dist = 'max.dist')

perf.splsda6$error.rate

plot(perf.splsda6)

#obtain results

name.var6 <- tax_table(d30)[, "Species"]

plotLoadings(splsda6, 
             comp = 1,
             method = 'mean', 
             contrib = 'max', 
             ndisplay = 10, 
             name.var = name.var6, 
             size.title = 1, 
             size.name = 0.5, 
             size.legend = 0.5, 
             title = "FMT vs Saline at 30 Days: \nOTUs Selected on sPLSDA Comp1")

plotLoadings(splsda6, 
             comp = 2,
             method = 'mean', 
             contrib = 'max', 
             ndisplay = 10, 
             name.var = name.var6, 
             size.title = 1, 
             size.name = 0.5, 
             size.legend = 0.5, 
             title = "FMT vs Saline at 30 Days: \nOTUs Selected on sPLSDA Comp2")





