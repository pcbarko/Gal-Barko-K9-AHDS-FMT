
##Initialize Session====

library(phyloseq)
library(tidyverse)
library(ggplot2)
library(ape)
library(dplyr)
library(vegan)
library(microbiome)
library(ranacapa)
library(ggpubr)

#load OTU table
otu <- read.csv("FMT_OTU_table.csv", check.names = F, skipNul = T)

otu[2:37] <- lapply(otu[2:37], as.numeric)

#load taxa table
tax <- read.csv("FMT_taxa_table.csv", check.names = F, skipNul = T)


#confirm that taxa IDs are identical in otu and taxa tables
x <- (tax[[1]] == otu[[1]])
length(which(x == TRUE))
length(which(x == FALSE))

#load metadata table
meta <- read.csv("FMT_sample_data.csv",  check.names = F, skipNul = T)

meta <- meta[, -1]

meta$DogID <- as.factor(meta$DogID)
meta$`Donor_ID` <- as.factor(meta$`Donor_ID`)


#confirm that sample IDs are identical in otu and meta tables
colnames(otu[2:37]) == meta[[1]]

#convert tables to matrices

otu_mat <- as.matrix(otu[, -1])
class(otu_mat[[16]])
rownames(otu_mat) <- otu[, 1]


tax_mat <- as.matrix(tax[, -1])
rownames(tax_mat) <- tax[, 1]
colnames(tax_mat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

meta2 <- meta[, -1]
rownames(meta2) <- meta[, 1]

meta2$Collection_Time <-  as.factor(sub(" ", "", meta2$Collection_Time))


OTU <- otu_table(otu_mat, taxa_are_rows = T)
TAX <- tax_table(tax_mat)
META <- sample_data(meta2)

physeq <- phyloseq(OTU, TAX, META)

#filter out OTUs with 0 reads
physeq1 <- prune_taxa(taxa_sums(physeq) > 0, physeq)

#add phylogenetic tree
random_tree = rtree(ntaxa(physeq1), rooted=TRUE, tip.label=taxa_names(physeq1))
physeq1 = merge_phyloseq(physeq, random_tree)

#summary of phyloseq object
summarize_phyloseq(physeq1)

#Are there any mitochondrial or diet-related OTUs?
table(tax_table(physeq1)[,'Domain'])

#Sequencing depth ====

SeqDepth = colSums(otu_table(physeq1))
sample_data(physeq1)$SeqDepth <- SeqDepth

summary(SeqDepth)

#are sequence counts normally distributed?
qplot(SeqDepth, geom = "histogram")

shapiro.test(SeqDepth)

##Looks like sequence counts are normally distributed

#assess rarefaction plot to determine sequancing depth is sufficient to estimate alpha diversity
p_rare <- suppressMessages(ggrare(physeq1, step = 1000, 
                                  color = "Group", 
                                  label = "DogID", 
                                  se = FALSE,
                                  plot = FALSE
))

p_rare <- p_rare + ggtitle("Alpha Rarefaction")

##looks like sequencing depth is adequate
##There are a few samples with lower counts, but this does not seem to impact estimation of alpha diversity

#what is the range of reads per OTU?
range(taxa_sums(physeq1))

#distribution of reads per OTU
hist(log2(taxa_sums(physeq1)), 1000)

#distribution of reads per sample
ggplot(data = data.frame(
  SampleSums = sample_sums(physeq1),
  Names = factor(sample_names(physeq1), 
                 ordered = TRUE, 
                 levels = sample_names(physeq1)),
  Group = factor(sample_data(physeq1)$Treatment_Group, 
                 ordered = TRUE)), 
  aes(y = SampleSums, x = Names, fill = Group)) + 
  geom_bar(stat = 'identity' ) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#How many OTUs are unassigned at Phlyum level
table(tax_table(physeq1)[,"Phylum"], exclude = NULL)

#remove OTUs with unknown phlyum
physeq2 <- subset_taxa(physeq1, !is.na(Phylum))

table(tax_table(physeq2)[,"Phylum"], exclude = NULL)

#Perform simple prevalence filtering with a threshold of 5% ==== 

prevdf <- apply(otu_table(physeq2),
                MARGIN = ifelse(taxa_are_rows(physeq2), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)}
)

prevdf <- data.frame(Prevalence =  prevdf, 
                     TotalAbundance = taxa_sums(physeq2),
                     tax_table(physeq2))

pthresh <- 0.05
prevThreshold <- pthresh * nsamples(physeq2)

keepTaxa <- rownames(prevdf)[(prevdf$Prevalence >= prevThreshold)]
physeq3 <- prune_taxa(keepTaxa, physeq2)

table(tax_table(physeq3)[,"Phylum"], exclude = NULL)

#How many OTUs were retained after filtering?
physeq3
summarize_phyloseq(physeq3)

#Agglomerate OTUs on Genus:

physeq_glom <- tax_glom(physeq3, taxrank = "Genus", NArm = TRUE)

#How many OTUs are retained in the filtered, agglomerated phyloseq object?
physeq_glom
summarize_phyloseq(physeq_glom)

#Export data ====

#raw, unfiltered data
levels(physeq1@sam_data$Diagnosis)
physeq1@sam_data$Diagnosis  <- factor(physeq1@sam_data$Diagnosis, 
                                      levels(physeq1@sam_data$Diagnosis)[c(2, 1)])

levels(physeq1@sam_data$Collection_Time)
physeq1@sam_data$Collection_Time  <- factor(physeq1@sam_data$Collection_Time, 
                                            levels(physeq1@sam_data$Collection_Time)[c(2, 3, 1)])

levels(physeq1@sam_data$Treatment_Group)
physeq1@sam_data$Treatment_Group  <- factor(physeq1@sam_data$Treatment_Group, 
                                            levels(physeq1@sam_data$Treatment_Group)[c(1, 3, 2)])

saveRDS(physeq1, file = "FMT_physeq_raw.RDS")

#filtered, agglomerated
levels(physeq_glom@sam_data$Diagnosis)
physeq_glom@sam_data$Diagnosis  <- factor(physeq_glom@sam_data$Diagnosis, 
                                          levels(physeq_glom@sam_data$Diagnosis)[c(2, 1)])

levels(physeq_glom@sam_data$Collection_Time)
physeq_glom@sam_data$Collection_Time  <- factor(physeq_glom@sam_data$Collection_Time, 
                                                levels(physeq_glom@sam_data$Collection_Time)[c(2, 3, 1)])

levels(physeq_glom@sam_data$Treatment_Group)
physeq_glom@sam_data$Treatment_Group  <- factor(physeq_glom@sam_data$Treatment_Group, 
                                                levels(physeq_glom@sam_data$Treatment_Group)[c(1, 3, 2)])
levels(physeq_glom@sam_data$Treatment_Group)

saveRDS(physeq_glom, file = "FMT_physeq_glom.RDS")

#transformed to relative abundance
physeq_glom_ra <- transform_sample_counts(physeq_glom, function(x) x/sum(x))

#confirm transformation; column sums should = 1 
otu <-as.matrix(otu_table(physeq_glom_ra))
colsums <- colSums(otu)
which(colsums != 1)

levels(physeq_glom@sam_data$Diagnosis)
levels(physeq_glom@sam_data$Collection_Time)
levels(physeq_glom@sam_data$Treatment_Group)


saveRDS(physeq_glom_ra, file = "FMT_physeq_glom_ra.RDS")













