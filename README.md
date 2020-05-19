# Gal-Barko-K9-AHDS-FMT
Analysis of the effect of fecal microbiota transplantation (FMT) on the fecal microbiome of dogs with acute hemorrhagic diarrhea syndrome (AHDS)

This repository contains the data and R code to reproduce the results in the manuscript titled "Characterization of donor microbiota signature following fecal microbiota transplantation in dogs with acute hemorrhagic diarrhea syndrome" and submitted to Animal Microbiome.

The raw sequencing data were deposited in the National Center for Biotechnology Information (NCBI) Short Read Archive (SRA) under the BioProject accession PRJNA629659 with the BioSamples in the range SAMN14785809 to SAMN14785844. For full reproducibility, these data can accessed and preprocessed using the shell script: 

The demultiplexed, trimmed, and aligned reads, along with sample and taxonomy tbales have been deposited here for convenience.

Data:

1. FMT_OTU_table.csv
2. FMT_taxa_table.csv
3. FMT_sample_data.csv

Helper Functions:

1. SumarizeFit.R

To replicate our analysis, execute the R scripts in the following order:

1. Gal_Barko_FMT_Preprocessing.R
2. Gal_Barko_FMT_Diversity.R
3. Gal_Barko_FMT_Diff_Abundance.R
4. Gal_Barko_FMT_WGCNA.R

