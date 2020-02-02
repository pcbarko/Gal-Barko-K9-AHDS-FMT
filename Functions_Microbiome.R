
#FUNCTIONS FOR ANALYSIS OF 16S-rRNA SEQUENCING STUDIES

#Normalization ====

#TSS normalization (via mixOmics). 

#TSS normalisation accommodates for varying sampling and sequencing depth. 
#The read count is divided by the total number of read counts in each individual sample. 
#TSS results in compositional data (or proportions) that are restricted to a space 
#where the sum of all OTU proportions for a given sample sums to 1. 
#Using standard statistical methods on such data may lead to spurious results 
#and therefore the data must be further transformed.

TSS.divide = function(x){
  x/sum(x)
}

#Wrapper Functions for Differential Abundance Analysis Using DESeq2 ====

##function to extract results from specific comparisons and output to file 
generateContrastResults <- function(deseq, contrast, physeq, file = "results.txt") {
  res = results(deseq, 
                cooksCutoff = FALSE, 
                contrast = contrast,
                test = "Wald"
  )
  res = cbind(as(res, "data.frame"), as(tax_table(physeq)[rownames(res), ], "matrix"))
  
  write.table(res, file, row.names = TRUE, col.names = NA, sep = "\t")
  return(res)
}

##plotting functions ====


plotDESeq2Genus <- function(res) {
  sigtab <- res[which(res$padj < 0.05), ]
  scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
  }
  
  x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
  
  # Genus order
  x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
  p <- ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
  return(p)
}




plotDESeq2Family <- function(res) {
  sigtab <- res[which(res$padj < 0.05), ]
  scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
  }
  
  x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
  
  # Family order
  x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))
  p <- ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Family)) + geom_point(size=6) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
  return(p)
}

