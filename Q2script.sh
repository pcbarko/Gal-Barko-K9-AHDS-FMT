#!/bin/bash
#
#	A bash script to run QIIME2 for basic analysis using the manifest.txt file 
#	input method for QIIME2 data import
#	based on the metafile which has been previously prepared
#
#	The results folder is presumed to have been previously made and is a 
#	subfolder of the project folder
#	The manifest and metaData files are both in the project folder.
#
# 	the script is run from the created results folder  
#
#	created by pjb on: 		2017-05-25
#	last edited by pjb on:	2020-06-04
#
#	email: 	p.biggs@massey.ac.nz
#	GitHub:	https://github.com/pjbiggs/
#
#
# 	this is based on 2019.1 code
#
#
##################################################


## move to the results folder

cd /path_to_results_folder




## import the data

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path ../manifest.txt --output-path demux-paired-end.qza --input-format PairedEndFastqManifestPhred33

qiime demux summarize --i-data demux-paired-end.qza --o-visualization demux-paired-end.qzv




## perform QC and analyse after manual inspection of the QC data to choose trimming options  

qiime dada2 denoise-paired --i-demultiplexed-seqs demux-paired-end.qza --o-table table_tf10_tr20_truf250_trur250 --o-representative-sequences rep-seqs_tf10_tr20_truf250_trur250 --p-trim-left-f 10 --p-trim-left-r 20 --p-trunc-len-f 250 --p-trunc-len-r 250 --o-denoising-stats denoising-stats.qza

qiime feature-table tabulate-seqs --i-data rep-seqs_tf10_tr20_truf250_trur250.qza --o-visualization rep-seqs_tf10_tr20_truf250_trur250.qzv
 
qiime feature-table summarize --i-table table_tf10_tr20_truf250_trur250.qza --o-visualization table_tf10_tr20_truf250_trur250.qzv --m-sample-metadata-file ../metaData.txt

qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization denoising-stats.qzv




## generate the tree

qiime alignment mafft --i-sequences rep-seqs_tf10_tr20_truf250_trur250.qza --o-alignment aligned-rep-seqs.qza

qiime alignment mask --i-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza 

qiime phylogeny fasttree --i-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza

qiime phylogeny midpoint-root --i-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza




## perform the core diversity

qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table_tf10_tr20_truf250_trur250.qza --p-sampling-depth <<chosen_value>> --m-metadata-file ../metaData.txt --output-dir ../core-metrics-results-test 

qiime diversity alpha-group-significance --i-alpha-diversity ../core-metrics-results-test/faith_pd_vector.qza --m-metadata-file ../metaData.txt --o-visualization ../core-metrics-results-test/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance --i-alpha-diversity ../core-metrics-results-test/evenness_vector.qza --m-metadata-file ../metaData.txt --o-visualization ../core-metrics-results-test/evenness-group-significance.qzv




## taxonomy 

qiime feature-classifier classify-sklearn --i-classifier /path_to_greengenes_folder/gg-13-8-99-515-806-nb-classifier.qza --i-reads rep-seqs_tf10_tr20_truf250_trur250.qza --o-classification taxonomy.qza

qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv
  
qiime taxa barplot --i-table table_tf10_tr20_truf250_trur250.qza --i-taxonomy taxonomy.qza --m-metadata-file ../metaData.txt --o-visualization taxa-bar-plots.qzv    


