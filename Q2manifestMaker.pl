#!/usr/bin/perl
#
#	A Perl script to make a manifest.txt file for QIIME2 data import
#	based on the metafile which has been previously prepared
#
#	3 variables need to be set for the script to work
#	
#	$root:		base folder for data analysis
#	$seqRoot:	folder where fastq sequences are found
#	$inMeta:	QIIME2-ready metafile with the first column containing 
#					the sample names as per the fastq sequence files 
#
#	created by pjb on: 		2017-05-25
#	last edited by pjb on:	2020-06-04
#
#	email: 	p.biggs@massey.ac.nz
#	GitHub:	https://github.com/pjbiggs/
#
##################################################

use warnings;
use strict;

#######################

# the 3 variables to be set for the analysis #

my $root 		= ("/path_to_data_folder/");
my $seqRoot		= ($root . "subfolder_where_sequences_are_located/");
my $inMeta		= ($root . "your_metafile.txt");


## make a log file ##
my $manifest	= ($root . "manifest_inGuts.txt");
my $log 		= ($root . "makingManifest_file.txt");

open (LOG, ">$log") or die ("couldn't open $log: $!\n");

print ("Process started at " . scalar(localtime) . ".\n");
print LOG ("Process started at " . scalar(localtime) . ".\n");


## do the work ##

open (LIST, "<$inMeta") or die ("couldn't open $inMeta: $!\n");
open (OUT, ">$manifest") or die ("couldn't open $manifest: $!\n");

print OUT ("sample-id,absolute-filepath,direction\n");

while (<LIST>) {
	chomp;
	my @data	= split;
	my $sample	= $data[0];

	if ($sample ne '#SampleID') {	
		print OUT ($sample . "," . $seqRoot . $sample . "_L001_R1_001.fastq.gz,forward\n");
		print OUT ($sample . "," . $seqRoot . $sample . "_L001_R2_001.fastq.gz,reverse\n");		
	}	
}

close LIST;
close OUT;

print ("All finished at " . scalar(localtime) . ".\n");
print LOG ("All finished at " . scalar(localtime) . ".\n");

close LOG;
