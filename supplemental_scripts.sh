#!/bin/bash

# Written between Friday March 1 2019 and February 2020
# Prepared for publication October 17-28
# Not executible, but includes a shebang for ease of reading in a text editor 
# This document should contain 5 sections. 

# Supplemental files can be found in this google drive folder
https://drive.google.com/drive/folders/1EImtHkYpOoUk7LOqPE__94r6CZSxk_fw?usp=sharing

## R packages required for many scripts
library(tidyverse)
library(stats)
library(cutr)
library(VennDiagram)

########################################################
## Section one: break down VCF and get a column of each relevent part 
## This section follows a theoretical example file, minimap_pacbio.vcf, as it is reshaped into a usable file, minimap2_pacbio_processed.csv
## Taking the files in the google drive folder Raw VCFs through the steps in sections 1 and 2 will yeild the results found in Processed VCFs

#### in BASH
	# remoive header and hash
	sed -i.bak '' '/##/d' minimap_pacbio.vcf
	sed -i.bak ''  's/^#\(.*\)/\1/' minimap_pacbio.vcf

	# this snippet splits by info and then by semicolon
	cut -f 8 minimap_pacbio.vcf > info_minimap_pacbio.vcf
	awk  '{gsub(";","\t",$0); print;}' info_minimap_pacbio.vcf > tab_info_minimap_pacbio.vcf
	cut -f 1-13 tab_info_minimap_pacbio.vcf > paste_cut_tab_info_minimap_pacbio.vcf

	# pull out remaining columns of interst and store in a new vcf
	cut -f 1,2,5,10 minimap_pacbio.vcf > paste_cut_minimap_pacbio.vcf

	paste paste_* > minimap2_pacbio_processed.vcf

	# remove header
	tail -n +2 "minimap2_pacbio_processed.vcf" > "minimap2_38_pacbio_processed.vcf.tmp" && mv "minimap2_pacbio_processed.vcf.tmp" "minimap2_pacbio_processed.vcf"

	# add column for provenance filename
	# sed -i "s///"

	sed -i '' 's/$/ pacbio_minimap_sniffles/' minimap2_pacbio_processed.vcf

	# change tsv to csv

	awk  '{gsub("\t",",",$0); print;}' sniffles_sorted_minimap2_38_pacbio_processed.vcf > minimap2_pacbio_processed.csv

	# add new header

	sed -i.bak '1i\
	CHROMPOS, ALT, minimap2pb38, INFO, SVMethod, CHR_2, POS_END, STD_quant_start, STD_quant_stop, Kurtosis_quant_start, Kurtosis_quant_stop, SVTYPE, SUPTYPE, SVLEN, STRANDS, RE, 	Pipeline\
	' minimap2_pacbio_processed.csv

	# results in a CSV with separate columns for SVTYPE and SVLEN, as well as a new column for alignment tool
	
	# repeat for each .vcf and cat  *.csv > [platform].all.csv

########################################################
## Section two: change RefSeq to UCSC style name (human readable) and filter for varients found in the main assembly. Bin by SVLEN
## This section follows the processing of a theoretical file containing all_processed.vcfs for one given platform's technology, in our case data from ONT or SMRT cells

#### in R ####

	#read in files
	PB = read.csv('[platform].all.csv', header=T)
	ref = read.csv('chr_key.csv', header=T)
	
	#add human readable chromosome names and save
	PB$Data_set = colnames()[10]
	PB$UCSC_style__name = merge(PB, ref, by="CHROM")
	write.table(PB, "PB_processed.txt", sep="\t")
	 
	#only include main assembly and m variants
	PB_main_assembly = filter(PB, UCSC_style__name == "chr1" | UCSC_style__name == "chr2" | UCSC_style__name == "chr3"  UCSC_style__name == "chr4" | UCSC_style__name == "chr5" | 	 	 	UCSC_style__name == "chr6" | UCSC_style__name == "chr7" | UCSC_style__name == "chr8" | UCSC_style__name == "chr9" | UCSC_style__name == "chr10" | UCSC_style__name == "chr11" | 	 	 	UCSC_style__name == "chr12" | UCSC_style__name == "chr13" | UCSC_style__name == "chr14" | UCSC_style__name == "chr15" | UCSC_style__name == "chr16" | UCSC_style__name == "chr17" | 	 	UCSC_style__name == "chr18" | UCSC_style__name == "chr19" | UCSC_style__name == "chr20" | UCSC_style__name == "chr21" | UCSC_style__name == "chr22" | UCSC_style__name == "chrX" | 	 	UCSC_style__name == "chrY" | UCSC_style__name == "chrM")
	
	#sort by chromosome
	sorted_PB_main_assembly = PB_main_assembly[order(as.numeric(as.character(PB_main_assembly$UCSC_style_name)))]
	
	#take absolute value of SVLEN
	sorted_PB_main_assembly$ABS_SVLEN = abs(sorted_PB_main_assembly$SVLEN)
	
	#add bins of ABS_SVLEN
	sorted_PB_main_assembly$BINS = smart_cut(sorted_PB_main_assembly$ABS_SVLEN,c(50, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 20000, 30000, 40000, 50000, 60000, 	70000, 80000, 90000, 100000),labels = ~paste(sep="~",.y[1],round(mean(.x),2),.y[2]), "breaks", simplify = TRUE)

########################################################
## Section three: Prepare truthset for comparison
## This section follows the downloadable file from the DGV consortium as it is prepared for comparison to data from vcfs
## Taking GRCh38_hg38_variants_2016-08-31.txt in the google drive folder through these steps will yeild a file equivalent to dgv_truth_na12878.txt

#### in BASH ####
	#filter only variants present in NA12878
	grep -o  '^[NA12878]*$' GRCh38_hg38_variants_2016-08-31.txt > dgv_truth_na12878.txt

#### in R ####
	truth = read.csv("dgv_truth_na12878.txt", header=T, sep="\t")
	
	#add columns equivalent to those created for VCFs
	truth$UCSC_style_name = "chr"
	truth$UCSC_style_name = paste(truth$UCSC_style_name, truth$chr, sep='')
	truth$Data_set = "DGV_dataset"
	truth$SVLEN = mutate(truth, truth$end - truth$start)
	truth$BINS = smart_cut(truth$SVLEN,c(50, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000),labels = 	~paste(sep="~",.y[1],round(mean(.x),2),.y[2]), "breaks", simplify = FALSE)
	
	#include only main assembly and mt variants
	final_na12878_dgv_truth = subset(truth, truth$chr == "1" |truth$chr == "15" | truth$chr == "16" | truth$chr == "17" | truth$chr == "18" | truth$chr == "19" | truth$chr == "2" | truth$chr == "20" | truth$chr == "21" | truth$chr == "22" | truth$chr == "3" | truth$chr == truth$chr == "4" | truth$chr == "5" | truth$chr == "6" | truth$chr == "7" | truth$chr == "8" | truth$chr == truth$chr == "9" | truth$chr == "X" | truth$chr == "Y")
	
	#check for the presence of alt loci variants in truthset, should return empty dataframe
	altloci_truth = subset(truth, truth$chr == "17_KI270857v1_alt" | truth$chr == "4_GL000008v2_random" |truth$chr == "8_KI270821v1_alt" | 	truth$chr == "14_GL000009v2_random" | truth$chr == "Un_KI270742v1" | truth$chr == "19_KI270938v1_alt")

########################################################
## Section four: generate flat contigency table for variant frequency, and graph
## This section shows scripts which generate a frequency table based on the fields required to probe alignment tool, variant size, and type
## Completion of these steps will yeild a file equivalen to all_freq_table.txt in the google drive folder

#### in R ####

	test_pb_table = table(sorted_ONT_main_assembly$Alignment_tool, sorted_ONT_main_assembly$BINS, sorted_ONT_main_assembly$SVTYPE)
	ftable(test_pb_table)
	write.table(test_pb_table, "pb_freq_table.txt", sep="\t")
	
	# repeat for each file

#### in bash ####
	# cat all frequenct tables and inspect, graph in excel

	cat *freq_table.txt > all_freq_table.txt

	# open in excel for inspection and graphing with pivot charts

########################################################
## Section five: overlapping unmapped reads
## This section examines the lists of readnames from unaligned reads in a given bam from a single dataset to find overlap for representation in a Venn diagram

## BAMs are available upon request

	# Extract unmapped reads with samtools
	samtools view -f4 query.bam > unmapped.bam | cut -f1 > platform_alignmenttool.qname.txt
	
	#compare unmapped read names to determin overlap
	#### in R ####
	# a, b, and c are lists of qnames from three different bams containing unaligned read from an alignmen experiment
	a_b = intersect(a,b)
	a_c = intersect(a,c)
	b_c = intersect(b,c)
	
	a_b_c = intersect(a_b, c)
	b_c_a = intersect(b_c, a)
	
	all_intersect = intersect(a_b_c, b_c_a)
	
	len(a_b)
	len(a_c)
	len(b_c)
	len(a_b_c)
	len(b_c_a)
	len(all_intersect)
	
	Venn_one = draw.triple.venn(2.637, 2.675, 3.288, 2.021, 2.124, 2.096, 1.610, 
																				scaled = TRUE, 
																				fill = c("#3B9AB2", "#EBCC2A", "#F21A00"), 
																				alpha = 0.5, filename = "venn_test.tiff", 
																				col = "transparent", 
																				fontfamily = "sans", 
																				main.cex = 40)
