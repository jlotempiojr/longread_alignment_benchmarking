#!/bin/bash
mkdir trash
for vcf in *.vcf
do
	# remoive header and hash
	sed -i'.orig' '/##/d' $vcf
	sed -i'.proc.orig' 's/^#\(.*\)/\1/' $vcf
	# this snippet splits by info and then by semicolon
	cut -f 8 $vcf > info.$vcf
	awk  '{gsub(";","\t",$0); print;}' info.$vcf > tab_info.$vcf
	cut -f 2,10,12 tab_info.$vcf > paste_cut_tab_info.$vcf
	# pull out remaining columns of interst and store in a new vcf
	cut -f 1,2,7 $vcf> paste_cut_$vcf
	paste paste_* > processed.$vcf
	# remove header
	sed 1d processed.$vcf > processed.$vcf.tmp
	# add column for provenance filename
	sed -i '.temp' "s/$/,$vcf/" processed.$vcf.tmp
	# change tsv to csv
	awk  '{gsub("\t",",",$0); print;}' processed.$vcf.tmp > processed.$vcf.csv
	# add new header, new lines are a required character, not aesthetic
	sed -i '.bak' '1i\
CHROM,POS,FILTER,SVMETHOD,SVTYPE,SVLEN,PROV\
' processed.$vcf.csv
	#cleanup
	mv *.orig trash
	mv info* trash
	mv tab* trash
	mv paste* trash
	mv *.tmp trash
	mv *.temp trash
	mv *.bak trash
	mv processed*.vcf trash
	mv $vcf trash
done


