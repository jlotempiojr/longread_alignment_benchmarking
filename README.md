# Longread alignment benchmarking

Scripts relevent to benchmarking longread sequence alignment tools.

## Citation

LoTempio, Jonathan; Délot, Emmanuèle; Vilain, Eric. Benchmarking long-read genome sequence alignment tools for human genomics applications. 11 July 2021. biorXiv. https://doi.org/10.1101/2021.07.09.451840

https://www.biorxiv.org/content/10.1101/2021.07.09.451840v1.full.pdf

## Data

https://github.com/jlotempiojr/longread_alignment_benchmarking/tree/main/data

## Scripts from the June 2022 version

### Processing the vcfs
These scripts will take sniffles vcfs generated on long read data and extract importand svlen and svtype fields for comparison to eachother or a standard:

#### PacBio:
https://github.com/jlotempiojr/longread_alignment_benchmarking/blob/main/vcf_processing_scripts.sh

#### Nanopore:
https://github.com/jlotempiojr/longread_alignment_benchmarking/blob/main/process_sniffles_nanopore_vcf.sh

### Wrangling the data for visualization
This example shows a pattern that one could use to prepare the vcfs and truths sets for vizualization
https://github.com/jlotempiojr/longread_alignment_benchmarking/tree/main/python_wrangling_example

### Scripts from versions 1, 2
These scripts are from earlier versions of this project, but retained for transperncy. The scripts are commented to outline the steps taken in VCF processing and utilize functionality in shell scripting and R. You can find them here:
https://github.com/jlotempiojr/longread_alignment_benchmarking/blob/db099c783099716704ff7ee5e53f68accea19253/supplemental_scripts.sh

**Contact**
evilain@childrensnational.org
