#!/bin/bash

## we performed 5'-3' gene body coverage analysis using RSeQC package \
## make sure it's properly installed in your system before running this.


source ~anaconda3/bin/activate rseqc
geneBody_coverage.py -r protein_coding_genes.bed -i Vasa_trim5000_sort.bam,Smartseq3_trim5000_sort.bam,MW3_trim5000_sort.bam,MW2_trim5000_sort.bam,flash_trim5000_sort.bam,drop_trim5000_sort.bam,10X_trim5000_sort.bam -o trim5000


## genes.bed is in bed12 format. We converted from filtered gtf by UCSC tools.
gtfToGenePred genes.gtf genes.genePred
genePredToBed genes.genePred genes.bed
## protein_coding_genes.bed file is constructed by greping the protein coding genes.