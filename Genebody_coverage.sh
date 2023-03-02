#!/bin/bash

## we performed 5'-3' gene body coverage analysis using RSeQC package \
## make sure it's properly installed in your system before running this.

geneBody_coverage.py -r genes.bed -i 10x.bam,mw3.bam,smartseq.bam,vasa.bam -o genebody

## genes.bed is in bed12 format. We converted from filtered gtf by UCSC tools.
gtfToGenePred genes.gtf genes.genePred
genePredToBed genes.genePred genes.bed
## using genes.bed for RSeQC