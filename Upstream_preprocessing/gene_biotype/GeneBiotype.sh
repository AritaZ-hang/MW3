#!/bin/bash

## We use TagReadWithGeneFunction in Drop-seq tools to annotate reads with gene biotype \
## and summarize the file using customed R script.

~/Drop-seq_tools/TagReadWithGeneFunction I=MW3.bam O=MW3_gene.bam ANNOTATIONS_FILE=~/genome/hg19_mm10_transgenes.gtf

## then execute this to extract annotated biotype information from bam:
samtools view MW3_gene.bam | awk '{ for (i=1;i<=NF;++i){if($i ~"^gn:Z:"){print $i}}}'> MW3.txt

