#!/bin/bash

## Before running this, please make sure SCRIP is properly installed in your laptop \
## and processing the extracted tile matrix into h5ad format.
## To obtain enrichment score matrix:
SCRIP enrich -i use.h5ad -s mm -p Brain

## To investigate potential target genes of any interested TF:
SCRIP impute -i use.h5ad -s mm -p nkx2_1 --factor Nkx2-1
SCRIP target -i Nkx2-1_imputed.h5ad -s mm -o Nkx2_1.h5ad

## using Nkx2_1.h5ad for SCRIP_target.R