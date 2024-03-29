#!/bin/bash

cd `pwd`

workdir=YOURWORKDIR
scriptdir=YOURSCRIPTDIR
jaccard=${workdir}/2w_jaccard.csv ## generated by pool_jaccard_sample.py
pipeline=${workdir}/Jaccard_pipeline.R
outputDir=output

Rscript $pipeline $jaccard $outputDir $scriptdir