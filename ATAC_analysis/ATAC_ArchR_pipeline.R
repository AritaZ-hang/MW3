## SET UP WORKING ENVIRONMENT ##
rm(list = ls()); gc()
setwd('YOURWORKDIR')
source('ArchR_Wrapped_func.R')
library(tidyverse)
library(data.table)
library(stringr)
library(ArchR)
addArchRGenome('mm10')
library(dplyr)
library(magrittr)

## STEP 01 Create Arrow Files ##
file.list = list.files(pattern = '.bed.gz$')
sample.names = gsub('.bed.gz', '', file.list)

### Default Paramaters for QC ###
ArrowFiles = createArrowFiles(inputFiles = file.list,sampleNames = sample.names,minTSS = 7, TileMatParams = list(tileSize = 5000),GeneScoreMatParams = list(tileSize = 2000), addTileMat = TRUE, addGeneScoreMat = TRUE,minFrags = 1000)

## STEP 02: Create ArchR project, do not set 'copyArrows' = TRUE cuz it's extremely time-consuming! ##

ArrowFiles = list.files(pattern = '*.arrow')
proj = ArchRProject(ArrowFiles = ArrowFiles,copyArrows = F)

## STEP 03: Perform LSI & clustering & UMAP in one row using the wrapped function ##

proj = FunctionsForConvenience(proj, resolution = 1, maxClusters = NULL, maxDim = 30, iterations = 2, outlierQuantiles = c(0.98, 0.02), varFeatures = 25000, TileOrPeak = 'Tile')

## STEP 04: Obtain markers for each identified cluster ## 
markers = GetMarkerList(proj, res = 0.5, FDR = 0.1, Log2FC = 1, type = 'Tile')
