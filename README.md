# MW3
This is the code repository for the bioinformatical analyses used in MW3 project. 
TODO: Files will be organized into seperate folders. 
## To merge beads
Pipeline are adapted from bap2.
1. Jaccard_pipeline.R
2. Jaccard_plotting.R
3. Cutoff.R
4. Knee.R
5. pool_jaccard_sample.py
## To perform basic analysis
We've provided simplified software pipelines to perform basic scATAC-seq/scRNA-seq analysis. For the use of parameters, please check our methods & the softwares' official tutorials.
1. ATAC_ArchR_pipeline.R
2. RNA_Seurat_pipeline.R
## Benchmark with other methods
Benchmark with 3'-seq data (Mouse Cell Atlas, 10X Genomics, Tabula Muris). 
1. Genebody_coverage.sh
2. GeneBiotype.sh
3. Brain_diffGenes.R
4. Brain_uniqueGenes.R
## To construct gene regulatory network
We've employed SCRIP & SCENIC+ to conduct this part's analysis. For the use of parameters, please check our methods & the softwares' official tutorials.
1. SCRIP.sh
2. SCRIP_enrich.R
3. SCENICplus_preprocessing.ipynb
4. SCENICplus_pipeline.ipynb
5. Oligo_pseudotime.R
6. SCRIP_target.R
## To infer copy number alterations
InferCNV, CopyscAT & Alleloscope were used to infer large-scale copy number alterations in scRNA-seq & scATAC-seq respectively. 
1. inferCNV_pipeline.R
2. IterativeHclust.R
3. CopyscAT_refinedFunctions.R
4. CopyscAT_pipeline.R
5. Alleloscope_preprocessing.sh
6. Alleloscope_heatmap.R
## Wrapped function
1. ArchR_Wrapped_func.R
