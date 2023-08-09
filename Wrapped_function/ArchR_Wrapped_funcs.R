## functions are wrapped for simplifying ArchR's pipeline

## FunctionsForConvenience: this function could complete LSI, clustering and UMAP in one row. Change the parameters to meet realistic requirements.
FunctionsForConvenience = function(proj, resolution, maxClusters, maxDim, iterations, outlierQuantiles, varFeatures, TileOrPeak = 'Tile')
{
  res = resolution
  
  if(TileOrPeak == 'Tile')
  {
    addIterativeLSI_name = paste0('IterativeLSI_bin_res',res)
    addClusters_name = paste0('Clusters_bin_res',res)
    addUMAP_name = paste0('UMAP_bin_res',res)
    use.matrix = 'TileMatrix'
  }
  else if(TileOrPeak == 'Peak')
  {
    addIterativeLSI_name = paste0('IterativeLSI_peak_res',res)
    addClusters_name = paste0('Clusters_peak_res',res)
    addUMAP_name = paste0('UMAP_peak_res',res)
    use.matrix = 'PeakMatrix'
  }
  tmp = proj
  # ArchR pipeline
  tmp = addIterativeLSI(ArchRProj = tmp, useMatrix = use.matrix, name = addIterativeLSI_name, iterations = iterations, clusterParams = list(resolution = res, sampleCells = NULL, maxClusters = maxClusters, n.start = 10), varFeatures = varFeatures, dimsToUse = 1:maxDim,force = T,
                        corCutOff = 0.5, outlierQuantiles = outlierQuantiles,
                        excludeChr = c('chrX','chrY')
  )
  tmp = addClusters(input = tmp, reducedDims = addIterativeLSI_name, name = addClusters_name, method = 'Seurat', resolution = res, force = TRUE, maxClusters = maxClusters, corCutOff = 0.5, dimsToUse = 1:maxDim, scaleDims = T, nOutlier = 70)
  tmp = addUMAP(ArchRProj = tmp, reducedDims = addIterativeLSI_name, name =addUMAP_name, nNeighbors = 30, minDist = 0.5, metric = "cosine",force = TRUE, corCutOff = 0.5, dimsToUse = 1:maxDim, scaleDims = T)
  return(tmp)
}

## GetMarkerList: This function will return a list of marker genes arranged in descending order by Log2FC. Significance is measured by FDR. Change necessary parameters to meet realistic requirements.
GetMarkerList = function(proj, res, FDR, Log2FC, type = 'Tile')
{
  if(type == 'Tile')
  {
    addClusters_name = paste0('Clusters_bin_res', res)
  }
  else if(type == 'Peak')
  {
    addClusters_name = paste0('Clusters_peak_res',res)
  }
  else
    addClusters_name = type
  filter = paste('FDR <=',FDR, '&' ,'Log2FC >=',Log2FC, sep = ' ')
  temp = getMarkerFeatures(ArchRProj = proj, useMatrix = 'GeneScoreMatrix', groupBy = addClusters_name, bias = c('TSSEnrichment', 'log10(nFrags)'), testMethod = 'wilcoxon')
  temp_list = getMarkers(temp, cutOff = filter)
  
  temp=temp_list@listData
  names(temp)=c(1:length(temp_list))
  use<-data.frame(matrix(ncol = 3,nrow = 0))
  colnames(use)=c('cluster','gene','Log2FC')
  for (i in 1:length(temp_list)) {
    if(dim(temp[[i]])[1]>0){
      tmp=data.frame(matrix(ncol = 3,nrow = dim(temp[[i]])[1]))
      colnames(tmp)=c('cluster','gene','Log2FC')
      tmp$cluster=i
      tmp$gene=temp[[i]]$name
      tmp$Log2FC=temp[[i]]$Log2FC
      tmp=tmp[order(-tmp$Log2FC),]}
    if(dim(temp[[i]])[1]==0){
      tmp=data.frame(matrix(ncol = 3,nrow = 1))
      colnames(tmp)=c('cluster','gene','Log2FC')
      tmp$cluster=i
      tmp$gene='null'
      tmp$Log2FC='NA'}
    use = rbind(use,tmp)}
  return(use)
}

## LabelTransfer: This function will return an ArchRproject with imputed gene integration matrix. Parameters @pbmc is your scRNA-seq seurat object, @anno_name is the metadata name which store the cell type annotations, @res is the resolution selected for the ArchRproject in previous analysis.
LabelTransfer = function(proj, pbmc, anno_name, res, TileorPeak = 'Tile')
{
  if(TileorPeak == 'Tile')
  {
    addReducedDims_name = paste0('IterativeLSI_bin_res', res)
  }
  else
    addReducedDims_name = paste0('IterativeLSI_peak_res', res)
  # if the thread number is not 1, there will be some problems
  proj_LabelTrans = proj
  proj_LabelTrans = addGeneIntegrationMatrix(
    ArchRProj = proj_LabelTrans, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = addReducedDims_name,
    seRNA = pbmc,
    addToArrow = F,
    groupRNA = anno_name,
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un",
    threads = 1,
    force = T
  )
  return(proj_LabelTrans)
}