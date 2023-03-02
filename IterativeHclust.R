## SET UP WORKING ENVIRONMENT ##
rm(list = ls()); gc()
set.seed(1)
setwd(dir = 'YOURWORKDIR')

## ref_table & obs_table: generated in the 1st round of infercnv with cells from tumor-adjacent tissues setting as reference.

ref_table = read.table('./result/infercnv.references.txt', header = T)
obs_table = read.table('./result/infercnv.observations.txt', header = T)

processMatrix = function(table)
{
  cnv_score_mat = as.matrix(table)
  cnv_score_table = as.matrix(table)
  
  ## measure variation
  cnv_score_table[cnv_score_mat >= 0.85 & cnv_score_mat < 0.95] = 1 
  cnv_score_table[cnv_score_mat >= 0.95 & cnv_score_mat < 1.05] = 0
  cnv_score_table[cnv_score_mat >= 1.05 & cnv_score_mat <= 1.15] = 1
  
  res = list(cnv_score_mat, cnv_score_table); names(res) = c('mat', 'table')
  return(res)
}

MeanExpect = function(table)
{
  t1 = processMatrix(table)
  ## calculate mean cnv score 
  expect = floor(sum(rowSums(t1$table)/dim(t1$table)[2]))
  print(paste('mean cnv score expect:' ,expect))
  return(expect)
}

norm_expect = MeanExpect(ref_table)
mali_expect = MeanExpect(obs_table)

CellDiscern = function(table, k_groups, iteration = 2)
{
  mat_table = processMatrix(table)
  
  ## hclust to define less variable groups
  t1 = Reclustering(mat_table$mat, cells=colnames(mat_table$mat), clust_m = 'ward.D', k_groups = k_groups)
  
  ## judgment step
  res = OutlierCapture(mat_table$table, t1[[1]], t1[[2]])
  norm_cells = res$norm
  mali_cells = res$mali
  int_cells_tmp = res$int
  
  if(length(int_cells_tmp) > 0)
  {
    count = 0
    while(count < iteration)
    {
      print(paste('not converge, iteration step:', as.character(count+1)))
      t2 = Reclustering(mat_table$mat, cells = int_cells_tmp, clust_m = 'ward.D', k_groups = k_groups)
      tmp = OutlierCapture(mat_table$table, t2[[1]], out.id = t2[[2]])
      print(paste(as.character(count+1), "round discerned normal cells:", as.character(length(tmp$norm))))
      print(paste(as.character(count+1)," round discerned malignant cells:", as.character(length(tmp$mali)))); 
      print(paste(as.character(count+1)," round discerned temporary intermediate cells:", as.character(length(tmp$int))))
      norm_cells = c(norm_cells, tmp$norm)
      mali_cells = c(mali_cells, tmp$mali)
      if(length(tmp$int) == 0)
      {
        final_res = list(norm_cells, mali_cells); names(final_res) = c('pdNormal', 'pdMalignant')
        return(final_res)
      }
      else if(length(tmp$norm) == 0 & length(tmp$mali) == 0)
      {
        final_res = list(norm_cells, mali_cells, int_cells_tmp); names(final_res) = c('pdNormal', 'pdMalignant', 'pdIntermediate')
        return(final_res)
      }
      else
      {count = count+1
      int_cells_tmp = tmp$int}
    }
    final_res = list(norm_cells, mali_cells, int_cells_tmp); names(final_res) = c('pdNormal', 'pdMalignant', 'pdIntermediate')
    return(final_res)
  }
  else
  {
    final_res = list(norm_cells, mali_cells); names(final_res) = c('pdNormal', 'pdMalignant')
    return(final_res)
  }
}
## hclust step
Reclustering = function(cnv_score_mat, cells, clust_m = 'ward.D', k_groups = 10)
{
  mat = cnv_score_mat[, cells]
  dists = parallelDist::parallelDist(t(mat), method = 'euclidean')
  t = fastcluster::hclust(dists, method = clust_m)
  out.id = cutree(t, k = k_groups)
  new.df = as.data.frame(out.id)
  return(list(new.df, out.id))
}
## @params: type: nonmalignant or malignant
OutlierCapture = function(cnv_score_table, new.df, out.id)
{
  cnv_score_table = cnv_score_table[, rownames(new.df)]
  ## compute average 
  cell_scores_CNV = as.data.frame(colSums(cnv_score_table)) ## total cnv score for each cell
  colnames(cell_scores_CNV) = 'cnv_score'
  new.df = new.df[rownames(cell_scores_CNV),]
  cell_scores_CNV$group = new.df
  a1 = aggregate(cell_scores_CNV$cnv_score, by = list(cell_scores_CNV$group), sum) ## sum of cnv scores for each hclust group
  groups = as.data.frame(table(out.id))
  group_avg = a1$x / groups$Freq ## average cnv score for each hclust group
  print(group_avg)
  
  idx_norm = seq(1:length(group_avg))[group_avg <= norm_expect]
  idx_mali = seq(1:length(group_avg))[group_avg >= mali_expect]
  idx_int_tmp = seq(1:length(group_avg))[group_avg > norm_expect & group_avg < mali_expect]
  
  norm_cells = rownames(cell_scores_CNV[cell_scores_CNV$group %in% idx_norm,])
  mali_cells = rownames(cell_scores_CNV[cell_scores_CNV$group %in% idx_mali,])
  int_tmp_cells = rownames(cell_scores_CNV[cell_scores_CNV$group %in% idx_int_tmp,])
  
  res = list(norm_cells, mali_cells, int_tmp_cells);names(res) = c('norm', 'mali', 'int')
  return(res)
}

ref_discern_K20 = CellDiscern(ref_table, k_groups = 20, iteration = 10)
obs_discern_K20 = CellDiscern(obs_table, k_groups = 20, iteration = 10)


K20_pdNormal = c(ref_discern_K20$pdNormal, obs_discern_K20$pdNormal); K20_pdMalignant = c(ref_discern_K20$pdMalignant, obs_discern_K20$pdMalignant); K20_pdIntermediate = c(ref_discern_K20$pdIntermediate, obs_discern_K20$pdIntermediate) ## ..

K20 = list(K20_pdNormal, K20_pdMalignant, K20_pdIntermediate); names(K20) = c('pdNormal', 'pdMalignant', 'pdIntermediate')

saveRDS(K20, file = './K20.rds')
