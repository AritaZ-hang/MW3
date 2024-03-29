#' @param meta The metadata of the cells.
#' @param kn_wanted The number of pseudobulks to group.
#' @param annotation_key The key that stores cell type info in the metadata.
getTargetPseudobulk_num <- function(meta, kn_wanted, annotation_key)
{
  ncells_per_celltype <- as.data.frame(table(meta[[annotation_key]]))
  colnames(ncells_per_celltype) <- c("celltypes", "ncells")
  rownames(ncells_per_celltype) <- ncells_per_celltype[["celltypes"]]
  ceiling_pseudobulk_ncells <- as.integer(max(ncells_per_celltype[["ncells"]]) / kn_wanted)
  print("Following cell types are excluded for pseudobulk construction...")
  for(i in seq_len(ncells_per_celltype)[1])
  {
    if(ncells_per_celltype[i, "ncells"] < ceiling_pseudobulk_ncells)
      print(rownames(ncells_per_celltype)[i])
  }
  retained <- rownames(ncells_per_celltype)[which(ncells_per_celltype[["ncells"]] >= ceiling_pseudobulk_ncells)]
  target_pseudobulk_nums <- list()
  for(i in retained)
  {
    target_pseudobulk_nums[i] <- floor(ncells_per_celltype[i, "ncells"] / ceiling_pseudobulk_ncells)
  }
  return(target_pseudobulk_nums)
}
#' @param dge The cells dge that will construct pseudo-bulks.
#' @param kn The number of pseudobulks to group. 
hclust_meta <- function(dge, kn)
{
  eucli_d <- parallelDist::parallelDist(t(as.matrix(dge)), method = "euclidean")
  t <- fastcluster::hclust(eucli_d, method = "ward.D")
  out.id <- cutree(t, k = kn)
  new.df <- as.data.frame(out.id)
  return(new.df)
}
#' @param dge The cells dge that will construct pseudo-bulks.
#' @param kn The number of pseudobulks to group.
createPseudobulk <- function(dge, kn)
{ 
  new.df <- hclust_meta(dge, kn)
  pseudobulk <- c()
  for(j in seq_len(kn))
  {
    cells_to_use <- rownames(new.df)[which(new.df[, 1] == j)]
    tmp_counts <- dge[, cells_to_use]
    if(length(cells_to_use) > 1)
    {
      tmp <- as.vector(Matrix::rowSums(tmp_counts))
      raw.10e6 <-  as.vector(t(t(tmp) / sum(tmp)) * 1e6)
    }
    else
      raw.10e6 <- as.vector(t(t(tmp_counts)/sum(tmp_counts))*1e6)
    pseudobulk <- cbind(pseudobulk, raw.10e6)
  }
  pseudobulk <- as.data.frame(pseudobulk)
  rownames(pseudobulk) <- rownames(dge)
  names(pseudobulk) <- paste0("Pseudobulk_", 1:kn)
  res <- list(new.df, pseudobulk)
  names(res) <- c("final.index", "pseudobulk")
  return(res)
}
#' @param dge Merged obs & ref CNV score matrix generated by inferCNV. # nolint
#' @param meta The metadata of the cells.
#' @param annotation_key The key that stores cell type info in the metadata.
#' @param target_num_list The target pseudobulk number generated by "getTargetPseudobulk_num". # nolint
#' @param type The tissue origin of pseudobulks. Default is "adj".
generatePseudobulk <- function(dge, meta, annotation_key, target_num_list, type = "adj")
{
  pseudobulk_target <- list()
  for(i in names(target_num_list))
  {
    print(paste0("Processing:", i))
    target_cells <- rownames(meta)[which(meta[[annotation_key]] == i)]
    dge_target <- dge[, target_cells]
    target_pseudobulks <- createPseudobulk(dge_target, kn = target_num_list[[i]])
    colnames(target_pseudobulks[["pseudobulk"]]) <- paste0(i, "_", type, "_", colnames(target_pseudobulks[['pseudobulk']]))
    pseudobulk_target[[i]] <- target_pseudobulks
  }
  return(pseudobulk_target)
}
