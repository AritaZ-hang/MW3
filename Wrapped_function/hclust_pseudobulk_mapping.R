#' @param target_pseudobulks The pseudobulks lists generated by hclust_pseudobulk_building.R & create_pseudobulks.R # nolint
#' @param type The tissue origin of pseudobulks. Default is "adj".
generateIdMapping <- function(target_pseudobulks, type = "adj")
{
  pseudobulks_out.id.mapping <- list()
  for(i in seq_along(target_pseudobulks))
  {
    celltype <- names(target_pseudobulks)[i]
    out.id.mapping <- data.frame(out.id = 1:max(target_pseudobulks[[i]][["final.index"]][["out.id"]]),
                                pseudobulks_name = paste0(celltype, "_", type, "_Pseudobulk_", 1:max(target_pseudobulks[[i]][["final.index"]][["out.id"]])))
    pseudobulks_out.id.mapping[[celltype]] <- out.id.mapping
  }
  return(pseudobulks_out.id.mapping)
}
#' @param target_pseudobulks The pseudobulks lists generated by hclust_pseudobulk_building & create_pseudobulks.R # nolint
#' @param pseudobulks_out.id.mapping The mapping list contains the pseudo-bulks names and the real cell names, generated by "generateIdMapping". # nolint
generateTargetIdx <- function(target_pseudobulks, pseudobulks_out.id.mapping)
{
  target_idx <- list()
  for(i in seq_along(target_pseudobulks))
  {
    final.index <- target_pseudobulks[[i]][["final.index"]]
    final.index[["cells"]] <- rownames(final.index)
    id.mapping <- pseudobulks_out.id.mapping[[i]]
    
    final.index <- dplyr::inner_join(final.index, id.mapping, by = "out.id")
    target_idx[[names(target_pseudobulks)[i]]] <- final.index
  }
  return(target_idx)
}