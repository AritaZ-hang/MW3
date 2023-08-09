setwd("YOURWORKDIR")
library(data.table)
library(stringr)

cytoband <- as.data.frame(fread(file = "./cytoBand.txt.gz"))[,1:4]
names(cytoband) <- c("chr", "start", "end", "cytoband")
cytoband$cytoband <- str_remove(cytoband$cytoband, "[0-9.]+")
cytoband$merge <- paste0(cytoband$chr, cytoband$cytoband)
cytoband <- cytoband[mixedorder(cytoband$merge), ]

merged_cytoband_name <- unique(cytoband$merge)
cytoband_range <- list()
for(i in merged_cytoband_name)
{
  range_a <- cytoband[which(cytoband$merge == i),]
  start <- min(range_a[["start"]])
  end <- max(range_a[["end"]])
  res <- list()
  res[["cytoband"]] <- i
  res[["chr"]] <- unique(range_a[["chr"]])
  res[["start"]] <- start
  res[["end"]] <- end
  cytoband_range[[i]] <- res
}
cytoband_range_ls <- as.data.frame(do.call(rbind, cytoband_range))
cytoband_range_ls[["start"]] <- as.numeric(cytoband_range_ls[["start"]])
cytoband_range_ls[["end"]] <- as.numeric(cytoband_range_ls[["end"]])
cytoband_range_ls[["cytoband"]] <- as.character(cytoband_range_ls[["cytoband"]])
cytoband_range_ls[["chr"]] <- gsub("chr", "", cytoband_range_ls[["chr"]])
saveRDS(cytoband_range_ls, file = "./cytoband_range.rds")