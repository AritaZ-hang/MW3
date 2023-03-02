## SET UP WORKING ENVIRONMENT ##

rm(list = ls()); gc()
setwd('YOURWORKDIR')
library(tidyverse)
library(data.table)
library(stringr)
library(ArchR)
addArchRGenome('mm10')
library(dplyr)
library(magrittr)
library(Alleloscope)
dir_path = './Alleloscope'; 

## STEP 01: Data Preparation ## 

proj_cnv = readRDS(file = 'proj_cnv.rds')
TileMatrix = getMatrixFromProject(proj_cnv, useMatrix = 'TileMatrix', binarize = T)
tilemat = TileMatrix@assays@data@listData$TileMatrix
rownames(tilemat) <- paste(TileMatrix@elementMetadata@listData$seqnames,"-",
                           TileMatrix@elementMetadata@listData$start,"-",
                           TileMatrix@elementMetadata@listData$start+4999,
                           sep = "") 


size = read.table(file='mm10.chrom.sizes')
size = size[size$V1 %in% paste0('chr', 1:19),]

# SNPs by cell matrices for ref and alt alleles.
barcodes = fread('lung_bc.tsv', header = F)
alt_all = readRDS(file = 'alt_all.rds')
ref_all = readRDS(file = 'ref_all.rds')
var_all = readRDS(file = 'var_all.rds')


## so basically the colnames of these dgcMatrix are in accordance with the sequence of barcodes.
alt_all@Dimnames[2] = barcodes[,1]
ref_all@Dimnames[2] = barcodes[,1]

raw_counts = tilemat
meta_info = as.data.frame(proj_cnv@cellColData)
meta_info$Barcode = cells.new
tissue_info = meta_info[, c('Barcode', 'tissue')]
rownames(tissue_info) = tissue_info$Barcode

## randomly select 3000 cells to plot
selected.barcodes = sample(colnames(raw_counts), 3000)
raw_counts.sep = raw_counts[, selected.barcodes]
tissue_info.sep = tissue_info[selected.barcodes,]

## STEP 02: Pheatmap plotting
## adapted from 'Plot_scATAC_snv' function

raw_mat = raw_counts.sep
cell_type = tissue_info.sep
normal_lab = 'WT'
size = size
window_w = 10000000
window_step = 2000000
plot_path = NULL
nclust = 4
var.filter = FALSE

if(is.null(plot_path)){
  plot_path=paste0("./plots/CNV_cov_w",window_w,"_s",window_step,"_sub.pdf")
  dir.create("./plots/")
}
## normlize by cell size
if(var.filter==TRUE){
  vars=apply(raw_mat,1, var)
  cellsize=Matrix::colSums(raw_mat[which(vars<quantile(vars,0.99)),])
}else{
  cellsize=Matrix::colSums(raw_mat)
}
cellsize=Matrix::Matrix(rep(cellsize, nrow(raw_mat)), byrow =T, ncol=ncol(raw_mat))
raw_mat=raw_mat/cellsize
if(grepl("chr",size[1,1])){
  size=size
}else{
  size[,1]=paste0("chr",size[,1])
}
size=size[ (size[,1] %in% paste0('chr',1:19)),] ## Mus musculus
cnv_bin0=GRanges(size[,1], IRanges(1,as.numeric(size[,2])))
sw=slidingWindows(x = cnv_bin0, width = window_w, step = window_step)
cnv_bin=sw@unlistData
  
subject=GRanges(sapply(strsplit(rownames(raw_mat),":|_|-"),'[',1),
                IRanges(as.numeric(sapply(strsplit(rownames(raw_mat),":|_|-"),'[',2))+1,
                        as.numeric(sapply(strsplit(rownames(raw_mat),":|_|-"),'[',3))))
ov=findOverlaps(cnv_bin, subject )
ov=as.matrix(ov)
  
smooth_mat=matrix(ncol=ncol(raw_mat), nrow=length(cnv_bin))
for(ii in 1:length(cnv_bin)){
  ri=colSums(raw_mat[ov[which(ov[,1]==ii),2],, drop=F])
  smooth_mat[ii,]=ri
}
colnames(smooth_mat)=colnames(raw_mat)
rownames(smooth_mat)=paste0(as.character(seqnames(cnv_bin)),'-',start(cnv_bin),'-', end(cnv_bin))
rownames(cell_type)=cell_type[,1]
mat_celltype=cell_type[match(colnames(smooth_mat), rownames(cell_type)),2]
nontumor_mat=smooth_mat[,which((mat_celltype %in% normal_lab))]###
bin_ind=which(rowMedians(nontumor_mat)!=0)
norm_matrix=matrix(rep(rowMedians(nontumor_mat),ncol(smooth_mat)), ncol=ncol(smooth_mat), byrow=F)
smooth_mat_norm=smooth_mat[bin_ind,]/(norm_matrix[bin_ind,])
smooth_mat_norm=apply(smooth_mat_norm, c(1,2), function(x) min(x,5))
plot_matrix=t(smooth_mat_norm)
chrgap=(table(sapply(strsplit(rownames(smooth_mat[bin_ind,]),'-'),'[',1))[paste0('chr',1:nrow(size))])
chrgap[is.na(chrgap)]=0
col_lab=rep(" ", ncol(plot_matrix))
col_lab[c(0, cumsum(chrgap)[which(chrgap!=0)])[1:length(which(chrgap!=0))]+chrgap[which(chrgap!=0)]/2]=names(chrgap)[!is.na(names(chrgap))]
celltype=cell_type
celltype=as.data.frame(cell_type)
rownames(celltype)=celltype[,1]
## break for coloring
plot_matrix=plot_matrix*2
plot_matrix=apply(plot_matrix, c(1,2), function(x) if(x<=2.5 & x>=1.5){x=2}else{x=x})
  
breaklength = 50
setcolor = colorRampPalette(c('#6898C5' , 'white', '#C9676F'))(breaklength)
setbreaks = c(seq(min(plot_matrix), 1.7, length.out=ceiling(breaklength/2) + 1), 
              c(2.3,seq((max(plot_matrix)-2.3)/breaklength+2.3, max(plot_matrix), 
                        length.out=floor(breaklength/2)))[1:(breaklength/2)])
  
  

ann_colors = list(tissue = c(Tumor = '#D95F02', Adj = '#31A354', WT = '#1F78B4'))
tmp=pheatmap::pheatmap(plot_matrix,
                       cluster_cols = F, cluster_rows = TRUE,
                       show_rownames = F,
                       show_colnames = T,
                       color=setcolor,
                       breaks = setbreaks,
                       labels_col=col_lab,
                       clustering_distance_rows = "correlation",
                       clustering_method = "ward.D2",
                       gaps_col=cumsum(chrgap),
                       cutree_rows = nclust,
                       annotation_row=celltype[, ncol(celltype),drop=F],
                       filename = './Alleloscope.pdf',
                       width = 15, 
                       height = 6,
                       annotation_colors = ann_colors)
  



