## adjust functions in CopySCAT

## 1. normalizeMatrixN
normalizeMatrixN <- function(inputMatrix,logNorm=FALSE,maxZero=2000,imputeZeros=FALSE,blacklistProp=0.8, priorCount=1,blacklistCutoff=100,dividingFactor=1e6,upperFilterQuantile=0.95)
{
  sc_t<-data.table(t(inputMatrix))
  #sc_t
  cellReads<-data.table::transpose(sc_t[,lapply(.SD,sum)],keep.names="Cell")
  pdf(str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"_signal_distribution.pdf"),width=6,height=4)
  hist(cellReads$V1, breaks=50,main=scCNVCaller$outPrefix,xlab = "Signal")
  abline(v=quantile(cellReads$V1,upperFilterQuantile),col=c("red"),lty=2)
  dev.off()
  
  readsCells=cellReads[,mean(V1)]
  #apply quantile filter
  sc_t[,cellReads[cellReads[,V1>quantile(cellReads$V1,upperFilterQuantile)]]$Cell:=NULL]
  nCells<-nrow(cellReads)
  print(str_c("Total number of starting cells: ",nCells," Average reads per cell: ",mean(readsCells)))
  scCNVCaller$meanReadsPerCell<-readsCells
  blacklistPropCutoff=blacklistProp*nCells
  #good
  #find bad columns
  sc_lines<-data.table(inputMatrix)
  # head(sc_lines)
  #blacklistCutoff = 500
  sc_lines<-sc_lines[,lapply(.SD,function(x) x<blacklistCutoff)][,lapply(.SD,sum)]
  sc_pos<-data.table::transpose(sc_lines,keep.names = "Pos")
  blacklistRegions<-sc_pos[which(sc_pos[,V1>=blacklistPropCutoff]),]$Pos
  print(sprintf("Blaclist regions: %d",length(blacklistRegions)))
  print(sprintf("Total bins: %d",nrow(sc_pos)))
  if (length(blacklistRegions)==nrow(sc_pos))
  {
    print("Error: no regions meet cutoff criteria")
    return(NULL)
  }
  sc_t2<-sc_t[which(sc_pos[,V1<blacklistPropCutoff]),lapply(.SD,function(x) as.numeric(x<blacklistCutoff))][,lapply(.SD,sum)]
  print(sc_t[which(sc_pos[,V1<blacklistPropCutoff]),lapply(.SD,function(x) as.numeric(x<blacklistCutoff))][1:10,1:10])
  #print(sc_t2[,1:100])
  
  sc_zeros<-data.table::transpose(sc_t2,keep.names="Cells")
  #print(sc_zeros)
  #zero_cutoff=zero_cutoff
  zero_list<-sc_zeros[which(sc_zeros[,V1<maxZero])]$Cells
  zero_list
  #print(zero_list)
  #zeros<-rownames_to_column(tmp1,var="Chrom") %>% dplyr::filter(!(Chrom %in% blacklistRegions$Chrom)) %>% summarise_if(is.numeric,funs(sum(.==0))) %>% gather(Cell,Value) %>% mutate(cellPass=(Value<maxZero))
  scCNVCaller$cellsPassingFilter<-length(zero_list)
  print(str_c(scCNVCaller$cellsPassingFilter, " cells passing filter"))
  #get low confidence regions
  #message(blacklistRegions$Chrom)
  #hist(zeros$Value,breaks=100)
  #zeros$Cell[zeros$cellPass==TRUE]
  tmp1 <- as.data.frame(sc_t,stringsAsFactors=FALSE) %>% dplyr::select(zero_list)
  raw_medians<-t(data.table::transpose(sc_t[,..zero_list])[,lapply(.SD,median)])
  #tmp1
  tmp1<-cbind.data.frame(tmp1,raw_medians,stringsAsFactors=FALSE)
  #impute zeros in cells passing filter
  if (imputeZeros==TRUE)
  {
    tmp3 <- tmp1 %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),~ if_else(. == 0,raw_medians,.)) %>% dplyr::select(-raw_medians)
  }
  else
  {
    tmp3 <- tmp1  # %>% dplyr::select(-raw_medians)
  }
  #now normalize quantiles to account for differences in coverage
  #scData_n<- normalize.quantiles(scData_t)
  scData_n<-cpm(tmp3,log = logNorm,prior.count = priorCount)
  
  rownames(scData_n)<-colnames(inputMatrix)
  colnames(scData_n)<-colnames(tmp3)
  #scData_n
  head(scData_n)
  #var_int<-apply(scData_n[1:nrow(scData_n),], 1, var, na.rm=TRUE)
  #scData_nv<-cbind(scData_n,nvariance=var_int)
  scData_nc_split <- rownames_to_column(as.data.frame(scData_n),var = "Loc") %>% mutate(blacklist=(Loc %in% blacklistRegions)) %>% separate(col=Loc,into=c("chrom","pos"))
  scData_k<- scData_nc_split %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),list(~ (./ dividingFactor)))
  return(scData_k)
}

## 2.collapseChrom3N
collapseChrom3N<-function(inputMatrix,minimumSegments=40,summaryFunction=cutAverage,logTrans=FALSE,binExpand=1,minimumChromValue=2,tssEnrich=5,logBase=2,minCPG=300,powVal=0.73)
{
  #remove centromeric bins
  #these are ordered in the same fashion so should be accurate
  #divide signal by cpg density + 1 (To avoid dividing by zero)
  #remove blacklist regions
  # %>% dplyr::filter(cpg!=0)
  scData_k_norm <- inputMatrix %>% mutate(chromArm=scCNVCaller$cytoband_data$V4,cpg=scCNVCaller$cpg_data$cpg_density) %>% mutate(chrom=str_c(chrom,chromArm)) %>% dplyr::filter(cpg>minCPG)
  if (binExpand>1)
  {
    #collapse bins
    scData_k_norm <- scaleMatrixBins(scData_k_norm,binSizeRatio = binExpand,"pos")
    #rename column
    scData_k_norm <- scData_k_norm %>% rename(pos = pos_b)
  }
  # %>% dplyr::filter(cpg!=0)
  #gender determination
  sckn<-data.table(scData_k_norm)
  total_cpg<-sckn[,summaryFunction(cpg),by="chrom"]
  sckn[,c("blacklist","raw_medians","chromArm","cpg","pos"):=NULL]
  tail(colnames(sckn))
  setkey(sckn,chrom)
  chromXName="chrXq"
  chromYName="chrYq"
  x_tmp<-(tail(scCNVCaller$cytoband_data[which(scCNVCaller$cytoband_data$V1=="chrX"),c(1,4)],n=1))
  chromXName<-paste(x_tmp$V1,x_tmp$V4,sep="")
  y_tmp<-(tail(scCNVCaller$cytoband_data[which(scCNVCaller$cytoband_data$V1=="chrY"),c(1,4)],n=1))
  chromYName<-paste(y_tmp$V1,y_tmp$V4,sep="")
  ## check if the chrom name is in accordance with the input data;
  if(!chromXName %in% unique(scData_k_norm[,1]))
    chromXName=tail(unique(scData_k_norm[,1])[stringr::str_detect(unique(scData_k_norm[,1]), pattern='chrX')])[length(tail(unique(scData_k_norm[,1])[stringr::str_detect(unique(scData_k_norm[,1]), pattern='chrX')]))]
  if(!chromYName %in% unique(scData_k_norm[,1]))
    chromYName=tail(unique(scData_k_norm[,1])[stringr::str_detect(unique(scData_k_norm[,1]), pattern='chrY')])[length(tail(unique(scData_k_norm[,1])[stringr::str_detect(unique(scData_k_norm[,1]), pattern='chrY')]))]
  
  sckn<-sckn[c(chromXName,chromYName),lapply(.SD,quantile,probs=0.8),by="chrom"]
  nrow(sckn)
  total_cpg[chrom %in% c(chromXName,chromYName)]
  sckn[,cpg:=total_cpg[chrom %in% c(chromXName,chromYName)]$V1]
  #sckn[,print(.SD)]
  xy_signal<-data.table::transpose(sckn[,lapply(.SD,"/",1+cpg),by="chrom"][,cpg:=NULL],make.names = "chrom")[,lapply(.SD,quantile,0.8)]
  
  #t(sckn)[2:ncol(sckn),]
  #xy_signal <- as_tibble(apply(t(sckn)[2:ncol(sckn),],2,as.numeric))
  #xy_signal
  #xy_signal <-scData_k_norm %>% dplyr::filter(chrom %in% c("chrXq","chrYq")) %>% group_by(chrom) %>% summarise_if(is.numeric,list(quantile),probs=0.8) 
  sexCutoff=5e-7
  logTrans=FALSE
  if (logTrans)
  {
    xy_signal<-data.table::transpose(sckn[,lapply(.SD,"/",1+log(cpg,base=2)),by="chrom"][,cpg:=NULL],make.names = "chrom")[,lapply(.SD,quantile,0.8)]
    
    #xy_signal <- xy_signal %>% rowwise() %>% mutate(totalcpg=summaryFunction(cpg)) %>% dplyr::select(-cpg) %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),funs(./(1+log(totalcpg,base=2)))) 
  } else
  {
    xy_signal<-data.table::transpose(sckn[,lapply(.SD,"/",1+cpg),by="chrom"][,cpg:=NULL],make.names = "chrom")[,lapply(.SD,quantile,0.8)]
    
    #xy_signal <- xy_signal %>% rowwise() %>% mutate(totalcpg=summaryFunction(cpg)) %>% dplyr::select(-cpg) %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),funs(./(1+totalcpg))) 
    sexCutoff=0.25
  }
  #print(xy_signal)
  diffXY<-xy_signal$chrXq-xy_signal$chrYq
  # print(abs(diffXY))
  
  if (abs(diffXY)<sexCutoff)
  {
    scCNVCaller$isMale=TRUE
  }
  #top part works awesome
  
  #now to convert rest
  sckn<-data.table(scData_k_norm)
  #print(nrow(sckn))
  sckn <- sckn[chromArm!="cen"][blacklist==FALSE]
  #print(nrow(sckn))
  #total_cpg <- sckn[,summaryFunction(cpg),by="chrom"]
  #total_cpg
  #setnames(total_cpg,"V1","totalcpg")
  #remove raw_medians,blacklist,chromArm,cpg
  sckn<-sckn[,c("blacklist","raw_medians","chromArm","pos"):=NULL]
  #mm I see the issue - the cutAverage upper quantile busts / overcorrects the cpg - use a median instead
  cpg_density<-sckn[,cpg]
  median_chromosome_density<-sckn[1:(nrow(sckn)),lapply(.SD,summaryFunction),by="chrom"]
  med_temp<-data.table::transpose(median_chromosome_density[,c("cpg"):=NULL],make.names="chrom")
  #head(med_temp)
  data.table::transpose(med_temp[,lapply(.SD,summaryFunction)],keep.names="chrom")
  median_chrom_signal<-data.table::transpose(med_temp[,lapply(.SD,summaryFunction)],keep.names="chrom")
  #print(median_chrom_signal)
  
  
  #median_chromosome_density[,1:5]
  #]
  #median_chromosome_density<-scData_k_norm %>%
  #  group_by(chrom) %>% summarise_if(is.numeric,list(summaryFunction)) 
  #message(nrow(scData_k_norm))
  #now append 
  #summaryFunction
  
  
  
  #TODO: cpg of 200 and medianNorm cutoffs only work for 1e6 bins; need to change for others
  if (logTrans==TRUE)
  {
    #5.1e-6 old 
    #5.1e-06, totalcpg 130 for 200k
    #THIS ONE
    #transpose(sckn[,lapply(.SD,"/",1+log(cpg,base=2)),by="chrom"][,cpg:=NULL]
    sckn<-sckn[,lapply(.SD,"/",1+log(tssEnrich*cpg^powVal,base=logBase)),by="chrom"][,cpg:=NULL][,lapply(.SD,summaryFunction),by="chrom"]
    #scData_prechrom<-scData_k_norm %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),funs(./(1+log(cpg,base=2)))) %>% group_by(chrom) %>% summarise_if(is.numeric,list(summaryFunction)) %>% 
    #  dplyr::select(-cpg) %>% dplyr::filter(chrom %in% median_chrom_signal$chrom) %>% 
    #  mutate(medianNorm=median_chrom_signal$Value) %>% dplyr::filter(medianNorm>minimumChromValue) 
    #+log(totalcpg,base=2)
    #was 200
  }
  else
  {
    #THIS ONE
    sckn<-sckn[,lapply(.SD,"/",1+tssEnrich*cpg^powVal),by="chrom"][,cpg:=NULL][,lapply(.SD,summaryFunction),by="chrom"]
    #sckn<-sckn[,lapply(.SD,"/",1+tssEnrich*log(cpg,base=logBase)),by="chrom"][,cpg:=NULL][,lapply(.SD,sum),by="chrom"]
    #scData_prechrom<-scData_k_norm %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),funs(./(1+sqrt(cpg)))) %>% group_by(chrom) %>% summarise_if(is.numeric,list(summaryFunction)) %>% 
    #  dplyr::select(-cpg) %>% dplyr::filter(chrom %in% median_chrom_signal$chrom) %>%
    #  mutate(medianNorm=median_chrom_signal$Value) %>% dplyr::filter(medianNorm>minimumChromValue) 
  }
  
  sckn<-sckn[,medianNorm:=median_chrom_signal$V1]
  sckn<-sckn[which(sckn[,medianNorm>minimumChromValue])]
  #return the appropriate matrices
  return(as.data.frame(sckn))
}

## 3.annotateCNV4
annotateCNV4 <- function(cnvResults,saveOutput=TRUE,maxClust2=4,outputSuffix="_1",sdCNV=0.6,filterResults=TRUE,filterRange=0.8,minAlteredCells=40)
{
  ## cnvResults = candidate_cnvs_clean
  ## outputSuffix='_1'
  ## sdCNV = 0.9
  ## filterResults = FALSE
  ## maxClust2 = 10
  
  cell_assignments<-cnvResults[[1]]
  cell_assignments<-column_to_rownames(cell_assignments,var = "Cells")
  chrom_clusters_final<-cnvResults[[2]]
  #colnames(chrom_clusters_final)=c("Chrom",seq(1:(maxClust2-2)),"0")
  shift_val=0
  if (length(which(chrom_clusters_final$V2==0))>1)
  {
    shift_val = mean(chrom_clusters_final$V1[which(chrom_clusters_final$V2==0)])
  }
  #print(shift_val)
  #identify possible CNVs - can use bins to test
  possCNV <- cell_assignments %>% summarise_if(is.numeric,list(max)) %>% gather(Chrom,maxClust) %>% dplyr::filter(maxClust>1)
  possCNV <- possCNV %>% mutate(Type="Unknown")
  # print(cell_assignments %>% summarise_if(is.numeric,list(max)))
  #%>% gather(Chrom,maxClust))
  # print(possCNV)
  #label by Z scores
  chrom_clusters_final<-chrom_clusters_final %>% gather("cluster","zscore",starts_with("V")) %>% mutate(cluster=str_remove(cluster,"V")) 
  #print(chrom_clusters_final)
  #print(which(is.nan(chrom_clusters_final$zscore)))
  chrom_clusters_final$zscore<-sapply(sapply(chrom_clusters_final$zscore-shift_val,pnorm,log.p=TRUE),qnorm,2,sdCNV,log.p=TRUE)
  #replace each column
  #print(chrom_clusters_final)
  trimmedCNV<-vector()
  thresholdVal=filterRange
  for (chrom in possCNV$Chrom)
  {
    #  print(chrom)
    zscores<-chrom_clusters_final %>% dplyr::filter(Chrom==chrom) %>% dplyr::select(zscore)
    #  print(zscores)
    #print(head(cell_assignments[chrom]))
    #print(as.character(factor(cell_assignments[,chrom],levels=seq(from=0,to=6),labels = c("2",zscores$zscore))))
    cell_assignments[,chrom]<-as.numeric(as.character(factor(cell_assignments[,chrom],levels=seq(from=0,to=6),labels = c("2",zscores$zscore))))
    # print(zscores)
    #  print(chrom)
    #    print(range(zscores))
    if (diff(range(zscores))>=thresholdVal)
    {
      #     print(chrom)
      trimmedCNV<-append(trimmedCNV,chrom)
    }
  }
  #print(trimmedCNV)
  #output CNV list by cluster
  #cell_assignments %>% dplyr::filter()
  #WITHOUT BLANKS
  #consensus_CNV_clusters<-rownames_to_column(cell_assignments) %>% gather(Chrom,clust,2:(ncol(cell_assignments)+1)) %>% dplyr::filter(Chrom %in% possCNV$Chrom) %>% spread(Chrom,clust) %>% arrange(rowname)
  
  #WITH BLANKS
  if (!filterResults)
  {
    consensus_CNV_clusters<-rownames_to_column(cell_assignments) %>% dplyr::filter(str_detect(rowname,"X",negate=TRUE)) %>% gather(Chrom,clust,2:(ncol(cell_assignments)+1)) %>% dplyr::filter(Chrom %in% possCNV$Chrom) %>% spread(Chrom,clust) %>% arrange(rowname)
  }
  else
  {
    #identify chromosomes not used
    unused = cnvResults[[1]] %>% dplyr::filter(str_detect(Cells,"X",negate=TRUE)) %>% gather(chrom,cluster,starts_with("chr")) %>% group_by(chrom,cluster) %>% dplyr::filter(cluster!=0) %>% dplyr::summarise(num=n()) %>% group_by(chrom) %>% summarise(min=min(num)) %>% dplyr::filter(min<minAlteredCells)
    # print(unused)
    #print(trimmedCNV)
    #print(trimmedCNV %in% unused$chrom)
    trimmedCNV<-trimmedCNV[!(trimmedCNV %in% unused$chrom)]
    consensus_CNV_clusters<-rownames_to_column(cell_assignments) %>% dplyr::filter(str_detect(rowname,"X",negate=TRUE)) %>% gather(Chrom,clust,2:(ncol(cell_assignments)+1)) %>% dplyr::filter(Chrom %in% trimmedCNV) %>% spread(Chrom,clust) %>% arrange(rowname)
  }
  if (saveOutput==TRUE)
  {
    write.table(x=consensus_CNV_clusters,file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,outputSuffix,"_cnv_scores.csv"),quote=FALSE,row.names = FALSE,sep=",")
  }
  return(list(chrom_clusters_final,possCNV,consensus_CNV_clusters))
}

## 4.identifyNonNeoplastic
identifyNonNeoplastic <- function(inputMatrix,estimatedCellularity=0.8,nmfComponents=5,outputHeatmap=TRUE,cutHeight=0.6,methodHclust="ward.D")
{
  #uses NMF and fast hclust packages
  message("Running NMF")
  res <- nmf(column_to_rownames(inputMatrix,var="chrom"), c(nmfComponents), "lee", seed="nndsvd",.stop=nmf.stop.threshold(0.1),maxIter=2500)
  message("Computing clusters")
  #dist<-dist(t(coef(res)),method="euclidean")
  tscale<-scale(x=t(coef(res)),center=TRUE)
  dist=dist(tscale,method="euclidean")
  
  clusts<-fastcluster::hclust(dist,method = methodHclust)
  #average or ward
  #clusts<-agnes(t(coef(res)),diss = FALSE,metric="euclidean",method="ward")
  #YEAH THIS WORKS AWESOME
  #some tricksy samples may need 4
  
  #add PDF
  cell_assigns<-cutree(clusts,h=max(clusts$height)*cutHeight)
  
  if (outputHeatmap==TRUE)
  {
    pdf(file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"_nmf_heatmap.pdf"),width=6,height=6)
    #print(cell_assigns)
    #print(factor(cell_assigns))
    heatmap.2(coef(res),Rowv=FALSE,Colv=as.dendrogram(clusts),dendrogram="column",density.info="none",trace="none",scale="none",labCol=FALSE,col=colorRampPalette(viridis(5)),symkey=FALSE,useRaster=TRUE,ColSideColors=viridis(n=length(unique(as.character(cell_assigns))))[cell_assigns])
    legend("topright",fill=viridis(n=3),x.intersp = 0.8,y.intersp=0.8,legend=unique(as.character(cell_assigns)),horiz = FALSE,cex = 0.9,border=TRUE,bty="n")
    dev.off()
  }
  by_cluster<-left_join(rownames_to_column(data.frame(t(column_to_rownames(inputMatrix,"chrom"))),var="barcode"),
                        rownames_to_column(data.frame(cluster=cell_assigns),var="barcode"),by="barcode")
  #now estimate the residuals etc for each cluster
  if (outputHeatmap==TRUE)
  {
    pdf(file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"_nmf_violinplot.pdf"),width=12,height=6)
    print(ggplot(by_cluster %>% gather(starts_with("chr"),key="chrom",value="value"),aes(chrom,value)) + geom_violin() + theme(axis.text.x=element_text(angle=-90,vjust = 0.5,hjust = 0, color = "#000000",size = 6)) + facet_wrap(~cluster))
    dev.off()
  }
  target_cluster_1=data.frame(t(column_to_rownames(by_cluster %>% group_by(cluster) %>% summarise_if(is.numeric,list(mean)),var="cluster"))) %>% summarise_if(is.numeric,list(var)) %>% gather(key="cluster",value="var") %>% mutate(cluster=as.numeric(str_remove(cluster,"X"))) %>% arrange(var) %>% slice_head(n=1) %>% dplyr::select(cluster)
  
  #cluster_order<-data.frame(cluster=cell_assigns) %>% group_by(cluster) %>% summarise(number=n()) %>% arrange(number)
  tumor_cell_ids=names(which(cell_assigns!=target_cluster_1))
  normalCluster=unlist(target_cluster_1)
  normal_cell_ids=names(which(cell_assigns==unlist(normalCluster)))
  return(list(cellAssigns=cell_assigns,normalBarcodes=normal_cell_ids,clusterNormal=normalCluster))
}

## 5. identifyDoubleMinutes
identifyDoubleMinutes<-function(inputMatrix,minCells=100,qualityCutoff2=100,peakCutoff=5,lossCutoff=-1,doPlots=FALSE,imageNumber=1000,logTrans=FALSE,cpgTransform=FALSE,doLosses=FALSE,minThreshold=4)
{
  dm_per_cell<-data.frame(cellName=colnames(inputMatrix %>% dplyr::select(ends_with(scCNVCaller$cellSuffix))),stringsAsFactors=FALSE)
  #  print(head(dm_per_cell))
  dm_per_cell[,seq(from=2,to=5)]<-FALSE
  
  #initialize repeat index
  firstRepeatIndex=2
  #MINIMUM # of cells with alteration - 50 is defualt
  qualityCutoff=minCells
  #test - can remove centromeres
  #cytoband_data$V4
  scData_k_cpg<-inputMatrix %>% mutate(cpg=scCNVCaller$cpg_data$cpg_density) %>% mutate(arm=scCNVCaller$cytoband_data$V4) %>% dplyr::filter(arm!="cen", blacklist==0) %>% dplyr::select(-arm,-blacklist) #%>%  #cpg+
  if (cpgTransform==TRUE)
  {
    if (logTrans==FALSE)
    {
      scData_k_cpg <- scData_k_cpg %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),funs(./(cpg+1)))
    }
    else
    {
      scData_k_cpg <- scData_k_cpg %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),funs(./(log(cpg,base=2)+1)))
    }
  }
  for (m in 1:(nrow(scCNVCaller$chrom_sizes[scCNVCaller$chrom_sizes[,3],])-3)) # so we adjust the algorithm only searching on regular chromosomes.
  {
    
    #loop through chromosomes
    currentChrom = scCNVCaller$chrom_sizes$chrom[m]
    message(str_c("processing ",currentChrom))
    # print(nrow(scData_k_cpg))
    #scData_1h <- scData_k_cpg %>% dplyr::filter(chrom==currentChrom ) %>% mutate(loc1=paste(chrom,pos,sep="_")) %>% dplyr::select(-chrom,-pos)
    #NEW:
    scData_1h <- data.frame(scData_k_cpg %>% dplyr::filter(chrom==currentChrom ) %>% mutate(loc1=paste(chrom,pos,sep="_")) %>% dplyr::select(-chrom,-pos,-cpg))
    #  print(scData_1h)
    #seleect peak cutoff
    #unlistValues <- unlist((scData_1h %>% dplyr::select(-loc1)))
    #peakCutoff = quantile(scale(unlistValues,center=median(unlistValues),scale=IQR(unlistValues)),.99)
    
    #scData_1h[,2:ncol(scData_1h)]<-scale(scData_1h[,2:ncol(scData_1h)],center=TRUE,scale=TRUE)
    dm_list=data.frame(dm="null",stringsAsFactors=FALSE)
    if (m>1)
    {
      #rm(dm_per_cell_clean)
      rm(recurrentEvents)
    }
    #initialize temporary matrix
    dm_per_cell_temp<-as.data.frame(matrix(nrow = nrow(dm_per_cell),ncol=5))
    colnames(dm_per_cell_temp)<-c("cellName","V2","V3","V4","V5")
    dm_per_cell_temp[,1]<-dm_per_cell$cellName
    #blank out temporary file
    #need to be smarter with this bit here but I can't be bothered RN
    dm_per_cell_temp[,seq(from=2,to=5)]<-FALSE
    
    #START TESTING HERE
    #temp_dm<-getDoubleMinutes(scData_1h,17,doPlot=TRUE)
    #dm_index<-which(dm_list$dm==temp_dm_boundaries)
    #dm_index
    plottableCell=FALSE
    plotEachTime=imageNumber
    pb<-txtProgressBar(min=0,max=nrow(dm_per_cell)-1,style=3)
    for (cellNum in 1:(nrow(dm_per_cell)-1))
    {
      setTxtProgressBar(pb,cellNum)
      #  print(cellNum)
      # medianDensity<-median(scData_1h[,5])
      # if (medianDensity>0){
      #   tmp<-scData_1h[,5]
      #   #replace telomere values
      #   tmp[tmp>4e-6]<-medianDensity
      #   plot(tmp)
      #   #consider using IQR Here
      #   tmp<-as.numeric(scale(tmp,center=median(tmp),scale=IQR(tmp)),stringsAsFactors=FALSE)
      if ((cellNum %% plotEachTime)==0)
      {
        plottableCell=TRUE
        #  print(cellNum)
        #  print(plottableCell)
      }
      if (doPlots & plottableCell)
      {
        #print(plottableCell)
        pdf(str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"_cell",cellNum,"_",currentChrom,".pdf"),width=6,height=4)
        temp_dm<-getDoubleMinutes(scData_1h,cellNum,doLosses=doLosses,peakCutoff = peakCutoff,lossCutoff = lossCutoff,doPlot=TRUE,minThreshold=minThreshold)
        dev.off()
        plottableCell=FALSE
      }
      else
      {
        temp_dm<-getDoubleMinutes(scData_1h,cellNum,doLosses=doLosses,peakCutoff = peakCutoff,lossCutoff = lossCutoff,doPlot=FALSE,minThreshold=minThreshold)
      }
      
      if (!is.null(temp_dm))
      {
        #iterate through each DM identified
        for (j in 1:length(temp_dm))
        {
          #LENGTH CHECK on temp_dm
          #get boundaries - this will be our identifier for this location
          #OBSOLETE: temp_dm_boundaries<-sprintf(fmt="%s.%s",temp_dm[j,1],temp_dm[j,ncol(temp_dm)])
          #this works to ID if this vector already exists
          #,arr.ind=TRUE
          dm_index<-which(dm_list==temp_dm[j])
          #message(dm_index)
          #message("identify index")
          if (length(dm_index)==0)
          {
            #message("adding to list")
            dm_list<-rbind.data.frame(dm_list,data.frame(dm=temp_dm[j],stringsAsFactors = FALSE),stringsAsFactors=FALSE)
            dm_per_cell_temp[cellNum,nrow(dm_list)]=TRUE
            #TODO: get way to name this after the interval (name column name)
            # dm_per_cell<-cbind.data.fraame(dm_per_cell,data.frame(temp_dm_boundaries=rep(FALSE,times=nrow(dm_per_cell))))
          } else {
            
            # dm_index2<-which(dm_list==temp_dm_boundaries,arr.ind=TRUE)
            #because cellName is a column, and first value of dm is always blank, so don't add one
            dm_per_cell_temp[cellNum,dm_index]=TRUE
          }
        }
      }
    }
    close(pb)
    names(dm_per_cell_temp)[2:nrow(dm_list)]<-dm_list$dm[2:nrow(dm_list)]
    #dm_per_cell
    #dm_per_cell %>% rename_at(2:(nrow(dm_list)-1),~ dm_list$dm)
    #test for overlaps - start index 2
    #start one smaller
    
    dm_list_test2=dm_list %>% mutate(newIndex=1)
    
    chrom_name=currentChrom
    #now loop through all intervals - if detected greater than one
    if (nrow(dm_list_test2)>1)
    {
      for (j1 in 2:(nrow(dm_list_test2)))
      {
        if(dm_list_test2$newIndex[j1]==1)
        {
          dm_list_test2$newIndex[j1]=j1
        }
        #(j1+1)
        for (j2 in 2:(nrow(dm_list_test2)))
        {
          if (dm_list_test2$newIndex[j2]!=j1)
          {
            #message(j2)
            interval.1<-dm_list_test2$dm[j1]
            interval.2<-dm_list_test2$dm[j2]
            overlapResult<-testOverlap(interval.1,interval.2)
            #message(overlapResult[3])
            if (overlapResult[3]=="TRUE")
            {
              #update table
              dm_list_test2$dm[j1]=str_c(chrom_name,"_",overlapResult[1],".",chrom_name,"_",overlapResult[2])
              #message(dm_list_test2$dm[j2])
              dm_list_test2$dm[j2]=str_c(chrom_name,"_",overlapResult[1],".",chrom_name,"_",overlapResult[2])
              #message(dm_list_test2$dm[j2])
              dm_list_test2$newIndex[j2]=dm_list_test2$newIndex[j1]
            }
          }
        }
      }
      #NOTE: still some minor kinks in the overlap function; it works but the intervals aren't coming up perfect
      #TODO: fudge indices to empty blanks
      tmp_fact<-factor(dm_list_test2$newIndex)
      levels(tmp_fact)<-seq(1,length(levels(tmp_fact)))
      dm_list_test2$newIndex<-as.numeric(levels(tmp_fact)[tmp_fact])
      #option: can clean up dataframe here and replace NA with FALSE, and then proceed instead of making a giant one up front
      dm_per_cell_temp[is.na(dm_per_cell_temp)]<-FALSE
      #OK that worked, now to amalgamate cell-wide data
      for (j3 in 2:nrow(dm_list_test2))
      {
        #combine and rename columns appropriately
        dm_per_cell_temp[,dm_list_test2$newIndex[j3]]<-(dm_per_cell_temp[,j3] | dm_per_cell_temp[,dm_list_test2$newIndex[j3]])
      }
      #now rename colnames
      for (j3 in 2:nrow(dm_list_test2))
      {
        colnames(dm_per_cell_temp)[dm_list_test2$newIndex[j3]]<-dm_list_test2$dm[j3]
      }
      dm_per_cell_temp<-dm_per_cell_temp[,1:(max(dm_list_test2$newIndex))]
      
      #now summarise and count each; remove any with 5 or fewer events
      
      recurrentEvents <- dm_per_cell_temp  %>% summarize_if(is.logical,sum) %>% select_if(function(x) any(x>qualityCutoff))
      #select only these columns that meet cutoff
      #check that recurrent events have been identified
      if (dim(recurrentEvents)[2]!=0)
      {
        message("copying recurrent events")
        dm_per_cell_clean <- dm_per_cell_temp %>% dplyr::select(cellName,colnames(recurrentEvents))
        #copy over columns and column names and increment column counter - ncol-2
        dm_per_cell[,firstRepeatIndex:(firstRepeatIndex+ncol(dm_per_cell_clean)-2)]=dm_per_cell_clean[,2:(ncol(dm_per_cell_clean))]
        #check indexes in following line
        colnames(dm_per_cell)[firstRepeatIndex:(firstRepeatIndex+ncol(dm_per_cell_clean)-2)]<-colnames(dm_per_cell_clean[2:(ncol(dm_per_cell_clean))])
        firstRepeatIndex=firstRepeatIndex + ncol(dm_per_cell_clean)-1
      }
    }
    #beautiful
  }
  #now clean this
  if (firstRepeatIndex>2)
  {
    dm_per_cell_clean<-dm_per_cell[,1:(firstRepeatIndex-1)]
    #  dm_per_cell_clean %>% summarize_if(is.logical,sum)
    
    #can cutoff here too (def 50-75)
    secondThreshold=qualityCutoff2
    keepAlterations <- dm_per_cell_clean %>% summarize_if(is.logical,sum) %>% gather(Chrom,Value) %>% dplyr::filter(Value>secondThreshold)
    dm_per_cell_clean <- dm_per_cell_clean %>% dplyr::select(cellName, keepAlterations$Chrom)
    return(dm_per_cell_clean)
  }
  else
  {
    print("No double minutes identified")
  }
  return(NULL)
}