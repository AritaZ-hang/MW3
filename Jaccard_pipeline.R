
## Do not execute, strictly for test.
if(F)
{
  setwd('~/Desktop/MW3.0/jaccard/')
  source('./cutoff.R')
  source('./Jaccard_plotting.R')
  source('./knee.R')
  dat = fread('./test_jaccard.csv')
  output = './'
}


args = commandArgs(trailingOnly=TRUE)
input = args[1]
output = args[2]
script = args[3]


source(paste0(script, 'knee.R'))
source(paste0(script, 'cutoff.R'))
source(paste0(script, 'Jaccard_plotting.R'))

if(!dir.exists(output))
{dir.create(output)}

suppressPackageStartupMessages(library('data.table'))
suppressPackageStartupMessages(library('ggplot2'))
suppressPackageStartupMessages(library('dplyr'))

dat = fread(input)
dat = dat[order(dat$jaccard,decreasing = T),]

## plot Jaccard plot

dat1=dat[dat$jaccard>0,]

## compute the knee of jaccards.
head(dat)

jaccard_threshold = get_density_threshold_CL(dat$jaccard,type='jaccard')
plotJaccard(df = dat, outD = output, jaccard_threshold = jaccard_threshold[1])

dat2=dat1[dat1$jaccard> jaccard_threshold,]

ovdf=dat2
bead=union(unique(ovdf$match_a),unique(ovdf$match_b)) 
length(bead)

##eat up beads and get translate file
barcode_translate_group=data.frame(matrix(ncol = 2,nrow = 0))
while(dim(ovdf)[1]>0){
  barcode=ovdf$match_a[1]
  friends1 <-NULL
  friends2 <-NULL
  barcode_combine<-NULL
  friendsRow1 <- which(barcode ==  ovdf[,"match_a", drop = TRUE])
  if(length(friendsRow1) > 0){
    friends1 <- ovdf[friendsRow1,]}
  friendsRow2 <- which(barcode ==  ovdf[,"match_b", drop = TRUE])
  if(length(friendsRow2) > 0){
    friends2 <- ovdf[friendsRow2,]}
  if(length(friendsRow2) > 0 ||length(friendsRow1) > 0){
    friends=rbind(friends1,friends2)
    friends=friends[order(friends$jaccard,decreasing=T),]
    threshold=cutoff(friends$jaccard,0.9)
    friends=friends[1:threshold,]
    if(dim(friends)[1]>0)  {
      friends1 <- as.character(friends[friends$match_a==barcode,]$match_b)
      friends2 <- as.character(friends[friends$match_b==barcode,]$match_a)
      barcode_combine <- c(barcode, friends1,friends2)}
  }
  if(length(friendsRow2) == 0 && length(friendsRow1)==0){
    barcode_combine=barcode}
  ovdf <- ovdf[!ovdf$match_a %in% barcode_combine,] 
  ovdf <- ovdf[!ovdf$match_b %in% barcode_combine,]
  if(length(barcode_combine)>1){
    barcode_translate=data.frame(old=barcode_combine,new=barcode_combine[1])
    barcode_translate_group=rbind(barcode_translate_group,barcode_translate)
  }
}

x=as.data.frame(table(barcode_translate_group$new))
table(x$Freq)
#hist(x$Freq,xlim = c(0,10),breaks = 100)

x1=x[x$Freq<=6,]
table(x1$Freq)

## Freq plotting

beads.freq = as.data.frame(table(x1$Freq))
write.csv(beads.freq, file = paste0(output,'beads_freq.csv'), row.names = F)

message('merging beads num:',sum(x1$Freq))
xx=x1$Var1
barcode_translate_group1=barcode_translate_group[barcode_translate_group$new%in%xx,]
write.csv(barcode_translate_group1,file = paste0(output, 'translate.csv'),row.names = F,quote = F)
