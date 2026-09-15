library(dplyr) 
library(magrittr)
library(edgeR)
library(knitr)
library(ggplot2)
library(stringr)
library(ggthemes)
library(ggrepel)
library(ggalt)
library(reshape2)
library(ggsignif)
library(gridExtra)
library(grid)
library(data.table)
library(tidyr)
library(ggnewscale)

#####################
#primary_folder=[primary_folder]
primary_folder="/analysisdata/fantom6/Interactome/single_cell_wallace/Figure/"
SCAFE_path=paste0(primary_folder,"Data_and_code/SCAFE/sc5end/")

#===============================================================================
# running SCAFE
## SCAFE solo
setwd(paste0(SCAFE_path,"solo"))
system("perl 2023070722.batch_run_scafe.F6_neuron_single_cell_solo.pl")
# -> this output CTSS with and without un-encoded G per replicate per cell-type

## SCAFE aggregate
setwd(paste0(SCAFE_path,"aggregate/run_full"))
system("perl 20230711.batch_run_scafe_interactome.scafe_aggregate.sc5end.pl")
# -> this output all the annotation and coordination of of high-confidence TSS clusters and tCREs across the three cell-types
# -> note: tCRE coordinates were only used in a comparison with ATAC-derived aCRE but not in downstream anaylses   

## SCAFE count
setwd(paste0(SCAFE_path,"count"))
system("perl 20230711.batch_run_scafe_interactome.scafe_count.pl")
# -> this script count the CTSS from each single cell to ATAC-derived aCRE regions, limited to only CTSS found inside high-confidence TSS clusters  
# -> same counting was also performed into tCRE regions, but not used in this manuscript

#===============================================================================
# basic summary of the SCAFE output
##% of UMI used in tCRE

path8=paste0(SCAFE_path,"solo/")
files0=list.files(pattern="read_count.txt", path=path8, recursive = T)
names0=sapply(strsplit(files0,"/"),"[",3)
data0=read.delim(paste0(path8,files0[1]),header=T, stringsAsFactors = F, check.names = F)
data0$lib=names0[1]
for (i in 2:length(files0)){
  data00=read.delim(paste0(path8,files0[i]),header=T, stringsAsFactors = F, check.names = F)
  data00$lib=names0[i]
  data0=rbind(data0,data00)}
sum(data0$total_read_num) #892286652
sum(data0$passed_read_num) #882581955
sum(data0$passed_read_num)/sum(data0$total_read_num) #0.9891238

data0$CTSSUMI_total=0
data0$CTSScell_total=0
data0$CTSSUMI_unG=0
data0$CTSScell_unG=0

files=list.files(pattern="collapse.ctss.bed.gz", path=path8, recursive = T)
files=files[-grep("tbi",files)]
filesa=files[grep("unencoded_G",files)]
filesb=files[-grep("unencoded_G",files)]

files.names=sapply(strsplit(filesa, "\\/"),"[", 1)
for (i in 1:length(filesa)){
  data1=read.delim(paste0(path8,filesa[i]),header=F, stringsAsFactors = F, check.names = F)
  data2=read.delim(paste0(path8,filesb[i]),header=F, stringsAsFactors = F, check.names = F)
  data0$CTSSUMI_total[which(data0$lib == files.names[i])]=sum(data2$V5)
  data0$CTSScell_total[which(data0$lib == files.names[i])]=sum(data2$V4)
  data0$CTSSUMI_unG[which(data0$lib == files.names[i])]=sum(data1$V5)
  data0$CTSScell_unG[which(data0$lib == files.names[i])]=sum(data1$V4)}

write.table(data0, gzfile(paste0(SCAFE_path,"SCAFE_read_summary.tsv.gz")), col.names=T, row.names=F, quote=F, sep="\t")
sum(data0$CTSSUMI_total) #841699059 
sum(data0$CTSSUMI_unG) #426907683 
sum(data0$CTSSUMI_unG)/sum(data0$CTSSUMI_total) #0.5071975

cluster=read.delim(paste0(SCAFE_path,"aggregate/run_full/out/annotate/sc5end.iPSC_NSC_Neuron/bed/sc5end.iPSC_NSC_Neuron.cluster.annot.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
sum(cluster$V5) #575969385
# number of UMI in cluster = 575969385/841699059 = 0.6842937
