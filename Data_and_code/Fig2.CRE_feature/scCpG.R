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

###color#####
library("ggsci")
mypal = pal_npg("nrc", alpha = 1)(9)
mypal
library("scales")
show_col(mypal)

#####################
primary_folder="/analysisdata/fantom6/Interactome/single_cell_wallace/Figure/"
CRE_bed=paste0(primary_folder,"Data_and_code/Fig2.CRE_feature/CRE_bed/")
CpG_folder=paste0(primary_folder,"Data_and_code/Fig2.CRE_feature/CRE_bed/CpG_result/")


#CpG island
#===================================================
cd /analysisdata/fantom6/Interactome/resources/UCSC
wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cpgIslandExt.txt.gz \
| gunzip -c \
| awk 'BEGIN{ OFS="\t"; }{ print $2, $3, $4, $5$6, $7, $8, $9, $10, $11, $12 }' \
| sort-bed - \
> cpgIslandExt.hg38.bed
# This file is placed in [primary_folder]/Data_and_code/Fig2.CRE_feature/CRE_bed/

#===============================================================================
ump=read.delim(paste0(CRE_bed,"untranscribed_aCRE.midpoint_summit5000.bed.gz", header=F, stringsAsFactors = F)
tmp=read.delim(paste0(CRE_bed,"transcribed_aCRE.midpoint_summit5000.bed.gz", header=F, stringsAsFactors = F)
tt=read.delim(paste0(CRE_bed,"transcribed_aCRE.t_summit5000.stranded.bed.gz", header=F, stringsAsFactors = F)
ump$bimodal="untranscrib_midpoint"
tmp$bimodal="transcrib_midpoint"
tt$bimodal="transcrib_stranded_summit"
all=rbind(ump,tmp,tt)
all%>%group_by(bimodal)%>%dplyr::summarise(count=n())

#===============================================================================
#bedtools
setwd(CRE_bed)
system("bedtools intersect -wb -a untranscribed_aCRE.midpoint_summit5000.bed.gz -b cpgIslandExt.hg38.main_chr.bed | gzip > untranscribed_aCRE.midpoint_summit5000.CpG.bed.gz")
system("bedtools intersect -wb -a transcribed_aCRE.t_summit5000.stranded.bed.gz -b cpgIslandExt.hg38.main_chr.bed | gzip > transcribed_aCRE.t_summit5000.stranded.CpG.bed.gz")
system("bedtools intersect -wb -a transcribed_aCRE.midpoint_summit5000.bed.gz -b cpgIslandExt.hg38.main_chr.bed | gzip > transcribed_aCRE.midpoint_summit5000.CpG.bed.gz")

#===============================================================================
setwd(CRE_bed)
CPGump=read.delim("untranscribed_aCRE.midpoint_summit5000.CpG.bed.gz", header=F, stringsAsFactors = F)
CPGtmp=read.delim("transcribed_aCRE.midpoint_summit5000.CpG.bed.gz", header=F, stringsAsFactors = F)
CPGtt=read.delim("transcribed_aCRE.t_summit5000.stranded.CpG.bed.gz", header=F, stringsAsFactors = F)
length(unique(CPGump$V4))#34329
length(unique(CPGtmp$V4))#22210
length(unique(CPGtt$V4))#22213

CPGump$bimodal="untranscrib_midpoint"
CPGtmp$bimodal="transcrib_midpoint"
CPGtt$bimodal="transcrib_stranded_summit"
CPG=rbind(CPGump,CPGtmp,CPGtt)

CPG=left_join(CPG, all[,c(2,4,7)], by=c("V4","bimodal"), copy=F, suffix=c("","_ori"))
CPG$locS=CPG$V2-CPG$V2_ori-5000
CPG$locE=CPG$V3-CPG$V2_ori-5000
CPG1=CPG[which(CPG$V6 != "-"),c(1,19,20,4,10,12,16,17)]
CPG2=CPG[which(CPG$V6 == "-"),c(1,20,19,4,10,12,16,17)]

CPG2$locE=CPG2$locE * (-1)
CPG2$locS=CPG2$locS * (-1)
CPG2$V1="chr1"
CPG1$V1="chr1"
colnames(CPG2)=colnames(CPG1)
CPG1$locS=CPG1$locS+5000
CPG1$locE=CPG1$locE+5000
CPG2$locS=CPG2$locS+5001
CPG2$locE=CPG2$locE+5001

CPG=rbind(CPG1,CPG2)
colnames(CPG)[c(4,5,6,7)]=c("CREID","name","cpgNum","obsExp")
write.table(CPG[order(CPG$locS),],gzfile(paste0(CpG_folder,"ll_summit_5kb_extend.CpG.fakebed.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

CREanno=read.delim(paste0(primary_folder,"scATAC/merged_peak/all_peaks.merged.all.markalone_k16.p.e.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
np.t=CREanno$V4[which(CREanno$promoter_type2 == "promoter-like" & CREanno$tCRE =="Yes")]
ne.t=CREanno$V4[which(CREanno$promoter_type2 == "enhancer-like" & CREanno$tCRE =="Yes")]
nu.t=CREanno$V4[which(CREanno$promoter_type2 == "unclassed" & CREanno$tCRE =="Yes")]
nc.t=CREanno$V4[which(CREanno$promoter_type2 == "CTCF-alone" & CREanno$tCRE =="Yes")]
np.u=CREanno$V4[which(CREanno$promoter_type2 == "promoter-like" & CREanno$tCRE =="No")]
ne.u=CREanno$V4[which(CREanno$promoter_type2 == "enhancer-like" & CREanno$tCRE =="No")]
nu.u=CREanno$V4[which(CREanno$promoter_type2 == "unclassed" & CREanno$tCRE =="No")]
nc.u=CREanno$V4[which(CREanno$promoter_type2 == "CTCF-alone" & CREanno$tCRE =="No")]

CPGp=CPG[which(CPG$CREID %in% c(np.t, np.u)),]
CPGe=CPG[which(CPG$CREID %in% c(ne.t, ne.u)),]
CPGu=CPG[which(CPG$CREID %in% c(nu.t, nu.u)),]
CPGc=CPG[which(CPG$CREID %in% c(nc.t, nc.u)),]
need=c(20,40,60,80)
need2=unique(CPG$bimodal)

for (i in 1:length(need)){
  for (j in 1: length(need2)){
  CPGc1=CPGc[which(CPGc$cpgNum >= need[i] & CPGc$bimodal == need2[j]),]
  CPGp1=CPGp[which(CPGp$cpgNum >= need[i] & CPGp$bimodal == need2[j]),]
  CPGe1=CPGe[which(CPGe$cpgNum >= need[i] & CPGe$bimodal == need2[j]),]
  CPGu1=CPGu[which(CPGu$cpgNum >= need[i] & CPGu$bimodal == need2[j]),]
  path2=paste0(CpG_folder,"n",need[i],"/",need2[j],"_")
  write.table(CPGc1[order(CPGc1$locS),c(1:4,8)],gzfile(paste0(path2,"summit_5kb_extend.CpGc.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGp1[order(CPGp1$locS),c(1:4,8)],gzfile(paste0(path2,"summit_5kb_extend.CpGp.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGe1[order(CPGe1$locS),c(1:4,8)],gzfile(paste0(path2,"summit_5kb_extend.CpGe.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGu1[order(CPGu1$locS),c(1:4,8)],gzfile(paste0(path2,"summit_5kb_extend.CpGu.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)}}

fake2=cbind("chr1",c(0:10000),c(1:10001),c(-5000:5000))
write.table(fake2,paste0(CpG_folder,"CpG.location10001.fakebed.bed"), row.names=F, col.names=F, sep="\t", quote=F)

#===============================================================================
#bedtools count
setwd(paste0(CpG_folder,"n20"))
system("for file in *summit_5kb_extend.CpG*.bed.gz; do bedtools intersect -wa -a ../CpG.location10001.fakebed.bed -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done")
setwd(paste0(CpG_folder,"n40"))
system("for file in *summit_5kb_extend.CpG*.bed.gz; do bedtools intersect -wa -a ../CpG.location10001.fakebed.bed -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done")
setwd(paste0(CpG_folder,"n60"))
system("for file in *summit_5kb_extend.CpG*.bed.gz; do bedtools intersect -wa -a ../CpG.location10001.fakebed.bed -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done")
setwd(paste0(CpG_folder,"n80"))
system("for file in *summit_5kb_extend.CpG*.bed.gz; do bedtools intersect -wa -a ../CpG.location10001.fakebed.bed -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done")

#===============================================================================
size1=rep(c(length(nc.t),length(ne.t),length(np.t),length(nu.t)),2)
size2=c(length(nc.u),length(ne.u),length(np.u),length(nu.u))
size=rep(c(size1,size2),4)
files=list.files(path=CpG_folder, pattern=".result.bed.gz", recursive=T)
files=files[-grep("/n",files)]
files.index=sapply(strsplit(files, "\\/"),"[",1)
files.names=sapply(strsplit(files, "extend."),"[",2)
files.names=gsub(".result.bed.gz","",files.names)
files.bimodal=sapply(strsplit(files, "\\/"),"[",2)
files.bimodal=sapply(strsplit(files.bimodal, "_summit"),"[",1)
data=read.delim(paste0(CpG_folder, files[1]), header=F, stringsAsFactors = F)
data$group=files.names[1]
data$signalID=files.index[1]
data$bimodal=files.bimodal[1]
data$V5=data$V5/size[1]
for (i in 2:length(files)){
  data1=read.delim(paste0(CpG_folder, files[i]), header=F, stringsAsFactors = F)
  data1$group=files.names[i]
  data1$signalID=files.index[i]
  data1$bimodal=files.bimodal[i]
  data1$V5=data1$V5/size[i]
  data=rbind(data,data1)}

data$anno_region="promoter-like"
data$anno_region[grep("CpGe",data$group)]="enhancer-like"
data$anno_region[grep("CpGu",data$group)]="unclassed"
data$anno_region[grep("CpGc",data$group)]="CTCF-alone"
data$group="CpG_island"
write.table(data, gzfile(paste0(CpG_folder,"CpG_island_result.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
data=read.delim(paste0(CpG_folder,"CpG_island_result.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
