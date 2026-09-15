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
TATA_folder=paste0(primary_folder,"Data_and_code/Fig2.CRE_feature/CRE_bed/TATA/")


#==TATA box and initiator=================
#cd /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/SCAFE/20231126/Neuron_THP1_S3/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1/bed/TATAbox
#/home/yip/homer/bin/scanMotifGenomeWide.pl MA0108.3.motif hg38 -bed -keepAll -p 10 >hg38.TATAbox.bed
# the motif file need to be generated, all separators are tab 

#cd /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/SCAFE/20231126/Neuron_THP1_S3/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1/bed/initiator
#/home/yip/homer/bin/scanMotifGenomeWide.pl INR_oldJaspar_10100.motif hg38 -bed -keepAll -p 10 >hg38.INR.bed
# the motif file need to be generated, all separators are tab 
#======================================
setwd(CRE_bed)
tata=fread("/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/SCAFE/20231126/Neuron_THP1_S3/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1/bed/TATAbox/hg38.TATAbox.bed", header=F, stringsAsFactors = F)
tata=tata[which(nchar(tata$V1)<=5),] #main chromosome only
tata$V2=tata$V2-1
tata$V4=paste0("TBP",1:nrow(tata))
write.table(tata, gzfile("hg38.TATAbox.main.bed.gz"), row.names=F, col.names=F, sep="\t", quote=F)

#===============================================================================
setwd(CRE_bed)
ump=read.delim("untranscribed_aCRE.midpoint_summit50.bed.gz", header=F, stringsAsFactors = F)
tmp=read.delim("transcribed_aCRE.midpoint_summit50.bed.gz", header=F, stringsAsFactors = F)
tt=read.delim("transcribed_aCRE.t_summit50.stranded.bed.gz", header=F, stringsAsFactors = F)
ump$bimodal="untranscrib_midpoint"
tmp$bimodal="transcrib_midpoint"
tt$bimodal="transcrib_stranded_summit"
all=rbind(ump,tmp,tt)
all%>%group_by(bimodal)%>%dplyr::summarise(count=n())

#===============================================================================
#bedtools
setwd(CRE_bed)
system("bedtools intersect -wb -a untranscribed_aCRE.midpoint_summit50.bed.gz -b hg38.TATAbox.main.bed.gz | gzip > untranscribed_aCRE.midpoint_summit50.TATA.bed.gz")
system("bedtools intersect -wb -a transcribed_aCRE.t_summit50.stranded.bed.gz -b hg38.TATAbox.main.bed.gz -s | gzip > transcribed_aCRE.t_summit50.stranded.TATA.bed.gz")
system("bedtools intersect -wb -a transcribed_aCRE.midpoint_summit50.bed.gz -b hg38.TATAbox.main.bed.gz | gzip > transcribed_aCRE.midpoint_summit50.TATA.bed.gz")

#==========================
TATAump=read.delim("untranscribed_aCRE.midpoint_summit50.TATA.bed.gz", header=F, stringsAsFactors = F)
TATAtmp=read.delim("transcribed_aCRE.midpoint_summit50.TATA.bed.gz", header=F, stringsAsFactors = F)
TATAtt=read.delim("transcribed_aCRE.t_summit50.stranded.TATA.bed.gz", header=F, stringsAsFactors = F)
length(unique(TATAump$V4)) #320247
length(unique(TATAtmp$V4)) #20538
length(unique(TATAtt$V4)) #16284

TATAump$bimodal="untranscrib_midpoint"
TATAtmp$bimodal="transcrib_midpoint"
TATAtt$bimodal="transcrib_stranded_summit"
TATA=rbind(TATAump,TATAtmp,TATAtt)

TATA=left_join(TATA, all[,c(2,4,7)], by=c("V4","bimodal"), copy=F, suffix=c("","_ori"))
TATA$locS=TATA$V2-TATA$V2_ori-50
TATA$locE=TATA$V3-TATA$V2_ori-50
TATA1=TATA[which(TATA$V6 != "-"),c(1,15,16,4,10,11,12,13)]
TATA2=TATA[which(TATA$V6 == "-"),c(1,16,15,4,10,11,12,13)]

TATA2$locE=TATA2$locE * (-1)
TATA2$locS=TATA2$locS * (-1)
TATA2$V1="chr1"
TATA1$V1="chr1"
colnames(TATA2)=colnames(TATA1)
TATA1$locS=TATA1$locS+50
TATA1$locE=TATA1$locE+50
TATA2$locS=TATA2$locS+51
TATA2$locE=TATA2$locE+51

TATA=rbind(TATA1,TATA2)
colnames(TATA)[c(4,5,6,7)]=c("CREID","name","motif_score","TATA_strand")
write.table(TATA[order(TATA$locS),],paste0(TATA_folder,"all_summit_50bb_extend.TBP.fakebed.tsv.gz"), row.names=F, col.names=T, sep="\t", quote=F)
TATA$locE=TATA$locS+4
#TATA$position=abs(TATA$locS+20)
#taking position will have bias, take the highest score
CREanno=read.delim("/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/JC_merged/all_peaks.merged.all.markalone_k16.p.e.tsv", header=T, stringsAsFactors = F, check.names = F)
np.t=CREanno$V4[which(CREanno$promoter_type2 == "promoter-like" & CREanno$tCRE =="Yes")]
ne.t=CREanno$V4[which(CREanno$promoter_type2 == "enhancer-like" & CREanno$tCRE =="Yes")]
nu.t=CREanno$V4[which(CREanno$promoter_type2 == "unclassed" & CREanno$tCRE =="Yes")]
nc.t=CREanno$V4[which(CREanno$promoter_type2 == "CTCF-alone" & CREanno$tCRE =="Yes")]
np.u=CREanno$V4[which(CREanno$promoter_type2 == "promoter-like" & CREanno$tCRE =="No")]
ne.u=CREanno$V4[which(CREanno$promoter_type2 == "enhancer-like" & CREanno$tCRE =="No")]
nu.u=CREanno$V4[which(CREanno$promoter_type2 == "unclassed" & CREanno$tCRE =="No")]
nc.u=CREanno$V4[which(CREanno$promoter_type2 == "CTCF-alone" & CREanno$tCRE =="No")]

TATAp=TATA[which(TATA$CREID %in% c(np.t, np.u)),]
TATAe=TATA[which(TATA$CREID %in% c(ne.t, ne.u)),]
TATAu=TATA[which(TATA$CREID %in% c(nu.t, nu.u)),]
TATAc=TATA[which(TATA$CREID %in% c(nc.t, nc.u)),]

need=c(1,2,3,4)
need2=unique(TATA$bimodal)
for (i in 1:length(need)){
  for (j in 1: length(need2)){
  TATAc1=TATAc[which(TATAc$motif_score >= need[i] & TATAc$bimodal == need2[j]),]
  TATAp1=TATAp[which(TATAp$motif_score >= need[i] & TATAp$bimodal == need2[j]),]
  TATAe1=TATAe[which(TATAe$motif_score >= need[i] & TATAe$bimodal == need2[j]),]
  TATAu1=TATAu[which(TATAu$motif_score >= need[i] & TATAu$bimodal == need2[j]),]
  #
  TATAc1=TATAc1%>%group_by(CREID)%>%dplyr::slice_max(motif_score)
  TATAp1=TATAp1%>%group_by(CREID)%>%dplyr::slice_max(motif_score)
  TATAe1=TATAe1%>%group_by(CREID)%>%dplyr::slice_max(motif_score)
  TATAu1=TATAu1%>%group_by(CREID)%>%dplyr::slice_max(motif_score)
  #if strand overlap for unstranded summit, just take the plus strand
  TATAc1=TATAc1%>%group_by(CREID)%>%dplyr::slice_max(TATA_strand)
  TATAp1=TATAp1%>%group_by(CREID)%>%dplyr::slice_max(TATA_strand)
  TATAe1=TATAe1%>%group_by(CREID)%>%dplyr::slice_max(TATA_strand)
  TATAu1=TATAu1%>%group_by(CREID)%>%dplyr::slice_max(TATA_strand)
  path2=paste0(TATA_folder,"n",need[i],"/",need2[j],"_")
  write.table(TATAc1[order(TATAc1$locS),c(1:8)],gzfile(paste0(path2,"summit_50bp_extend.TATAc.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(TATAp1[order(TATAp1$locS),c(1:8)],gzfile(paste0(path2,"summit_50bp_extend.TATAp.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(TATAe1[order(TATAe1$locS),c(1:8)],gzfile(paste0(path2,"summit_50bp_extend.TATAe.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(TATAu1[order(TATAu1$locS),c(1:8)],gzfile(paste0(path2,"summit_50bp_extend.TATAu.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)}}
 
fake2=cbind("chr1",c(0:100),c(1:101),c(-50:50))
write.table(fake2,paste0(TATA_folder,"TATA.location101.fakebed.bed"), row.names=F, col.names=F, sep="\t", quote=F)

#===============================================================================
#bedtools count
setwd(paste0(TATA_folder,"n1/"))
system("for file in *fakebed.bed; do bedtools intersect -wa -a ../TATA.location101.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done")
setwd(paste0(TATA_folder,"n2/"))
system("for file in *fakebed.bed; do bedtools intersect -wa -a ../TATA.location101.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done")
setwd(paste0(TATA_folder,"n3/"))
system("for file in *fakebed.bed; do bedtools intersect -wa -a ../TATA.location101.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done")
setwd(paste0(TATA_folder,"n4/"))
system("for file in *fakebed.bed; do bedtools intersect -wa -a ../TATA.location101.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done")

#===============================================================================

size1=rep(c(length(nc.t),length(ne.t),length(np.t),length(nu.t)),2)
size2=c(length(nc.u),length(ne.u),length(np.u),length(nu.u))
size=rep(c(size1,size2),4)
files=list.files(path=TATA_folder, pattern=".result.bed.gz", recursive=T)
files.index=sapply(strsplit(files, "\\/"),"[",1)
files.names=sapply(strsplit(files, "extend."),"[",2)
files.names=gsub(".result.bed.gz","",files.names)
files.bimodal=sapply(strsplit(files, "\\/"),"[",2)
files.bimodal=sapply(strsplit(files.bimodal, "_summit"),"[",1)
data=read.delim(paste0(TATA_folder, files[1]), header=F, stringsAsFactors = F)
data$group=files.names[1]
data$signalID=files.index[1]
data$bimodal=files.bimodal[1]
data$V5=data$V5/size[1]
for (i in 2:length(files)){
  data1=read.delim(paste0(TATA_folder, files[i]), header=F, stringsAsFactors = F)
  data1$group=files.names[i]
  data1$signalID=files.index[i]
  data1$bimodal=files.bimodal[i]
  data1$V5=data1$V5/size[i]
  data=rbind(data,data1)}

data$anno_region="promoter-like"
data$anno_region[grep("TATAe",data$group)]="enhancer-like"
data$anno_region[grep("TATAu",data$group)]="unclassed"
data$anno_region[grep("TATAc",data$group)]="CTCF-alone"
data$group="TATA_box"
write.table(data, gzfile(paste0(TATA_folder,"TATA_result.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)







