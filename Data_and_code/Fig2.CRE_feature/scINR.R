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
INR_folder=paste0(primary_folder,"Data_and_code/Fig2.CRE_feature/CRE_bed/INR/")

#==========================================
#INR
#prepare the genome wide INR file 
#file size too big, split into 10 files and bedtool intersect directly
pathtt="/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/SCAFE/20231126/Neuron_THP1_S3/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1/bed/initiator/"
for (i in 1:length(files)){
files=list.files(path=pathtt, pattern="part_*")
bed=fread(paste0(pathtt,files[i]), stringsAsFactors = F, header=F)
bed=bed[which(bed$V5>=48),]
write.table(bed,paste0(pathtt,files[i]),col.names=F, row.names=F, sep="\t", quote=F)}

#===============================================================================
#bedtools
setwd(INR_folder)
system("sh bed_intersect.sh")

#===============================================================================
setwd(INR_folder)
ump=read.delim(paste0(CRE_bed,"untranscribed_aCRE.midpoint_summit50.bed.gz"), header=F, stringsAsFactors = F)
tmp=read.delim(paste0(CRE_bed,"transcribed_aCRE.midpoint_summit50.bed.gz"), header=F, stringsAsFactors = F)
tt=read.delim(paste0(CRE_bed,"transcribed_aCRE.t_summit50.stranded.bed.gz"), header=F, stringsAsFactors = F)

ump$bimodal="untranscrib_midpoint"
tmp$bimodal="transcrib_midpoint"
tt$bimodal="transcrib_stranded_summit"
all=rbind(ump,tmp,tt)
all%>%group_by(bimodal)%>%dplyr::summarise(count=n())

files=list.files(path=INR_folder, pattern="bed.gz")
files=files[grep("summit50",files)]
files_bimodal=c(rep("transcrib_midpoint",10),rep("transcrib_stranded_summit",10),rep("untranscrib_midpoint",10))

INR=read.delim(paste0(INR_folder,files[1]), header=F, stringsAsFactors = F)
INR$bimodal=files_bimodal[1]
for (i in 2:length(files)){
  INR1=read.delim(paste0(INR_folder,files[i]), header=F, stringsAsFactors = F)
  INR1$bimodal=files_bimodal[i]
  INR=rbind(INR,INR1)}
INR=left_join(INR, all[,c(2,4,7)], by=c("V4","bimodal"), copy=F, suffix=c("","_ori"))
INR%>%group_by(bimodal)%>%dplyr::summarise(count=length(unique(V4)))

INR$locS=INR$V2-INR$V2_ori-50
INR$locE=INR$V3-INR$V2_ori-50
INR1=INR[which(INR$V6 != "-"),c(1,15,16,4,10,11,12,13)]
INR2=INR[which(INR$V6 == "-"),c(1,16,15,4,10,11,12,13)]

INR2$locE=INR2$locE * (-1)
INR2$locS=INR2$locS * (-1)
INR2$V1="chr1"
INR1$V1="chr1"
colnames(INR2)=colnames(INR1)
INR1$locS=INR1$locS+50
INR1$locE=INR1$locE+50
INR2$locS=INR2$locS+51
INR2$locE=INR2$locE+51

INR=rbind(INR1,INR2)
colnames(INR)[c(4,5,6,7)]=c("CREID","name","motif_score","INR_strand")
write.table(INR[order(INR$locS),],gzfile(paste0(INR_folder,"all_summit_50bb_extend.INR.fakebed.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)
INR$locE=INR$locS+2

CREanno=read.delim(paste0(primary_folder,"scATAC/merged_peak/all_peaks.merged.all.markalone_k16.p.e.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
np.t=CREanno$V4[which(CREanno$promoter_type2 == "promoter-like" & CREanno$tCRE =="Yes")]
ne.t=CREanno$V4[which(CREanno$promoter_type2 == "enhancer-like" & CREanno$tCRE =="Yes")]
nu.t=CREanno$V4[which(CREanno$promoter_type2 == "unclassed" & CREanno$tCRE =="Yes")]
nc.t=CREanno$V4[which(CREanno$promoter_type2 == "CTCF-alone" & CREanno$tCRE =="Yes")]
np.u=CREanno$V4[which(CREanno$promoter_type2 == "promoter-like" & CREanno$tCRE =="No")]
ne.u=CREanno$V4[which(CREanno$promoter_type2 == "enhancer-like" & CREanno$tCRE =="No")]
nu.u=CREanno$V4[which(CREanno$promoter_type2 == "unclassed" & CREanno$tCRE =="No")]
nc.u=CREanno$V4[which(CREanno$promoter_type2 == "CTCF-alone" & CREanno$tCRE =="No")]

INRp=INR[which(INR$CREID %in% c(np.t, np.u)),]
INRe=INR[which(INR$CREID %in% c(ne.t, ne.u)),]
INRu=INR[which(INR$CREID %in% c(nu.t, nu.u)),]
INRc=INR[which(INR$CREID %in% c(nc.t, nc.u)),]

need=c(48.5,49,49.5,50)
need2=unique(INR$bimodal)
for (i in 1:length(need)){
  for (j in 1: length(need2)){
    INRc1=INRc[which(INRc$motif_score >= need[i] & INRc$bimodal == need2[j]),]
    INRp1=INRp[which(INRp$motif_score >= need[i] & INRp$bimodal == need2[j]),]
    INRe1=INRe[which(INRe$motif_score >= need[i] & INRe$bimodal == need2[j]),]
    INRu1=INRu[which(INRu$motif_score >= need[i] & INRu$bimodal == need2[j]),]
    #
    INRc1=INRc1%>%group_by(CREID)%>%dplyr::slice_max(motif_score)
    INRp1=INRp1%>%group_by(CREID)%>%dplyr::slice_max(motif_score)
    INRe1=INRe1%>%group_by(CREID)%>%dplyr::slice_max(motif_score)
    INRu1=INRu1%>%group_by(CREID)%>%dplyr::slice_max(motif_score)
    #if strand overlap for unstranded summit, just take the plus strand
    INRc1=INRc1%>%group_by(CREID)%>%dplyr::slice_max(INR_strand)
    INRp1=INRp1%>%group_by(CREID)%>%dplyr::slice_max(INR_strand)
    INRe1=INRe1%>%group_by(CREID)%>%dplyr::slice_max(INR_strand)
    INRu1=INRu1%>%group_by(CREID)%>%dplyr::slice_max(INR_strand)
    path2=paste0(INR_folder,"n",need[i],"/",need2[j],"_")
    write.table(INRc1[order(INRc1$locS),c(1:8)],gzfile(paste0(path2,"summit_50bp_extend.INRc.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
    write.table(INRp1[order(INRp1$locS),c(1:8)],gzfile(paste0(path2,"summit_50bp_extend.INRp.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
    write.table(INRe1[order(INRe1$locS),c(1:8)],gzfile(paste0(path2,"summit_50bp_extend.INRe.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
    write.table(INRu1[order(INRu1$locS),c(1:8)],gzfile(paste0(path2,"summit_50bp_extend.INRu.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)}}

fake2=cbind("chr1",c(0:100),c(1:101),c(-50:50))
write.table(fake2,paste0(INR_folder,"INR.location101.fakebed.bed"), row.names=F, col.names=F, sep="\t", quote=F)

#===============================================================================
#bedtools count
setwd(paste0(INR_folder,"n48.5"))
system("for file in *fakebed.bed.gz; do bedtools intersect -wa -a ../INR.location101.fakebed.bed -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done")
setwd(paste0(INR_folder,"n49"))
system("for file in *fakebed.bed.gz; do bedtools intersect -wa -a ../INR.location101.fakebed.bed -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done")
setwd(paste0(INR_folder,"n49.5"))
system("for file in *fakebed.bed.gz; do bedtools intersect -wa -a ../INR.location101.fakebed.bed -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done")
setwd(paste0(INR_folder,"n50"))
system("for file in *fakebed.bed.gz; do bedtools intersect -wa -a ../INR.location101.fakebed.bed -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done")

#===============================================================================

size1=rep(c(length(nc.t),length(ne.t),length(np.t),length(nu.t)),2)
size2=c(length(nc.u),length(ne.u),length(np.u),length(nu.u))
size=rep(c(size1,size2),4)
files=list.files(path=INR_folder, pattern=".result.bed.gz", recursive=T)
files.index=sapply(strsplit(files, "\\/"),"[",1)
files.names=sapply(strsplit(files, "extend."),"[",2)
files.names=gsub(".result.bed.gz","",files.names)
files.bimodal=sapply(strsplit(files, "\\/"),"[",2)
files.bimodal=sapply(strsplit(files.bimodal, "_summit"),"[",1)
data=read.delim(paste0(INR_folder, files[1]), header=F, stringsAsFactors = F)
data$group=files.names[1]
data$signalID=files.index[1]
data$bimodal=files.bimodal[1]
data$V5=data$V5/size[1]
for (i in 2:length(files)){
  data1=read.delim(paste0(INR_folder, files[i]), header=F, stringsAsFactors = F)
  data1$group=files.names[i]
  data1$signalID=files.index[i]
  data1$bimodal=files.bimodal[i]
  data1$V5=data1$V5/size[i]
  data=rbind(data,data1)}

data$anno_region="promoter-like"
data$anno_region[grep("INRe",data$group)]="enhancer-like"
data$anno_region[grep("INRu",data$group)]="unclassed"
data$anno_region[grep("INRc",data$group)]="CTCF-alone"
data$group="INR"
write.table(data, gzfile(paste0(INR_folder,"INR_result.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

data=read.delim(paste0(INR_folder,"INR_result.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
flank.plot2 <- ggplot()+
  scale_color_npg()+
  geom_line(data = data, aes(x=V4, y=V5, color=signalID, group=signalID),alpha=0.75, size=0.25)+
  facet_wrap(vars(anno_region,bimodal), ncol=3, nrow=4, scales="free_y")+
  labs(color=NULL, x=NULL, y=NULL, title="TATA box")+
  scale_y_continuous(labels = scales::percent)+
  theme1+theme(legend.position="bottom")
pdf(paste0(path1,"INR_all.aCRE.pdf"), width = 4, height = 4)
print(flank.plot2)
dev.off()

#==================================