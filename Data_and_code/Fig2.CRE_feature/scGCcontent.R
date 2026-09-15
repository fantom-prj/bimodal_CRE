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
GC_folder=paste0(primary_folder,"Data_and_code/Fig2.CRE_feature/CRE_bed/GC/")

#===============================================
#sample the untranscribed_aCRE.midpoint_summit2000.bed.gz into 100000 event for random
setwd(CRE_bed)
bed=read.delim("untranscribed_aCRE.midpoint_summit2000.bed.gz", header=F, stringsAsFactors = F)
bed1=bed[sample(nrow(bed),100000),]
write.table(bed1[order(bed1$V1,bed1$V2),], gzfile("untranscribed_aCRE.midpoint_summit2000sample100000.bed.gz"), row.names=F, col.names=F, sep="\t", quote=F)
#==============================================================================
#GC content tCRE correct version

#cd /analysisdata/fantom6/Interactome/resources/UCSC
#~/bigWigToBedGraph hg38.gc5Base.bw | gzip > hg38.gc5Base.bed.gz


setwd(CRE_bed)
system("bedtools intersect -wb -a untranscribed_aCRE.midpoint_summit2000.bed.gz -b /analysisdata/fantom6/Interactome/resources/UCSC/hg38.gc5Base.bed.gz | gzip > untranscribed_aCRE.midpoint_summit2000.gc.bed.gz")
system("bedtools intersect -wb -a transcribed_aCRE.midpoint_summit2000.bed.gz -b /analysisdata/fantom6/Interactome/resources/UCSC/hg38.gc5Base.bed.gz | gzip > transcribed_aCRE.midpoint_summit2000.gc.bed.gz")
system("bedtools intersect -wb -a transcribed_aCRE.t_summit2000.stranded.bed.gz -b /analysisdata/fantom6/Interactome/resources/UCSC/hg38.gc5Base.bed.gz | gzip > transcribed_aCRE.t_summit2000.stranded.gc.bed.gz")

#random seq -> exclude the tCRE 4001 bp region, exclude all exons from gencode
setwd(CRE_bed)

system("zcat untranscribed_aCRE.midpoint_summit2000.bed.gz /analysisdata/fantom6/Interactome/resources/gencode.v39.annotation.all_exon.bed.gz | sort -k1,1 -k2,2n | gzip > untranscrib_midpoint_gencode_exon.bed.gz")
system("bedtools shuffle -i untranscribed_aCRE.midpoint_summit2000sample100000.bed.gz -g /analysisdata/fantom6/Interactome/resources/hg38.chrom.sizes.bed -noOverlapping -excl untranscrib_midpoint_gencode_exon.bed.gz | gzip > untranscrib_midpoint_random_2kb_extend.bed.gz")
system("bedtools intersect -wa -wb -a untranscrib_midpoint_random_2kb_extend.bed.gz -b /analysisdata/fantom6/Interactome/resources/UCSC/hg38.gc5Base.bed.gz | gzip > untranscrib_midpoint_random_2kb_extend.gc.bed.gz")

system("zcat transcribed_aCRE.midpoint_summit2000.bed.gz /analysisdata/fantom6/Interactome/resources/gencode.v39.annotation.all_exon.bed.gz | sort -k1,1 -k2,2n | gzip > transcrib_midpoint_gencode_exon.bed.gz")
system("bedtools shuffle -i transcribed_aCRE.midpoint_summit2000.bed.gz -g /analysisdata/fantom6/Interactome/resources/hg38.chrom.sizes.bed -noOverlapping -excl transcrib_midpoint_gencode_exon.bed.gz | gzip > transcrib_midpoint_random_2kb_extend.bed.gz")
system("bedtools intersect -wa -wb -a transcrib_midpoint_random_2kb_extend.bed.gz -b /analysisdata/fantom6/Interactome/resources/UCSC/hg38.gc5Base.bed.gz | gzip > transcrib_midpoint_random_2kb_extend.gc.bed.gz")

system("zcat transcribed_aCRE.t_summit2000.stranded.bed.gz /analysisdata/fantom6/Interactome/resources/gencode.v39.annotation.all_exon.bed.gz | sort -k1,1 -k2,2n | gzip > transcrib_stranded_summit_gencode_exon.bed.gz")
system("bedtools shuffle -i transcribed_aCRE.t_summit2000.stranded.bed.gz -g /analysisdata/fantom6/Interactome/resources/hg38.chrom.sizes.bed -noOverlapping -excl transcrib_stranded_summit_gencode_exon.bed.gz | gzip > transcrib_stranded_summit_random_2kb_extend.bed.gz")
system("bedtools intersect -wa -wb -a transcrib_stranded_summit_random_2kb_extend.bed.gz -b /analysisdata/fantom6/Interactome/resources/UCSC/hg38.gc5Base.bed.gz | gzip > transcrib_stranded_summit_random_2kb_extend.gc.bed.gz")

#=============================================================================
setwd(CRE_bed)
ump=read.delim("untranscribed_aCRE.midpoint_summit2000.bed.gz", header=F, stringsAsFactors = F)
tmp=read.delim("transcribed_aCRE.midpoint_summit2000.bed.gz", header=F, stringsAsFactors = F)
tt=read.delim("transcribed_aCRE.t_summit2000.stranded.bed.gz", header=F, stringsAsFactors = F)
ump$bimodal="untranscrib_midpoint"
tmp$bimodal="transcrib_midpoint"
tt$bimodal="transcrib_stranded_summit"
all=rbind(ump,tmp,tt)
all%>%group_by(bimodal)%>%dplyr::summarise(count=n())

GCump=fread("untranscribed_aCRE.midpoint_summit2000.gc.bed.gz", header=F, stringsAsFactors = F)
GCtmp=fread("transcribed_aCRE.midpoint_summit2000.gc.bed.gz", header=F, stringsAsFactors = F)
GCtt=fread("transcribed_aCRE.t_summit2000.stranded.gc.bed.gz", header=F, stringsAsFactors = F)
length(unique(GCump$V4)) #420169
length(unique(GCtmp$V4)) #43697
length(unique(GCtt$V4)) #43697

GCump$bimodal="untranscrib_midpoint"
GCtmp$bimodal="transcrib_midpoint"
GCtt$bimodal="transcrib_stranded_summit"
GC=rbind(GCump,GCtmp,GCtt)
rm(GCump)
rm(GCtmp)
rm(GCtt)

GC=left_join(GC, all[,c(2,4,7)], by=c("V4","bimodal"), copy=F, suffix=c("","_ori"))

GC$locS=GC$V2-GC$V2_ori-2000
GC$locE=GC$V3-GC$V2_ori-2000
GC=GC[,c(1,13,14,4,10,11,6)]
GC2=GC[which(GC$V6 == "-" ),c(1,3,2,4:6)]
GC=GC[which(GC$V6 != "-" ),c(1:6)]
GC2$locE=GC2$locE * (-1)
GC2$locS=GC2$locS * (-1)
GC2$V1="chr1"
GC$V1="chr1"
colnames(GC2)=colnames(GC)
GC$locS=GC$locS+2000
GC$locE=GC$locE+2000
GC2$locS=GC2$locS+2001
GC2$locE=GC2$locE+2001
GC=rbind(GC,GC2)
colnames(GC)[c(4,5)]=c("CREID","percent")

write.table(GC[order(GC$locS),],paste0(GC_folder,"all_summit_2kb_extend.GC.fakebed.tsv"), row.names=F, col.names=T, sep="\t", quote=F)
GC=fread(paste0(GC_folder,"all_summit_2kb_extend.GC.fakebed.tsv"), header=T, stringsAsFactors = F, check.names = F)

CREanno=read.delim(paste0(primary_folder,"scATAC/merged_peak/all_peaks.merged.all.markalone_k16.p.e.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
np.t=CREanno$V4[which(CREanno$promoter_type2 == "promoter-like" & CREanno$tCRE =="Yes")]
ne.t=CREanno$V4[which(CREanno$promoter_type2 == "enhancer-like" & CREanno$tCRE =="Yes")]
nu.t=CREanno$V4[which(CREanno$promoter_type2 == "unclassed" & CREanno$tCRE =="Yes")]
nc.t=CREanno$V4[which(CREanno$promoter_type2 == "CTCF-alone" & CREanno$tCRE =="Yes")]
np.u=CREanno$V4[which(CREanno$promoter_type2 == "promoter-like" & CREanno$tCRE =="No")]
ne.u=CREanno$V4[which(CREanno$promoter_type2 == "enhancer-like" & CREanno$tCRE =="No")]
nu.u=CREanno$V4[which(CREanno$promoter_type2 == "unclassed" & CREanno$tCRE =="No")]
nc.u=CREanno$V4[which(CREanno$promoter_type2 == "CTCF-alone" & CREanno$tCRE =="No")]

GCp=GC[which(GC$CREID %in% c(np.t, np.u)),]
GCe=GC[which(GC$CREID %in% c(ne.t, ne.u)),]
GCu=GC[which(GC$CREID %in% c(nu.t, nu.u)),]
GCc=GC[which(GC$CREID %in% c(nc.t, nc.u)),]

fake2=data.frame(cbind("chr1",c(0:4000),c(1:4001),c(-2000:2000)))
fake2$GCpercent=0
colnames(fake2)=c("chr","start","end","position","GCprecent")
need1=c("transcrib_stranded_summit","transcrib_midpoint","untranscrib_midpoint")
rm(GC)
for (j in 3){
GCpCount=fake2
for (i in 1:4001){GCpCount$GCprecent[i]=mean(GCp$percent[which(GCp$locS<i & GCp$locE>=i & GCp$bimodal == need1[j])],na.rm=T)}
GCeCount=fake2
for (i in 1:4001){GCeCount$GCprecent[i]=mean(GCe$percent[which(GCe$locS<i & GCe$locE>=i & GCe$bimodal == need1[j])],na.rm=T)}
GCuCount=fake2
for (i in 1:4001){GCuCount$GCprecent[i]=mean(GCu$percent[which(GCu$locS<i & GCu$locE>=i & GCu$bimodal == need1[j])],na.rm=T)}
GCcCount=fake2
for (i in 1:4001){GCcCount$GCprecent[i]=mean(GCc$percent[which(GCc$locS<i & GCc$locE>=i & GCc$bimodal == need1[j])],na.rm=T)}
GCpCount$promoter_type="promoter-like"
GCeCount$promoter_type="enhancer-like"
GCuCount$promoter_type="unclassed"
GCcCount$promoter_type="CTCF-alone"
GCCount=rbind(GCpCount,GCeCount,GCuCount,GCcCount)
GCCount$group="anno_region"
write.table(GCCount,gzfile(paste0(GC_folder,need1[j],"_GCcontent_result.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)}


#add random seq GC content
GCump=fread("untranscrib_midpoint_random_2kb_extend.gc.bed.gz", header=F, stringsAsFactors = F)
GCtmp=fread("transcrib_midpoint_random_2kb_extend.gc.bed.gz", header=F, stringsAsFactors = F)
GCtt=fread("transcrib_stranded_summit_random_2kb_extend.gc.bed.gz", header=F, stringsAsFactors = F)
length(unique(GCump$V4)) #84457
length(unique(GCtmp$V4)) #40541
length(unique(GCtt$V4)) #40563

GCump$bimodal="untranscrib_midpoint"
GCtmp$bimodal="transcrib_midpoint"
GCtt$bimodal="transcrib_stranded_summit"
GC=rbind(GCump,GCtmp,GCtt)
rm(GCump)
rm(GCtmp)
rm(GCtt)

GC$locS=GC$V8-GC$V2-2000
GC$locE=GC$V9-GC$V2-2000
GC=GC[,c(1,12,13,4,10,11,6)]
GC2=GC[which(GC$V6 == "-" ),c(1,3,2,4:6)]
GC=GC[which(GC$V6 != "-" ),c(1:6)]
GC2$locE=GC2$locE * (-1)
GC2$locS=GC2$locS * (-1)
GC2$V1="chr1"
GC$V1="chr1"
colnames(GC2)=colnames(GC)
GC$locS=GC$locS+2000
GC$locE=GC$locE+2000
GC2$locS=GC2$locS+2001
GC2$locE=GC2$locE+2001
GC=rbind(GC,GC2)
colnames(GC)[c(4,5)]=c("CREID","percent")

write.table(GC[order(GC$locS),],gzfile(paste0(GC_folder,"all_summit_2kb_extend.random_GC.fakebed.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)
GC=fread(paste0(GC_folder,"all_summit_2kb_extend.random_GC.fakebed.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
GC$locS[which(GC$locS<0)]=0

CREanno=read.delim(paste0(primary_folder,"scATAC/merged_peak/all_peaks.merged.all.markalone_k16.p.e.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
np.t=CREanno$V4[which(CREanno$promoter_type2 == "promoter-like" & CREanno$tCRE =="Yes")]
ne.t=CREanno$V4[which(CREanno$promoter_type2 == "enhancer-like" & CREanno$tCRE =="Yes")]
nu.t=CREanno$V4[which(CREanno$promoter_type2 == "unclassed" & CREanno$tCRE =="Yes")]
nc.t=CREanno$V4[which(CREanno$promoter_type2 == "CTCF-alone" & CREanno$tCRE =="Yes")]
np.u=CREanno$V4[which(CREanno$promoter_type2 == "promoter-like" & CREanno$tCRE =="No")]
ne.u=CREanno$V4[which(CREanno$promoter_type2 == "enhancer-like" & CREanno$tCRE =="No")]
nu.u=CREanno$V4[which(CREanno$promoter_type2 == "unclassed" & CREanno$tCRE =="No")]
nc.u=CREanno$V4[which(CREanno$promoter_type2 == "CTCF-alone" & CREanno$tCRE =="No")]

GCp=GC[which(GC$CREID %in% c(np.t, np.u)),]
GCe=GC[which(GC$CREID %in% c(ne.t, ne.u)),]
GCu=GC[which(GC$CREID %in% c(nu.t, nu.u)),]
GCc=GC[which(GC$CREID %in% c(nc.t, nc.u)),]
rm(GC)

fake2=data.frame(cbind("chr1",c(0:4000),c(1:4001),c(-2000:2000)))
fake2$GCpercent=0
colnames(fake2)=c("chr","start","end","position","GCprecent")

need1=c("transcrib_stranded_summit","transcrib_midpoint","untranscrib_midpoint")

for (j in 1:3){
  GCpCount=fake2
  for (i in 1:4001){GCpCount$GCprecent[i]=mean(GCp$percent[which(GCp$locS<i & GCp$locE>=i & GCp$bimodal == need1[j])],na.rm=T)}
  GCeCount=fake2
  for (i in 1:4001){GCeCount$GCprecent[i]=mean(GCe$percent[which(GCe$locS<i & GCe$locE>=i & GCe$bimodal == need1[j])],na.rm=T)}
  GCuCount=fake2
  for (i in 1:4001){GCuCount$GCprecent[i]=mean(GCu$percent[which(GCu$locS<i & GCu$locE>=i & GCu$bimodal == need1[j])],na.rm=T)}
  GCcCount=fake2
  for (i in 1:4001){GCcCount$GCprecent[i]=mean(GCc$percent[which(GCc$locS<i & GCc$locE>=i & GCc$bimodal == need1[j])],na.rm=T)}
  GCpCount$promoter_type="promoter-like"
  GCeCount$promoter_type="enhancer-like"
  GCuCount$promoter_type="unclassed"
  GCcCount$promoter_type="CTCF-alone"
  GCCount=rbind(GCpCount,GCeCount,GCuCount,GCcCount)
  GCCount$group="random"
  write.table(GCCount,gzfile(paste0(GC_folder,need1[j],"_random_GCcontent_result.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)}


#===============================
files=list.files(path=GC_folder, pattern="GCcontent_result.tsv.gz")
files=files[-1]
files.names=sapply(strsplit(files,"_GC"),"[",1)
files.names=gsub("_random","",files.names)
data1=read.delim(paste0(path1,files[1]), header=T, stringsAsFactors = F, check.names = F)
data1$bimodal=files.names[1]
for (i in 2:length(files)){
  data2=read.delim(paste0(path1,files[i]), header=T, stringsAsFactors = F, check.names = F)
  data2$bimodal=files.names[i]
  data1=rbind(data1,data2)}
data2=data1[which(data1$group=="random" & data1$promoter_type=="enhancer-like"),]
data3=rbind(data2,data2,data2,data2)
data3$promoter_type=c(rep("enhancer-like",8002),rep("promoter-like",8002),rep("CTCF-alone",8002),rep("unclassed",8002))
data4=data3[which(data3$bimodal=="transcrib_midpoint"),]
data4$bimodal="untranscrib_midpoint"
data1=rbind(data1[which(data1$group=="anno_region"),],data3,data4)

write.table(data1,gzfile(paste0(GC_folder,"both_GCcontent_result.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)
data1=read.delim(gzfile(paste0(GC_folder,"both_GCcontent_result.tsv.gz")), header=T, stringsAsFactors = F, check.names = F)
flank.plot2 <- ggplot()+
  scale_color_npg()+
  geom_line(data = data1, aes(x=as.numeric(position), y=GCprecent, color=group, group=group),alpha=0.75, size=0.25)+
  facet_grid( rows=vars(promoter_type), cols=vars(bimodal),scales="free_y")+
  labs(color=NULL, x=NULL, y=NULL, title="GC content")+
  #scale_y_continuous(labels = scales::percent)+
  theme1+theme(legend.position="bottom")
pdf(paste0(path1,"GCcontent_all.aCRE.pdf"), width = 4, height = 4)
print(flank.plot2)
dev.off()

data5=data1[which(data1$promoter_type %in% c("promoter-like","enhancer-like") & data1$bimodal %in% c("transcrib_midpoint","untranscrib_midpoint")),]

flank.plot3 <- ggplot()+
  scale_color_npg()+
  geom_line(data = data5, aes(x=as.numeric(position), y=GCprecent/100, color=group, group=group),alpha=0.75, size=0.25)+
  facet_grid( rows=vars(promoter_type), cols=vars(bimodal))+
  labs(color=NULL, x=NULL, y=NULL, title="GC content")+
  scale_y_continuous(labels = scales::percent)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path1,"GCcontent_selected.aCRE.pdf"), width = 2, height = 1.7)
print(flank.plot3)
dev.off()



#===============================