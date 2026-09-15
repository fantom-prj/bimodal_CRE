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

###color#####
library("ggsci")
mypal = pal_npg("nrc", alpha = 1)(9)
mypal
library("scales")
show_col(mypal)

##########################
#primary_folder=[primary_folder]
primary_folder="/analysisdata/fantom6/Interactome/single_cell_wallace/Figure/"
peak_folder=paste0(primary_folder,"Data_and_code/scATAC/merged_peak/")
CRE_folder=paste0(primary_folder,"Data_and_code/Fig2.CRE_feature/")
path_fig2_data=paste0(primary_folder,"Fig2/data/")
SCAFE_path=paste0(primary_folder,"Data_and_code/SCAFE/sc5end/aggregate/run_full/out/annotate/sc5end.iPSC_NSC_Neuron/")
#===============================================================================
setwd(CRE_folder)

#=====================================
# get the CRE regions defined by SCAFE using TSS
setwd(paste0(SCAFE_path,"bed/"))
CRE_bed=read.delim("sc5end.iPSC_NSC_Neuron.CRE.coord.bed.gz", header=F, stringsAsFactors = F)

# merge the + and - stranded tCREs
setwd(paste0(SCAFE_path,"bed/"))
system("bedtools merge -i sc5end.iPSC_NSC_Neuron.CRE.coord.bed.gz -c 5 -o sum | gzip > sc5end.iPSC_NSC_Neuron.CRE.coord.merged.bed.gz")
system("bedtools intersect -wa -wb -a sc5end.iPSC_NSC_Neuron.CRE.coord.bed.gz -b sc5end.iPSC_NSC_Neuron.CRE.coord.merged.bed.gz | gzip > sc5end.iPSC_NSC_Neuron.CRE.coord.merged.label.bed.gz")

CRE_mbed=read.delim("sc5end.iPSC_NSC_Neuron.CRE.coord.merged.label.bed.gz", header=F, stringsAsFactors = F)
CRE_mbed$CREmID=paste0(CRE_mbed$V13,"_",CRE_mbed$V14,"_",CRE_mbed$V15)
CRE_mbed1=CRE_mbed%>%group_by(V13, V14, V15, CREmID, V16)%>%dplyr::summarise(CREID=paste(V4,collapse=";"))
summary(CRE_mbed1$V15-CRE_mbed1$V14)

write.table(CRE_mbed1[order(CRE_mbed1$V13,CRE_mbed1$V14),],gzfile("sc5end.iPSC_NSC_Neuron.CRE.coord.merged.label.bed.gz"), col.names=F, row.names=F, sep="\t", quote=F)

# tCRE got ATAC support & aCRE got TSS support
setwd(paste0(SCAFE_path,"bed/"))
system(paste0("bedtools intersect -wa -wb -a sc5end.iPSC_NSC_Neuron.CRE.coord.merged.label.bed.gz -b ",peak_folder, "filtered_peaks.merged.all.bed.gz | gzip > sc5end.iPSC_NSC_Neuron.CRE.coord.merged.label.aCRE.bed.gz"))
system(paste0("bedtools intersect -wa -wb -a ",peak_folder, "filtered_peaks.merged.all.bed.gz -b sc5end.iPSC_NSC_Neuron.CRE.coord.merged.label.bed.gz | gzip > filtered_peaks.merged.all.tCRE.bed.gz"))

CRE_mbed1$length=CRE_mbed1$V15-CRE_mbed1$V14
CRE_mbed2=read.delim("sc5end.iPSC_NSC_Neuron.CRE.coord.merged.label.aCRE.bed.gz", header=F, stringsAsFactors = F, check.names = F)
CRE_mbed1$ATAC_support="No"
CRE_mbed1$ATAC_support[which(CRE_mbed1$CREmID %in% CRE_mbed2$V4)]="Yes"
length(which(CRE_mbed1$ATAC_support=="No")) #30632 / 111478 -> 27.5%

aCRE_bed =read.delim(paste0(peak_folder,"filtered_peaks.merged.all.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
aCRE_bed1=read.delim("filtered_peaks.merged.all.tCRE.bed.gz", header=F, stringsAsFactors = F, check.names = F)
aCRE_bed$length=aCRE_bed$V3-aCRE_bed$V2
aCRE_bed$TSS_support="No"
aCRE_bed$TSS_support[which(aCRE_bed$V4 %in% aCRE_bed1$V4)]="Yes"
length(which(aCRE_bed$TSS_support=="No")) #399673 / 463866 -> 86.2%

CRE_mbed1$label="tCRE overlapped with aCRE"
colnames(CRE_mbed1)=colnames(aCRE_bed)
aCRE_bed$label="aCRE overlapped with tCRE"

atCRE=rbind(CRE_mbed1[,-c(5:6)],aCRE_bed[,-c(5:6)])
colnames(atCRE)[6]="support"
write.table(atCRE, gzfile(paste0(path_fig2_data,"separate_identification_aCRE_tCRE.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#=====================================
#prepare summit of CRE for TATA, CGI, GC content and initiator
setwd(peak_folder)
CRE_bed =read.delim("filtered_peaks.merged.all.bed.gz", header=F, stringsAsFactors = F, check.names = F)
CRE_bed$length=CRE_bed$V3-CRE_bed$V2
CRE_bed$mid_point_start=CRE_bed$V2+round(CRE_bed$length/2)
CRE_bed$mid_point_end=CRE_bed$mid_point_start+1

SCAFE_cluster_bed2=read.delim(paste0(primary_folder,"Data_and_code/SCAFE/sc5end/merge/annotate/sc5end.iPSC_NSC_Neuron/bed/sc5end.iPSC_NSC_Neuron.cluster.coord.bed.gz"), header=F, stringsAsFactors = F)
collapse_CTSS=read.delim(paste0(primary_folder,"Data_and_code/SCAFE/sc5end/aggregate/run_full/out/aggregate/sc5end.iPSC_NSC_Neuron/bed/sc5end.iPSC_NSC_Neuron.aggregate.collapse.ctss.bed.gz"), header=F, stringsAsFactors = F, check.names = F)

#=====================================
#bedtools collect real CTSS
setwd(paste0(primary_folder,"Data_and_code/SCAFE/sc5end/aggregate/run_full/out/aggregate/sc5end.iPSC_NSC_Neuron/bed/"))
system(paste0("bedtools intersect -wa -a sc5end.iPSC_NSC_Neuron.aggregate.collapse.ctss.bed.gz -b ",primary_folder,"Data_and_code/SCAFE/sc5end/merge/annotate/sc5end.iPSC_NSC_Neuron/bed/sc5end.iPSC_NSC_Neuron.cluster.coord.bed.gz -s | gzip > sc5end.iPSC_NSC_Neuron.aggregate.collapse.within_cluster.ctss.bed.gz"))

#=====================================
collapse_CTSS2=read.delim(paste0(primary_folder,"Data_and_code/SCAFE/sc5end/aggregate/run_full/out/aggregate/sc5end.iPSC_NSC_Neuron/bed/sc5end.iPSC_NSC_Neuron.aggregate.collapse.within_cluster.ctss.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
sum(collapse_CTSS2$V5) #605,273,491
sum(collapse_CTSS2$V4) #381,562,448

#=====================================
#bedtools intersect to locate greatest signal position in CRE
setwd(peak_folder)
system(paste0("bedtools intersect -wa -wb -a ",primary_folder,"Data_and_code/SCAFE/sc5end/aggregate/run_full/out/aggregate/sc5end.iPSC_NSC_Neuron/bed/sc5end.iPSC_NSC_Neuron.aggregate.collapse.within_cluster.ctss.bed.gz -b filtered_peaks.merged.all.bed.gz | gzip > filtered_peaks.merged.all_withCTSS.bed.gz"))

#===============================================================================
setwd(peak_folder)
CRE_bedt =read.delim("filtered_peaks.merged.all_withCTSS.bed.gz", header=F, stringsAsFactors = F, check.names = F)
sum(CRE_bedt$V5) #586,483,304
sum(CRE_bedt$V4) #368,118,596 -> UMI
sum(CRE_bedt$V5)/sum(collapse_CTSS2$V5) #96.89559%
sum(CRE_bedt$V4)/sum(collapse_CTSS2$V4) #96.47663%
CRE_bedt=CRE_bedt%>%group_by(V10)%>%dplyr::mutate(UMI_count=sum(V4))
CRE_bedt1=CRE_bedt%>%group_by(V10)%>%dplyr::slice_max(V5)
CRE_bedt1=left_join(CRE_bedt1,CRE_bed[,c(4,8)], by=c("V10"="V4"),copy=F)
CRE_bedt1$close_mid=abs(CRE_bedt1$V2-CRE_bedt1$mid_point_start)
CRE_bedt2=CRE_bedt1%>%group_by(V10)%>%dplyr::slice_min(close_mid)
CRE_bedt2=CRE_bedt2%>%group_by(V10)%>%dplyr::slice_min(V2)
CRE_bedt2=CRE_bedt2[which(CRE_bedt2$V10 %in% k3),]
options(scipen=999)
CRE_bed=left_join(CRE_bed,CRE_bedt2[,c(10,2,3,6)], by=c("V4"="V10"),copy=F, suffix=c("","_tCRE_summit"))
CRE_bed1=CRE_bed[which(!is.na(CRE_bed$V3_tCRE_summit)),]
CRE_bed1a=CRE_bed1[,c(1,8,9,4:6)]
write.table(CRE_bed1a[order(CRE_bed1a$V1,CRE_bed1a$mid_point_start),],gzfile("/analysisdata/fantom6/Interactome/single_cell_wallace/CRE_bed/transcribed_aCRE.midpoint_summit.bed.gz"),col.names=F, row.names=F, sep="\t", quote=F)

CRE_bed1a50=CRE_bed1a
CRE_bed1a250=CRE_bed1a
CRE_bed1a2000=CRE_bed1a
CRE_bed1a5000=CRE_bed1a
CRE_bed1a50$mid_point_start=CRE_bed1a50$mid_point_start-50
CRE_bed1a50$mid_point_end=CRE_bed1a50$mid_point_end+50
CRE_bed1a50$mid_point_start[which(CRE_bed1a50$mid_point_start<0)]=0
write.table(CRE_bed1a50[order(CRE_bed1a50$V1,CRE_bed1a50$mid_point_start),],gzfile(paste0(CRE_folder,"CRE_bed/transcribed_aCRE.midpoint_summit50.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)
CRE_bed1a250$mid_point_start=CRE_bed1a250$mid_point_start-250
CRE_bed1a250$mid_point_end=CRE_bed1a250$mid_point_end+250
CRE_bed1a250$mid_point_start[which(CRE_bed1a250$mid_point_start<0)]=0
write.table(CRE_bed1a250[order(CRE_bed1a50$V1,CRE_bed1a50$mid_point_start),],gzfile(paste0(CRE_folder,"CRE_bed/transcribed_aCRE.midpoint_summit250.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)

CRE_bed1a2000$mid_point_start=CRE_bed1a2000$mid_point_start-2000
CRE_bed1a2000$mid_point_end=CRE_bed1a2000$mid_point_end+2000
CRE_bed1a2000$mid_point_start[which(CRE_bed1a2000$mid_point_start<0)]=0
write.table(CRE_bed1a2000[order(CRE_bed1a2000$V1,CRE_bed1a2000$mid_point_start),],gzfile(paste0(CRE_folder,"CRE_bed/transcribed_aCRE.midpoint_summit2000.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)
CRE_bed1a5000$mid_point_start=CRE_bed1a5000$mid_point_start-5000
CRE_bed1a5000$mid_point_end=CRE_bed1a5000$mid_point_end+5000
CRE_bed1a5000$mid_point_start[which(CRE_bed1a5000$mid_point_start<0)]=0
write.table(CRE_bed1a5000[order(CRE_bed1a5000$V1,CRE_bed1a5000$mid_point_start),],gzfile(paste0(CRE_folder,"CRE_bed/transcribed_aCRE.midpoint_summit5000.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)

CRE_bed1b=CRE_bed1[,c(1,10,11,4,5,12)]
write.table(CRE_bed1b[order(CRE_bed1b$V1,CRE_bed1b$V2_tCRE_summit),],gzfile(paste0(CRE_folder,"CRE_bed/transcribed_aCRE.t_summit.stranded.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)
CRE_bed1b50=CRE_bed1b
CRE_bed1b2000=CRE_bed1b
CRE_bed1b5000=CRE_bed1b
CRE_bed1b50$V2_tCRE_summit=CRE_bed1b50$V2_tCRE_summit-50
CRE_bed1b50$V3_tCRE_summit=CRE_bed1b50$V3_tCRE_summit+50
CRE_bed1b50$V2_tCRE_summit[which(CRE_bed1b50$V2_tCRE_summit<0)]=0
write.table(CRE_bed1b50[order(CRE_bed1b50$V1,CRE_bed1b50$V2_tCRE_summit),],gzfile(paste0(CRE_folder,"CRE_bed/transcribed_aCRE.t_summit50.stranded.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)
CRE_bed1b2000$V2_tCRE_summit=CRE_bed1b2000$V2_tCRE_summit-2000
CRE_bed1b2000$V3_tCRE_summit=CRE_bed1b2000$V3_tCRE_summit+2000
CRE_bed1b2000$V2_tCRE_summit[which(CRE_bed1b2000$V2_tCRE_summit<0)]=0
write.table(CRE_bed1b2000[order(CRE_bed1b2000$V1,CRE_bed1b2000$V2_tCRE_summit),],gzfile(paste0(CRE_folder,"CRE_bed/transcribed_aCRE.t_summit2000.stranded.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)
CRE_bed1b5000$V2_tCRE_summit=CRE_bed1b5000$V2_tCRE_summit-5000
CRE_bed1b5000$V3_tCRE_summit=CRE_bed1b5000$V3_tCRE_summit+5000
CRE_bed1b5000$V2_tCRE_summit[which(CRE_bed1b5000$V2_tCRE_summit<0)]=0
write.table(CRE_bed1b5000[order(CRE_bed1b5000$V1,CRE_bed1b5000$V2_tCRE_summit),],gzfile(paste0(CRE_folder,"CRE_bed/transcribed_aCRE.t_summit5000.stranded.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)

CRE_bed2=CRE_bed[which(is.na(CRE_bed$V3_tCRE_summit)),c(1,8,9,4:6)]
write.table(CRE_bed2[order(CRE_bed2$V1,CRE_bed2$mid_point_start),],gzfile(paste0(CRE_folder,"CRE_bed/untranscribed_aCRE.midpoint_summit.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)
CRE_bed50=CRE_bed2
CRE_bed250=CRE_bed2
CRE_bed2000=CRE_bed2
CRE_bed5000=CRE_bed2
CRE_bed50$mid_point_start=CRE_bed50$mid_point_start-50
CRE_bed50$mid_point_end=CRE_bed50$mid_point_end+50
CRE_bed50$mid_point_start[which(CRE_bed50$mid_point_start<0)]=0
write.table(CRE_bed50[order(CRE_bed50$V1,CRE_bed50$mid_point_start),],gzfile(paste0(CRE_folder,"CRE_bed/untranscribed_aCRE.midpoint_summit50.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)
CRE_bed250$mid_point_start=CRE_bed250$mid_point_start-250
CRE_bed250$mid_point_end=CRE_bed250$mid_point_end+250
CRE_bed250$mid_point_start[which(CRE_bed250$mid_point_start<0)]=0
write.table(CRE_bed250[order(CRE_bed250$V1,CRE_bed250$mid_point_start),],gzfile(paste0(CRE_folder,"CRE_bed/untranscribed_aCRE.midpoint_summit250.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)

CRE_bed2000$mid_point_start=CRE_bed2000$mid_point_start-2000
CRE_bed2000$mid_point_end=CRE_bed2000$mid_point_end+2000
CRE_bed2000$mid_point_start[which(CRE_bed2000$mid_point_start<0)]=0
write.table(CRE_bed2000[order(CRE_bed2000$V1,CRE_bed2000$mid_point_start),],gzfile(paste0(CRE_folder,"CRE_bed/untranscribed_aCRE.midpoint_summit2000.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)
CRE_bed5000$mid_point_start=CRE_bed5000$mid_point_start-5000
CRE_bed5000$mid_point_end=CRE_bed5000$mid_point_end+5000
CRE_bed5000$mid_point_start[which(CRE_bed5000$mid_point_start<0)]=0
write.table(CRE_bed5000[order(CRE_bed5000$V1,CRE_bed5000$mid_point_start),],gzfile(paste0(CRE_folder,"CRE_bed/untranscribed_aCRE.midpoint_summit5000.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)

#===============================================================================
#After running ./scINR.R, ./scGCcontent.R, ./scTATA.R, ./scCpG.R

CGI=read.delim(paste0(CRE_folder,"CRE_bed/CpG_result/CpG_island_result.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
TATA=read.delim(paste0(CRE_folder,"CRE_bed/TATA/TATA_result.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
INR=read.delim(paste0(CRE_folder,"CRE_bed/INR/INR_result.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
GC=read.delim(paste0(CRE_folder,"CRE_bed/GC/both_GCcontent_result.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
write.table(GC, gzfile(paste0(path_fig2_data,"both_GCcontent_result.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
# -> for fig 2e

#===============================================================================
# add feature to all CREs - ChromHMM
# run ChromHMM according to ../Fig2.ChromHMM/chromHMM.R

#==
# bash
# link ChromHMM results to CRE
setwd(peak_folder)
system(paste0("bedtools intersect -wa -wb -a filtered_peaks.merged.all.bed.gz -b ",chromHMM_folder,"zenbu/markalone_iPS_k16.bed.gz | gzip > filtered_peaks.merged.all.markalone_k16ips.bed.gz"))
system(paste0("bedtools intersect -wa -wb -a filtered_peaks.merged.all.bed.gz -b ",chromHMM_folder,"zenbu/markalone_NSC_k16.bed.gz | gzip > filtered_peaks.merged.all.markalone_k16nsc.bed.gz"))
system(paste0("bedtools intersect -wa -wb -a filtered_peaks.merged.all.bed.gz -b ",chromHMM_folder,"zenbu/markalone_Neuron_k16.bed.gz | gzip > filtered_peaks.merged.all.markalone_k16nrn.bed.gz"))

#==
setwd(peak_folder)
aCRE=read.delim("peaks.merged.all.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
aCRE.i=read.delim("filtered_peaks.merged.all.markalone_k16ips.bed.gz", header=F, stringsAsFactors = F, check.names = F)
aCRE.s=read.delim("filtered_peaks.merged.all.markalone_k16nsc.bed.gz", header=F, stringsAsFactors = F, check.names = F)
aCRE.n=read.delim("filtered_peaks.merged.all.markalone_k16nrn.bed.gz", header=F, stringsAsFactors = F, check.names = F)

aCRE.i1=aCRE.i%>%group_by(V4)%>%dplyr::summarise(type=paste(unique(V10), collapse=";"))
aCRE.s1=aCRE.s%>%group_by(V4)%>%dplyr::summarise(type=paste(unique(V10), collapse=";"))
aCRE.n1=aCRE.n%>%group_by(V4)%>%dplyr::summarise(type=paste(unique(V10), collapse=";"))

aCRE.i1$type2="CTCF"
aCRE.i1$type2[grep("enhancer", aCRE.i1$type)]="enhancer"
aCRE.i1$type2[grep("romoter", aCRE.i1$type)]="promoter"
aCRE.i1%>%group_by(type2)%>%dplyr::summarise(count=n())
aCRE.s1$type2="CTCF"
aCRE.s1$type2[grep("enhancer", aCRE.s1$type)]="enhancer"
aCRE.s1$type2[grep("romoter", aCRE.s1$type)]="promoter"
aCRE.s1%>%group_by(type2)%>%dplyr::summarise(count=n())
aCRE.n1$type2="CTCF"
aCRE.n1$type2[grep("enhancer", aCRE.n1$type)]="enhancer"
aCRE.n1$type2[grep("romoter", aCRE.n1$type)]="promoter"
aCRE.n1%>%group_by(type2)%>%dplyr::summarise(count=n())

aCRE=left_join(aCRE,aCRE.i1[,c(1,3)],by="V4",copy=F)
aCRE=left_join(aCRE,aCRE.s1[,c(1,3)],by="V4",copy=F)
aCRE=left_join(aCRE,aCRE.n1[,c(1,3)],by="V4",copy=F)
colnames(aCRE)[c(9:11)]=c("k16ips","k16nsc","k16nrn")

aCRE$promoter_type2="unclassed"
aCRE$promoter_type2[which(aCRE$k16ips == "CTCF" | aCRE$k16nsc == "CTCF" | aCRE$k16nrn == "CTCF")]="CTCF-alone"
aCRE$promoter_type2[which(aCRE$k16ips == "enhancer" | aCRE$k16nsc == "enhancer" | aCRE$k16nrn == "enhancer")]="enhancer-like"
aCRE$promoter_type2[which(aCRE$k16ips == "promoter" | aCRE$k16nsc == "promoter" | aCRE$k16nrn == "promoter")]="promoter-like"

aCRE%>%group_by(promoter_type)%>%dplyr::summarise(count=n())
aCRE%>%group_by(promoter_type2)%>%dplyr::summarise(count=n())
aCRE[which(aCRE$tCRE=="Yes"),]%>%group_by(promoter_type)%>%dplyr::summarise(count=n())
aCRE[which(aCRE$tCRE=="Yes"),]%>%group_by(promoter_type2)%>%dplyr::summarise(count=n())

#write.table(aCRE[which(aCRE$promoter_type2 != "unclassed"),c(1,2,3,12,5,6)],gzfile("filtered_peaks.merged.all.markalone_k16.p.e.bed.gz"), col.names=F, row.names=F, sep="\t", quote=F)
write.table(aCRE,gzfile("peaks.merged.all.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)
write.table(aCRE,gzfile(paste0(path_fig2_data,"peaks.merged.all.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#==============
#CRE chromatin state per cell type
aCREe=left_join(aCRE[,c(4,8)] ,aCRE.i1,by="V4",copy=F)
aCREe=left_join(aCREe,aCRE.s1,by="V4",copy=F)
aCREe=left_join(aCREe,aCRE.n1,by="V4",copy=F)
colnames(aCREe)[c(3:8)]=c("k16ips_all","k16ips","k16nsc_all","k16nsc","k16nrn_all","k16nrn")
aCREe$promoter_type2="unclassed"
aCREe$promoter_type2[which(aCREe$k16ips == "CTCF" | aCREe$k16nsc == "CTCF" | aCREe$k16nrn == "CTCF")]="CTCF-alone"
aCREe$promoter_type2[which(aCREe$k16ips == "enhancer" | aCREe$k16nsc == "enhancer" | aCREe$k16nrn == "enhancer")]="enhancer-like"
aCREe$promoter_type2[which(aCREe$k16ips == "promoter" | aCREe$k16nsc == "promoter" | aCREe$k16nrn == "promoter")]="promoter-like"

# select chromatin state of CRE when CRE overlap more than one state
aCREe$k16ips_all[grep("Promoter",aCREe$k16ips_all)]="Promoter"
aCREe$k16ips_all[grep("Flanking_promoter",aCREe$k16ips_all)]="Promoter"
aCREe$k16ips_all[grep("Bivalent_promoter",aCREe$k16ips_all)]="Bivalent_promoter"
aCREe$k16ips_all[grep("Active_enhancer",aCREe$k16ips_all)]="Active_enhancer"
aCREe$k16ips_all[grep("Repressed_enhancer",aCREe$k16ips_all)]="Repressed_enhancer"
aCREe$k16ips_all[grep("Primed_enhancer",aCREe$k16ips_all)]="Primed_enhancer"
aCREe$k16ips_all[grep("CTCF",aCREe$k16ips_all)]="CTCF"

aCREe$k16nsc_all[grep("Promoter",aCREe$k16nsc_all)]="Promoter"
aCREe$k16nsc_all[grep("Flanking_promoter",aCREe$k16nsc_all)]="Promoter"
aCREe$k16nsc_all[grep("Bivalent_promoter",aCREe$k16nsc_all)]="Bivalent_promoter"
aCREe$k16nsc_all[grep("Active_enhancer",aCREe$k16nsc_all)]="Active_enhancer"
aCREe$k16nsc_all[grep("Repressed_enhancer",aCREe$k16nsc_all)]="Repressed_enhancer"
aCREe$k16nsc_all[grep("Primed_enhancer",aCREe$k16nsc_all)]="Primed_enhancer"
aCREe$k16nsc_all[grep("CTCF",aCREe$k16nsc_all)]="CTCF"

aCREe$k16nrn_all[grep("Promoter",aCREe$k16nrn_all)]="Promoter"
aCREe$k16nrn_all[grep("Flanking_promoter",aCREe$k16nrn_all)]="Promoter"
aCREe$k16nrn_all[grep("Bivalent_promoter",aCREe$k16nrn_all)]="Bivalent_promoter"
aCREe$k16nrn_all[grep("Active_enhancer",aCREe$k16nrn_all)]="Active_enhancer"
aCREe$k16nrn_all[grep("Repressed_enhancer",aCREe$k16nrn_all)]="Repressed_enhancer"
aCREe$k16nrn_all[grep("Primed_enhancer",aCREe$k16nrn_all)]="Primed_enhancer"
aCREe$k16nrn_all[grep("CTCF",aCREe$k16nrn_all)]="CTCF"

write.table(aCREe, gzfile(paste0(peak_folder,"peaks.merged.all.markalone_k16.subclass.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(aCREe, gzfile(paste0(path_fig2_data,"peaks.merged.all.markalone_k16.subclass.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===================
# add CRE transcription per cell type
setwd(peak_folder)
aCRE_CT=aCREe[,c(1,3,5,7,9)]
colnames(aCRE_CT)[c(1:4)]=c("peakID","iPSC_CT","NSC_CT","Neuron_CT")

aCRE_cell=read.csv("scATAC_merged_peaks_in_each_sample.csv", header=T, stringsAsFactors = F)
patht=paste0(primary_folder,"Data_and_code/SCAFE/sc5end/count/aCRE_original/")
tacre=list.files(path=patht, pattern="genes.log.tsv", recursive=T)
tacre_name=sapply(strsplit(tacre,"\\/"),"[",1)
tacre1=read.delim(paste0(patht,tacre[1]), header=T, stringsAsFactors = F, check.names = F)
tacre1=tacre1[,c(1,3)]
colnames(tacre1)[2]=tacre_name[1]
for (i in 2:6){
  tacre2=read.delim(paste0(patht,tacre[i]), header=T, stringsAsFactors = F, check.names = F)
  tacre2=tacre2[,c(1,3)]
  colnames(tacre2)[2]=tacre_name[i]
  tacre1=left_join(tacre1,tacre2,by="countRegion_ID", copy=F)}
tacre1$iPS=tacre1$iPS1+tacre1$iPS2
tacre1$NSC=tacre1$NSC1+tacre1$NSC2
tacre1$Neuron=tacre1$Neuron1+tacre1$Neuron2
tacre1=tacre1[,c(1,8:10)]
tacre1$iPS[which(tacre1$iPS >= 1)]=1
tacre1$NSC[which(tacre1$NSC >= 1)]=1
tacre1$Neuron[which(tacre1$Neuron >= 1)]=1
colnames(tacre1)=c("peakID","iPSC_t","NSC_t","Neuron_t")
colnames(aCRE_cell)=c("peakID","iPSC_a","NSC_a","Neuron_a")

aCRE_cell=left_join(aCRE_cell,tacre1, by="peakID", copy=F)
aCRE_cell=left_join(aCRE_cell,aCRE_CT, by="peakID", copy=F)
aCRE_cell$chr=sapply(strsplit(aCRE_cell$peakID,"-"),"[",1)
aCRE_cell$analysis="included"
aCRE_cell$analysis[which(nchar(aCRE_cell$chr)>=6)]="excluded"
aCRE_cell=left_join(aCRE_cell,aCRE[,c(4,7,8)],by=c("peakID"="V4"),copy=F)
colnames(aCRE_cell)[c(12:13)]=c("promoter_type_CT","promoter_type_SCREEN")
write.table(aCRE_cell,gzfile("all_aCRE.TSS.access.markalone_cuttag.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)
# -> datafile for cell type details

setwd(paste0(primary_folder,"Data_and_code/scATAC/merged_peak"))
aCRE=read.delim("peaks.merged.all.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
scACREv2=aCRE_cell
scACREv2%>%group_by(promoter_type_CT,promoter_type_SCREEN)%>%dplyr::summarise(count=n())
scACREv3=scACREv2[,c(1,8:14)]
write.table(scACREv3,gzfile("all_aCRE.PE.final.tsv.gz"), col.names=T, row.names = F, sep="\t", quote = F)
# -> data file for overall CRE features

#===============
#bash
#add SE
setwd(peak_folder)
system(paste0("bedtools intersect -c -a filtered_peaks.merged.all.bed.gz -b ",primary_folder,"Data_and_code/SE_identification/iPS_k27ac5/iPS_27ac_SE.bed | gzip > filtered_peaks.merged.all.SE_ips27ac.bed.gz"))
system(paste0("bedtools intersect -c -a filtered_peaks.merged.all.bed.gz -b ",primary_folder,"Data_and_code/SE_identification/NSC_k27ac5/NSC_27ac_SE.bed | gzip > filtered_peaks.merged.all.SE_nsc27ac.bed.gz"))
system(paste0("bedtools intersect -c -a filtered_peaks.merged.all.bed.gz -b ",primary_folder,"Data_and_code/SE_identification/NRN_k27ac5/NRN_27ac_SE.bed | gzip > filtered_peaks.merged.all.SE_neur27ac.bed.gz"))

#==============================
setwd(peak_folder)
iPS.SE=read.delim("filtered_peaks.merged.all.SE_ips27ac.bed.gz", header=F, stringsAsFactors = F)
NSC.SE=read.delim("filtered_peaks.merged.all.SE_nsc27ac.bed.gz", header=F, stringsAsFactors = F)
NRN.SE=read.delim("filtered_peaks.merged.all.SE_neur27ac.bed.gz", header=F, stringsAsFactors = F)
scACREv3$iPSC_SE="No"
scACREv3$iPSC_SE[which(scACREv3$peakID %in% iPS.SE$V4[which(iPS.SE$V7>0)])]="Yes"
scACREv3$NSC_SE="No"
scACREv3$NSC_SE[which(scACREv3$peakID %in% NSC.SE$V4[which(NSC.SE$V7>0)])]="Yes"
scACREv3$Neuron_SE="No"
scACREv3$Neuron_SE[which(scACREv3$peakID %in% NRN.SE$V4[which(NRN.SE$V7>0)])]="Yes"
colnames(scACREv3)=gsub("NRN","Neuron",colnames(scACREv3))

#===============================================================================
#add CGI
setwd(peak_folder)
system(paste0("bedtools intersect -wo -a filtered_peaks.merged.all.bed.gz -b ",CRE_folder,"CRE_bed/cpgIslandExt.hg38.main_chr.bed.gz | gzip > filtered_peaks.merged.all.cgi.bed.gz"))

#===============================================================================
setwd(peak_folder)
CGI=read.delim("filtered_peaks.merged.all.cgi.bed.gz", header=F, stringsAsFactors = F)
scACREv3$CGI="No"
scACREv3$CGI[which(scACREv3$peakID %in% unique(CGI$V4[which(CGI$V17>200)]))]="Yes"
scACREv3%>%group_by(promoter_type_CT,CGI)%>%dplyr::summarise(count=n())

#===============================================================================
#add TATA box
setwd(peak_folder)
tata_info=read.delim(paste0(CRE_folder,"CRE_bed/TATA/all_summit_50bb_extend.TBP.fakebed.tsv.gz"), header=T, stringsAsFactors = F)
tata_info1=tata_info[which(tata_info$motif_score >=3),]
tata_info1=tata_info1[which(tata_info1$locS>15 & tata_info1$locS < 24),] #only take the 15bp region
tata_info1=tata_info1[which(tata_info1$bimodal =="transcrib_stranded_summit"),]
tata_info2=tata_info1%>%group_by(CREID)%>%dplyr::slice_max(motif_score)

scACREv3$TATA_box="No"
scACREv3$TATA_box[which(scACREv3$peakID %in% tata_info2$CREID)]="Yes"
scACREv3%>%group_by(promoter_type_CT,TATA_box)%>%dplyr::summarise(count=n())

write.table(scACREv3,gzfile("all_aCRE.PE.final.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)
write.table(scACREv3[,c(1:9,18,19,13,14,23:27)],gzfile("all_aCRE.PE.final.select.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
#add TATA box prediction for un-transcribed CRE
scACREv3=read.delim("all_aCRE.PE.final.tsv.gz", header=T, stringsAsFactors = F)
TATA_predict=read.delim(paste0(CRE_folder,"CRE_bed/TATA_INR_prediction_untranscribed_CRE_midpoint250_about_cutoff95.tsv.gz"), header=T, stringsAsFactors = F)

scACREv3$TATA_box[which(scACREv3$peakID %in% unique(TATA_predict$CREID))]="Yes"
write.table(scACREv3,gzfile("all_aCRE.PE.final.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
#CRE intersect with GENCODE v39
setwd(peak_folder)
system(paste0("bedtools intersect -wa -wb -a filtered_peaks.merged.all.bed.gz -b ",primary_folder,"Data_and_code/resources/gencode.v39.annotation.bed6.exon.bed.gz | gzip > filtered_peaks.merged_v39exon.bed.gz"))
system(paste0("bedtools intersect -wa -wb -a filtered_peaks.merged.all.bed.gz -b ",primary_folder,"Data_and_code/resources/gencode.v39.annotation.bed6.intron.bed.gz | gzip > filtered_peaks.merged_v39intron.bed.gz"))
system(paste0("bedtools intersect -wa -wb -a filtered_peaks.merged.all.bed.gz -b ",primary_folder,"Data_and_code/resources/gencode.v39.annotation.bed6.3UTR.bed.gz | gzip > filtered_peaks.merged_v393UTR.bed.gz"))
system(paste0("bedtools intersect -wa -wb -a filtered_peaks.merged.all.bed.gz -b ",primary_folder,"Data_and_code/resources/gencode.v39.annotation.bed6.5UTR.bed.gz | gzip > filtered_peaks.merged_v395UTR.bed.gz"))
system(paste0("bedtools intersect -wa -wb -a filtered_peaks.merged.all.bed.gz -b ",primary_folder,"Data_and_code/resources/gencode.v39.annotation.5n.bed.gz | gzip > filtered_peaks.merged_v395n.bed.gz"))

#==========
setwd(peak_folder)
exon=read.delim("filtered_peaks.merged_v39exon.bed.gz", header=F, stringsAsFactors = F)
intron=read.delim("filtered_peaks.merged_v39intron.bed.gz", header=F, stringsAsFactors = F)
UTR3=read.delim("filtered_peaks.merged_v393UTR.bed.gz", header=F, stringsAsFactors = F)
UTR5=read.delim("filtered_peaks.merged_v395UTR.bed.gz", header=F, stringsAsFactors = F)
end5=read.delim("filtered_peaks.merged_v395n.bed.gz", header=F, stringsAsFactors = F)

scACREv3=read.delim("all_aCRE.PE.final.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
scACREv3$gencode_end5="No"
scACREv3$gencode_end5[which(scACREv3$peakID %in% unique(end5$V4))]="Yes"
scACREv3$gencode_5UTR="No"
scACREv3$gencode_5UTR[which(scACREv3$peakID %in% unique(UTR5$V4))]="Yes"
scACREv3$gencode_3UTR="No"
scACREv3$gencode_3UTR[which(scACREv3$peakID %in% unique(UTR3$V4))]="Yes"
scACREv3$gencode_exon="No"
scACREv3$gencode_exon[which(scACREv3$peakID %in% unique(exon$V4))]="Yes"
scACREv3$gencode_intron="No"
scACREv3$gencode_intron[which(scACREv3$peakID %in% unique(intron$V4))]="Yes"

scACREv3$genomic_region="Intergenic"
scACREv3$genomic_region[which(scACREv3$gencode_intron == "Yes")]="Intronic"
scACREv3$genomic_region[which(scACREv3$gencode_exon == "Yes")]="Exonic"
scACREv3$genomic_region[which(scACREv3$gencode_3UTR == "Yes")]="3'UTR"
scACREv3$genomic_region[which(scACREv3$gencode_5UTR == "Yes")]="5'UTR"
scACREv3$genomic_region[which(scACREv3$gencode_end5 == "Yes")]="5'end"
scACREv3%>%group_by(promoter_type_CT,genomic_region)%>%dplyr::summarise(count=n())
write.table(scACREv3,gzfile("all_aCRE.PE.final.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)
write.table(scACREv3,gzfile(paste0(path_fig2_data,"all_aCRE.PE.final.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

# -> datafile for overall CRE features

#===============================================================================
#FE enrichment between transcription and the ChromHMM result, per cell type
setwd(peak_folder)
aCRE_cell=read.delim("all_aCRE.TSS.access.markalone_cuttag.tsv.gz", header=T, stringsAsFactors = F)
aCRE_cell=left_join(aCRE_cell, scACREv3[,c(1,30:40)], by="peakID", copy=F)

aCRE_cell$iPSC.Active="No"
aCRE_cell$iPSC.Active[which(aCRE_cell$iPSC_CT=="Active_enhancer" & aCRE_cell$promoter_type_CT != "promoter-like")]="Yes"
aCRE_cell$iPSC.Active[grep("Promoter",aCRE_cell$iPSC_CT)]="Yes"
aCRE_cell$iPSC.Primed="No"
aCRE_cell$iPSC.Primed[grep("Primed",aCRE_cell$iPSC_CT)]="Yes"
aCRE_cell$iPSC.Repressed="No"
aCRE_cell$iPSC.Repressed[grep("Repressed",aCRE_cell$iPSC_CT)]="Yes"
aCRE_cell$iPSC.Bivalent="No"
aCRE_cell$iPSC.Bivalent[grep("Bivalent",aCRE_cell$iPSC_CT)]="Yes"
aCRE_cell$iPSC.Null="No"
aCRE_cell$iPSC.Null[which(is.na(aCRE_cell$iPSC_CT))]="Yes"
aCRE_cell$iPSC.Null[which(aCRE_cell$iPSC_CT == "CTCF")]="Yes"

aCRE_cell$NSC.Active="No"
aCRE_cell$NSC.Active[which(aCRE_cell$NSC_CT=="Active_enhancer" & aCRE_cell$promoter_type_CT != "promoter-like")]="Yes"
aCRE_cell$NSC.Active[grep("Promoter",aCRE_cell$NSC_CT)]="Yes"
aCRE_cell$NSC.Primed="No"
aCRE_cell$NSC.Primed[grep("Primed",aCRE_cell$NSC_CT)]="Yes"
aCRE_cell$NSC.Repressed="No"
aCRE_cell$NSC.Repressed[grep("Repressed",aCRE_cell$NSC_CT)]="Yes"
aCRE_cell$NSC.Bivalent="No"
aCRE_cell$NSC.Bivalent[grep("Bivalent",aCRE_cell$NSC_CT)]="Yes"
aCRE_cell$NSC.Null="No"
aCRE_cell$NSC.Null[which(is.na(aCRE_cell$NSC_CT))]="Yes"
aCRE_cell$NSC.Null[which(aCRE_cell$NSC_CT == "CTCF")]="Yes"

aCRE_cell$Neuron.Active="No"
aCRE_cell$Neuron.Active[which(aCRE_cell$Neuron_CT=="Active_enhancer" & aCRE_cell$promoter_type_CT != "promoter-like")]="Yes"
aCRE_cell$Neuron.Active[grep("Promoter",aCRE_cell$Neuron_CT)]="Yes"
aCRE_cell$Neuron.Primed="No"
aCRE_cell$Neuron.Primed[grep("Primed",aCRE_cell$Neuron_CT)]="Yes"
aCRE_cell$Neuron.Repressed="No"
aCRE_cell$Neuron.Repressed[grep("Repressed",aCRE_cell$Neuron_CT)]="Yes"
aCRE_cell$Neuron.Bivalent="No"
aCRE_cell$Neuron.Bivalent[grep("Bivalent",aCRE_cell$Neuron_CT)]="Yes"
aCRE_cell$Neuron.Null="No"
aCRE_cell$Neuron.Null[which(is.na(aCRE_cell$Neuron_CT))]="Yes"
aCRE_cell$Neuron.Null[which(aCRE_cell$Neuron_CT == "CTCF")]="Yes"

write.table(aCRE_cell, "all_aCRE.TSS.access.markalone_cuttag.tsv.gz", col.names=T, row.names=F, sep="\t", quote=F)

#====
# prepare fisher exact between transcription and chromatin state
# divided into cell types, only take cell type specific accessible CREs

cell=c("iPSC","NSC","Neuron")
pro_type=c("enhancer-like","promoter-like")
pro_type2=c("enhancer","romoter")
summaryCRE=data.frame(matrix(nrow=0,ncol=7))
colnames(summaryCRE)=c("variable","No_No","No_Yes","Yes_No","Yes_Yes","promoter_type_CT","cell")

for (k in 1:length(pro_type)){
  for (i in 1: length(cell)){
    aCRE_cell1=aCRE_cell[which(aCRE_cell$analysis == "included" & aCRE_cell$promoter_type_CT == pro_type[k]),]
    aCRE_cell1=aCRE_cell1[,grep(cell[i],colnames(aCRE_cell1))]
    aCRE_cell1=aCRE_cell1[grep(pro_type2[k],aCRE_cell1[,3]),]
    aCRE_cell1=aCRE_cell1[which(aCRE_cell1[,1]==1),]
    aCRE_cell1$transcription="No"
    aCRE_cell1$transcription[which(aCRE_cell1[,2] == 1)]="Yes"
    aCRE_cell2=reshape2::melt(aCRE_cell1[,c(4:8,10)],id=6)
    aCRE_cell3=aCRE_cell2%>%group_by(variable, transcription, value)%>%dplyr::summarise(count=n())
    aCRE_cell3$id=paste0(aCRE_cell3$transcription,"_",aCRE_cell3$value)
    aCRE_cell4=spread(aCRE_cell3[,c(1,4,5)],key=3,value=2)
    aCRE_cell4$promoter_type_CT=pro_type[k]
    aCRE_cell4$cell=cell[i]
    summaryCRE=rbind(summaryCRE,aCRE_cell4)}}

summaryCRE=summaryCRE[which(!is.na(summaryCRE$No_Yes)),]
for(i in 1:nrow(summaryCRE)){
  GSEATasting <- matrix(c(summaryCRE$Yes_Yes[i], summaryCRE$Yes_No[i], summaryCRE$No_Yes[i], summaryCRE$No_No[i]), nrow = 2, dimnames = list(oligo1 = c("yes", "no"), oligo2 = c("yes", "no")))
  summaryCRE$p.val[i] = fisher.test(GSEATasting, alternative = "two.sided")$p.value
  summaryCRE$OR[i] = fisher.test(GSEATasting, alternative = "two.sided")$estimate}

summaryCRE$variable=sapply(strsplit(as.character(summaryCRE$variable),"\\."),"[",2)
summaryCRE=summaryCRE[-which(summaryCRE$variable == "SE" & summaryCRE$promoter_type_CT == "promoter-like"),]

write.table(summaryCRE, gzfile(paste0(path_fig2_data,"all_aCRE.TSS.access.markeralone_cuttag.FEresult.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
# -> for fig 2d

#===============================================================================
# CRE-TSS, CRE-ATAC and gene signal to metacell UMAP

library(Seurat)
set.seed(2022)
library(Signac)
library(future)

bigDATA=paste0(primary_folder,"Data_and_code/DATA/")

combined=readRDS(paste0(bigDATA,"metacell.combined.feature.Rds"))
metadata=data.frame(combined@meta.data)
umap1=data.frame(combined@reductions[["umap"]]@cell.embeddings)
umap1$cellID=rownames(umap1)
metadata=left_join(metadata, umap1, by=c("SEACell"="cellID"), copy=F)

DefaultAssay(combined) <- "aCRE"
combined <- FindTopFeatures(combined, min.cutoff = 'q75', verbose = TRUE) # use top 25% variable peaks 
combined <- RunSVD(combined, verbose = TRUE)
DepthCor(combined)
combined <- RunUMAP(object = combined, reduction = 'lsi', dims = 2:30) # using first LSI component
DimPlot(object = combined, group.by = "cluster")
umap2=data.frame(combined@reductions[["umap"]]@cell.embeddings)
umap2$cellID=rownames(umap2)
metadata=left_join(metadata, umap2, by=c("SEACell"="cellID"), copy=F, suffix=c("_RNA",""))

DefaultAssay(combined) <- "tCRE"
combined <- FindTopFeatures(combined, min.cutoff = 'q0', verbose = TRUE) # use top 100% variable peaks 
combined <- RunSVD(combined, verbose = TRUE)
DepthCor(combined)
combined <- RunUMAP(object = combined, reduction = 'lsi', dims = 2:30) # using first LSI component
DimPlot(object = combined, group.by = "cluster")
umap3=data.frame(combined@reductions[["umap"]]@cell.embeddings)
umap3$cellID=rownames(umap3)
metadata=left_join(metadata, umap3, by=c("SEACell"="cellID"), copy=F, suffix=c("_aCRE","_tCRE"))
write.table(metadata,gzfile(paste0(CRE_folder,"metacell/meta_metadata.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
#=
cluster_info=read.delim(paste0(primary_folder,"Fig1/cluster_info.tsv"),header=T, stringsAsFactors = F, check.names = F)
metadata=read.delim(paste0(CRE_folder,"metacell/meta_metadata.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

data1=reshape2::melt(metadata[,c(1:23,25,27)], id=c(1:22))
data1$variable=gsub("UMAP_1_","",data1$variable)
data2=reshape2::melt(metadata[,c(1:22,24,26,28)], id=c(1:22))
data2$variable=gsub("UMAP_2_","",data2$variable)
data1=left_join(data1, data2[,c(8,23,24)], by=c("SEACell","variable"),copy=F)
colnames(data1)[c(24,25)]=c("UMAP1","UMAP2")
data1$variable=gsub("RNA","Gene",data1$variable)
data1$variable=gsub("tCRE","TSS-signal",data1$variable)
data1$variable=gsub("aCRE","ATAC-signal",data1$variable)
data1=left_join(data1, cluster_info, by="cluster", copy=F)

write.table(data1,gzfile(paste0(path_fig2_data,"metacell_UMAP.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
# -> for fig ex2d

#===============================================================================
# extract sample based bimodal signals of all CRE
library(Seurat)
combined=readRDS(paste0(bigDATA,"metacell.combined.feature.Rds"))

TSS_table <- as.data.frame(AggregateExpression(combined, group.by = "major_sampleID", assays = "tCRE", slot = "data", return.seurat = FALSE))
ATAC_table <- as.data.frame(AggregateExpression(combined, group.by = "major_sampleID", assays = "aCRE", slot = "data", return.seurat = FALSE))
TSS_table$peakID=rownames(TSS_table)
ATAC_table$peakID=rownames(ATAC_table)
both_table=full_join(ATAC_table,TSS_table, by="peakID", copy=F)
both_table[is.na(both_table)]=0
colnames(both_table)=gsub("iPS","iPSC",colnames(both_table))
colnames(both_table)=gsub("tCRE","TSS",colnames(both_table))
colnames(both_table)=gsub("aCRE","ATAC",colnames(both_table))
write.table(both_table[,c(4,1,3,2,5,7,6)],gzfile(paste0(CRE_folder,"metacell/CRE_463866_sample_expression_ATAC_TSS.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)

aCRE_cell=read.delim(paste0(peak_folder,"all_aCRE.TSS.access.markalone_cuttag.tsv.gz"), header=T, stringsAsFactors = F)
both_table=left_join(both_table, aCRE_cell[,c(1,12,41,49:60,8:10)], by="peakID", copy=F)
both_tablei=both_table[,c(1,8,9,grep("iPSC",colnames(both_table)))]
colnames(both_tablei)=gsub("iPSC","",colnames(both_tablei))
colnames(both_tablei)=gsub("\\.","",colnames(both_tablei))
colnames(both_tablei)=gsub("_","",colnames(both_tablei))
both_tablei$sample="iPSC"
both_tables=both_table[,c(1,8,9,grep("NSC",colnames(both_table)))]
colnames(both_tables)=gsub("NSC","",colnames(both_tables))
colnames(both_tables)=gsub("\\.","",colnames(both_tables))
colnames(both_tables)=gsub("_","",colnames(both_tables))
both_tables$sample="NSC"
both_tablen=both_table[,c(1,8,9,grep("Neuron",colnames(both_table)))]
colnames(both_tablen)=gsub("Neuron","",colnames(both_tablen))
colnames(both_tablen)=gsub("\\.","",colnames(both_tablen))
colnames(both_tablen)=gsub("_","",colnames(both_tablen))
both_tablen$sample="Neuron"
both_table=rbind(both_tablei,both_tables,both_tablen)
mboth_table=reshape2::melt(both_table, id=c(1:5,10,11))
mboth_table=mboth_table[which(mboth_table$value == "Yes"),]
both_table=left_join(both_table,mboth_table[,c("peakID","sample","variable")], by=c("peakID","sample"),copy=F)
both_table=both_table%>%mutate(variable = replace_na(as.character(variable), "Genomic"))
both_table$sample_promoter_type="unclassed"
both_table$sample_promoter_type[grep("enhancer", both_table$CT)]="enhancer-like"
both_table$sample_promoter_type[grep("romoter", both_table$CT)]="promoter-like"
both_table$sample_promoter_type[grep("CTCF", both_table$CT)]="CTCF-alone"
colnames(both_table)[12]="state"
write.table(both_table,gzfile(paste0(path_fig2_data,"CRE_463866_sample_expression_ATAC_TSS.plot.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
# -> for fig ex2f & f2h

#===============================================================================
# correlation between aCRE and tCRE singal
library(Seurat)
combined=readRDS(paste0(bigDATA,"metacell.combined.feature.Rds"))

data1=data.frame(rowSums(data.frame(combined@assays[["tCRE"]]@counts)))
data1$peakID=rownames(data1)
colnames(data1)[1]="tCRE_count"
k3=data1$peakID[which(data1$tCRE_count>0)]

need.atac=which(rownames(combined[["aCRE"]]) %in% k3)
need.rna=which(rownames(combined[["tCRE"]]) %in% k3)
combined[["aCRE"]] <- subset(combined[["aCRE"]], features = rownames(combined[["aCRE"]])[need.atac])
combined[["tCRE"]] <- subset(combined[["tCRE"]], features = rownames(combined[["tCRE"]])[need.rna])
tCRE1=data.frame(combined@assays[["tCRE"]]@data)
aCRE1=data.frame(combined@assays[["aCRE"]]@data)

d <- DGEList(counts=tCRE1)
RLE <- calcNormFactors(d, method="RLE")
RLE.tCRE1=cpm(RLE, normalized.lib.sizes=TRUE)
RLE.tCRE1=as.data.frame(RLE.tCRE1)
d <- DGEList(counts=aCRE1)
RLE <- calcNormFactors(d, method="RLE")
RLE.aCRE1=cpm(RLE, normalized.lib.sizes=TRUE)
RLE.aCRE1=as.data.frame(RLE.aCRE1)

RLE.tCRE1$CREID=rownames(RLE.tCRE1)
RLE.tCRE2=reshape2::melt(RLE.tCRE1,id=392)
RLE.aCRE1$CREID=rownames(RLE.aCRE1)
RLE.aCRE2=reshape2::melt(RLE.aCRE1,id=392)
CRE=left_join(RLE.tCRE2,RLE.aCRE2, by=c("CREID","variable"),suffix=c("_t","_a"))
CRE_rho=CRE%>%group_by(CREID)%>%dplyr::summarise(count=n(),rho=cor.test(value_t,value_a,method="spearman")$estimate,pv=cor.test(value_t,value_a,method="spearman")$p.value)
aCREe=read.delim(paste0(path_fig2_data,"peaks.merged.all.markalone_k16.subclass.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
CRE_rho=left_join(CRE_rho, aCREe[,c(1,9)], by=c("CREID"="V4"),copy=F)

write.table(CRE_rho,gzfile(paste0(path_fig2_data,"CRE_43697_rho_p_e.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)
# -> for fig ex2g

#===============================================================================
#FE on ChromHMM result with CGI/TATA
aCRE_cell=read.delim(paste0(peak_folder,"all_aCRE.TSS.access.markalone_cuttag.tsv.gz"), header=T, stringsAsFactors = F)

# prepare fisher exact between transcription and chromatin state
# divided into cell types, only take cell type specific accessible CREs

cell=c("iPSC","NSC","Neuron")
pro_type=c("enhancer-like","promoter-like")
pro_type2=c("enhancer","romoter")
summaryCRE=data.frame()

for (k in 1:length(pro_type)){
  for (i in 1: length(cell)){
    aCRE_cell1=aCRE_cell[which(aCRE_cell$analysis == "included" & aCRE_cell$promoter_type_CT == pro_type[k]),]
    aCRE_cell1=aCRE_cell1[,c(1,grep(cell[i],colnames(aCRE_cell1)),18,19,26)] # add peakID,CGI, TATA_box and group2
    aCRE_cell1=aCRE_cell1[grep(pro_type2[k],aCRE_cell1[,4]),]
    aCRE_cell1=aCRE_cell1[which(aCRE_cell1[,2]==1),] #only accessible CRE in specific cell-type

    aCRE_cell1_reg=reshape2::melt(aCRE_cell1[,c("peakID","CGI","TATA_box")],id=1)
    aCRE_cell1_state=reshape2::melt(aCRE_cell1[,c(1,13,6:9)],id=c(1,2))
    aCRE_cell1_both=left_join(aCRE_cell1_reg,aCRE_cell1_state, by="peakID",copy=F)
    aCRE_cell1_both2=aCRE_cell1_both%>%group_by(group2,variable.x,variable.y,value.x,value.y)%>%dplyr::summarise(count=n())
    aCRE_cell1_both2$reg_state=paste0(aCRE_cell1_both2$value.x,"_",aCRE_cell1_both2$value.y)
    colnames(aCRE_cell1_both2)[c(2:3)]=c("Reg","State")
    aCRE_cell1_both3=spread(aCRE_cell1_both2[,c(1,2,3,6,7)],key=5, value=4)
    aCRE_cell1_both3=aCRE_cell1_both3[which(!is.na(aCRE_cell1_both3$No_Yes)),]
    aCRE_cell1_both3$Yes_Yes[which(is.na(aCRE_cell1_both3$Yes_Yes))]=0
    aCRE_cell1_both3$promoter_type_CT=pro_type[k]
    aCRE_cell1_both3$cell=cell[i]
    summaryCRE=rbind(summaryCRE,aCRE_cell1_both3)}}

for(i in 1:nrow(summaryCRE)){
  GSEATasting <- matrix(c(summaryCRE$Yes_Yes[i], summaryCRE$Yes_No[i], summaryCRE$No_Yes[i], summaryCRE$No_No[i]), nrow = 2, dimnames = list(oligo1 = c("yes", "no"), oligo2 = c("yes", "no")))
  summaryCRE$p.val[i] = fisher.test(GSEATasting, alternative = "two.sided")$p.value
  summaryCRE$OR[i] = fisher.test(GSEATasting, alternative = "two.sided")$estimate}

summaryCRE$State=sapply(strsplit(as.character(summaryCRE$State),"\\."),"[",2)
write.table(summaryCRE, gzfile(paste0(path_fig2_data,"aCRE_cell_Reg_State_FE.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
# -> for fig 2g & fig ex2k

#===============================================================================
# supp table
aCRE_cell=read.delim(paste0(peak_folder,"all_aCRE.TSS.access.markalone_cuttag.tsv.gz"), header=T, stringsAsFactors = F)
aCRE_cell=aCRE_cell[which(aCRE_cell$analysis == "included"),-c(11,13,14,27:48)]
aCRE_cell$group2=gsub("Without transcription","Non-transcribed",aCRE_cell$group2)
aCRE_cell$group2=gsub("With transcription","Transcribed",aCRE_cell$group2)
aCRE_cell <- aCRE_cell %>% mutate(across(12:22, ~ replace(.x, .x == "Yes", "Y")))
aCRE_cell <- aCRE_cell %>% mutate(across(12:22, ~ replace(.x, .x == "No", "N")))
aCRE_cell <- aCRE_cell %>% mutate(across(2:7, ~ replace(.x, .x == 1, "Y")))
aCRE_cell <- aCRE_cell %>% mutate(across(2:7, ~ replace(.x, .x == 0, "N")))
write.table(aCRE_cell, paste0(primary_folder,"supp_table/all_CRE_features_from_fig2.tsv"), col.names=T, row.names=F, sep="\t", quote=F)







