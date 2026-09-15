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
ChIP_folder=paste0(primary_folder,"Data_and_code/ChIP_seq/")
path_fig4_data=paste0(primary_folder,"Fig4/data/")
FP_output=paste0(FP_folder,"output/")
peak_folder=paste0(primary_folder,"Data_and_code/scATAC/merged_peak/")

setwd(ChIP_folder)

#===============================================================================
# mapping by BWA
system("sh mapping.sh")
# fastq and bam are maintained in DDBJ: DRA026820 (DRR959942 - DRR959949)

#===============================================================================
#peak calling
system("sh macs2.sh")
#===
#example
#macs2 callpeak \
#-t TEAD4-rep1.markdup.sort.nodup.bam \
#-c IgG2a-rep1.markdup.sort.nodup.bam \
#-f BAMPE \
#-g hs \
#-n ./macs2/TEAD4-rep1 \
#-B \
#-q 0.01
#===

#===============================================================================
#version 1, merge rep1 & rep2 and "count" the hit to find overlap
setwd(paste0(ChIP_folder,"final_data/TEAD4/"))
system("cat TEAD4-rep1_peaks.narrowPeak TEAD4-rep2_peaks.narrowPeak | sort -k1,1 -k2,2n | bedtools merge -d -1 -i stdin | gzip > TEAD4_merged_peaks.narrowPeak.bed.gz")
system("bedtools intersect -wa -c -a TEAD4_merged_peaks.narrowPeak.bed.gz -b TEAD4-rep1_peaks.narrowPeak | gzip > TEAD4-rep1_count.bed.gz")
system("bedtools intersect -wa -c -a TEAD4_merged_peaks.narrowPeak.bed.gz -b TEAD4-rep2_peaks.narrowPeak | gzip > TEAD4-rep2_count.bed.gz")
setwd(paste0(ChIP_folder,"final_data/ONECUT2/"))
system("cat ONECUT-rep1_peaks.narrowPeak ONECUT-rep2_peaks.narrowPeak | sort -k1,1 -k2,2n | bedtools merge -d -1 -i stdin | gzip > ONECUT_merged_peaks.narrowPeak.bed.gz")
system("bedtools intersect -wa -c -a ONECUT_merged_peaks.narrowPeak.bed.gz -b ONECUT-rep1_peaks.narrowPeak | gzip > ONECUT-rep1_count.bed.gz")
system("bedtools intersect -wa -c -a ONECUT_merged_peaks.narrowPeak.bed.gz -b ONECUT-rep2_peaks.narrowPeak | gzip > ONECUT-rep2_count.bed.gz")

#version 2, idr
system("sh idr.sh")

#===============================================================================
#prepare homer input for merge
setwd(paste0(ChIP_folder,"final_data/TEAD4/"))
TEAD4m=read.delim("TEAD4_merged_peaks.narrowPeak.bed.gz", header=F, stringsAsFactors = F)
TEAD4m1=read.delim("TEAD4-rep1_count.bed.gz", header=F, stringsAsFactors = F)
TEAD4m2=read.delim("TEAD4-rep2_count.bed.gz", header=F, stringsAsFactors = F)
colnames(TEAD4m1)[4]="rep1"
colnames(TEAD4m2)[4]="rep2"
TEAD4m=left_join(TEAD4m,TEAD4m1,by=c("V1","V2","V3"),copy=F)
TEAD4m=left_join(TEAD4m,TEAD4m2,by=c("V1","V2","V3"),copy=F)
TEAD4m=TEAD4m[which(nchar(TEAD4m$V1)<=5),]
TEAD4m$V6="."
write.table(TEAD4m[order(TEAD4m$V1,TEAD4m$V2),c(1,2,3,6,5,6)],"TEAD4_2_output_merge.bed", col.names=F, row.names=F, sep="\t", quote=F)

setwd(paste0(ChIP_folder,"final_data/ONECUT2/"))
ONECUTm=read.delim("ONECUT_merged_peaks.narrowPeak.bed.gz", header=F, stringsAsFactors = F)
ONECUTm1=read.delim("ONECUT-rep1_count.bed.gz", header=F, stringsAsFactors = F)
ONECUTm2=read.delim("ONECUT-rep2_count.bed.gz", header=F, stringsAsFactors = F)
colnames(ONECUTm1)[4]="rep1"
colnames(ONECUTm2)[4]="rep2"
ONECUTm=left_join(ONECUTm,ONECUTm1,by=c("V1","V2","V3"),copy=F)
ONECUTm=left_join(ONECUTm,ONECUTm2,by=c("V1","V2","V3"),copy=F)
ONECUTm=ONECUTm[which(nchar(ONECUTm$V1)<=5),]
ONECUTm$V6="."
write.table(ONECUTm[order(ONECUTm$V1,ONECUTm$V2),c(1,2,3,6,5,6)],"ONECUT_2_output_merge.bed", col.names=F, row.names=F, sep="\t", quote=F)

#homer
setwd(paste0(ChIP_folder,"final_data/TEAD4/"))
system("/home/yip/homer/bin/findMotifsGenome.pl TEAD4_2_output_merge.bed hg38 ./homer_TEAD4_merge -size given")
setwd(paste0(ChIP_folder,"final_data/ONECUT2/"))
system("/home/yip/homer/bin/findMotifsGenome.pl ONECUT_2_output_merge.bed hg38 ./homer_ONECUT_merge -size given")

#===============================================================================
#prepare homer input for idr
setwd(paste0(ChIP_folder,"final_data/TEAD4/idr/"))
TEAD4_idr=read.delim("TEAD4_2c_output.txt", header=F, stringsAsFactors = F)
length(which(TEAD4_idr$V5 >=540)) #2931
TEAD4_idr1=TEAD4_idr[which(TEAD4_idr$V5 >=540),]
TEAD4_idr2=TEAD4_idr1[which(nchar(TEAD4_idr1$V1)<=5),]
write.table(TEAD4_idr2[order(TEAD4_idr2$V1,TEAD4_idr2$V2),],"TEAD4_2_output_p005.bed", col.names=F, row.names=F, sep="\t", quote=F)

setwd(paste0(ChIP_folder,"final_data/ONECUT2/idr/"))
ONECUT_idr=read.delim("ONECUT_2_output.txt", header=F, stringsAsFactors = F)
length(which(ONECUT_idr$V5 >=540)) #804
ONECUT_idr1=ONECUT_idr[which(ONECUT_idr$V5 >=540),]
ONECUT_idr2=ONECUT_idr1[which(nchar(ONECUT_idr1$V1)<=5),]
write.table(ONECUT_idr2[order(ONECUT_idr2$V1,ONECUT_idr2$V2),],"ONECUT_2_output_p005.bed", col.names=F, row.names=F, sep="\t", quote=F)

#homer
setwd(paste0(ChIP_folder,"final_data/TEAD4/idr/"))
system("/home/yip/homer/bin/findMotifsGenome.pl ONECUT_2_output_p005.bed hg38 ./homer_ONECUT_idr -size given")

setwd(paste0(ChIP_folder,"final_data/ONECUT2/idr/"))
system("/home/yip/homer/bin/findMotifsGenome.pl TEAD4_2_output_p005.bed hg38 ./homer_TEAD4_idr -size given")


#===============================================================================
#intersect the chip-seq with aCRE
setwd(paste0(ChIP_folder,"final_data/TEAD4/"))
system(paste0("bedtools intersect -wa -wb -a TEAD4_2_output_merge.bed -b ",primary_folder,"Data_and_code/scATAC/merged_peak/filtered_peaks.merged.all.bed.gz > TEAD4_2_output_merge_aCRE.bed"))
setwd(paste0(ChIP_folder,"final_data/TEAD4/idr/"))
system(paste0("bedtools intersect -wa -wb -a TEAD4_2_output_p005.bed -b ",primary_folder,"Data_and_code/scATAC/merged_peak/filtered_peaks.merged.all.bed.gz > TEAD4_2_output_p005_aCRE.bed"))
#
setwd(paste0(ChIP_folder,"final_data/ONECUT2/"))
system(paste0("bedtools intersect -wa -wb -a ONECUT_2_output_merge.bed -b ",primary_folder,"Data_and_code/scATAC/merged_peak/filtered_peaks.merged.all.bed.gz > ONECUT_2_output_merge_aCRE.bed"))
setwd(paste0(ChIP_folder,"final_data/ONECUT2/idr/"))
system(paste0("bedtools intersect -wa -wb -a ONECUT_2_output_p005.bed -b ",primary_folder,"Data_and_code/scATAC/merged_peak/filtered_peaks.merged.all.bed.gz > ONECUT_2_output_p005_aCRE.bed"))

TEAD4m=read.delim(paste0(ChIP_folder,"final_data/TEAD4/TEAD4_2_output_merge_aCRE.bed"), header=F, stringsAsFactors = F)
TEAD4m1=TEAD4m%>%group_by(V10)%>%dplyr::summarise(ChIP_merge=n())
TEAD4idr=read.delim(paste0(ChIP_folder,"final_data/TEAD4/idr/TEAD4_2_output_p005_aCRE.bed"), header=F, stringsAsFactors = F)
TEAD4idr1=TEAD4idr%>%group_by(V24)%>%dplyr::summarise(ChIP_idr=n())
TEAD4final=full_join(TEAD4m1,TEAD4idr1, by=c("V10"="V24"),copy=F)

TEAD4final=left_join(TEAD4final, TEAD_aCRE, by=c("V10"="ID"), copy=F)
colnames(TEAD4final)[1]="peakID"
write.table(TEAD4final,gzfile(paste0(ChIP_folder,"final_data/TEAD4_ChIP_aCRE_final.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
TEAD4final$chr=sapply(strsplit(TEAD4final$peakID,"-"),"[",1)
TEAD4final$start=sapply(strsplit(TEAD4final$peakID,"-"),"[",2)
TEAD4final$end=sapply(strsplit(TEAD4final$peakID,"-"),"[",3)
write.table(TEAD4final[order(TEAD4final$chr,TEAD4final$start),c(5:7,1)],paste0(ChIP_folder,"final_data/TEAD4/TEAD4_ChIP_aCRE.bed"), col.names=F, row.names=F, sep="\t", quote=F)
#
ONECUTm=read.delim(paste0(ChIP_folder,"final_data/ONECUT2/ONECUT_2_output_merge_aCRE.bed"), header=F, stringsAsFactors = F)
ONECUTm1=ONECUTm%>%group_by(V10)%>%dplyr::summarise(ChIP_merge=n())
ONECUTidr=read.delim(paste0(ChIP_folder,"final_data/ONECUT2/idr/ONECUT_2_output_p005_aCRE.bed"), header=F, stringsAsFactors = F)
ONECUTidr1=ONECUTidr%>%group_by(V24)%>%dplyr::summarise(ChIP_idr=n())
ONECUTfinal=full_join(ONECUTm1,ONECUTidr1, by=c("V10"="V24"),copy=F)

ONECUTfinal=left_join(ONECUTfinal, ONECUT_aCRE, by=c("V10"="ID"), copy=F)
colnames(ONECUTfinal)[1]="peakID"
write.table(ONECUTfinal,gzfile(paste0(ChIP_folder,"final_data/ONECUT_ChIP_aCRE_final.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
ONECUTfinal$chr=sapply(strsplit(ONECUTfinal$peakID,"-"),"[",1)
ONECUTfinal$start=sapply(strsplit(ONECUTfinal$peakID,"-"),"[",2)
ONECUTfinal$end=sapply(strsplit(ONECUTfinal$peakID,"-"),"[",3)
write.table(ONECUTfinal[order(ONECUTfinal$chr,ONECUTfinal$start),c(5:7,1)],paste0(ChIP_folder,"final_data/ONECUT2/ONECUT_ChIP_aCRE.bed"), col.names=F, row.names=F, sep="\t", quote=F)

#===============================================================================
#add back promoter typing by CUT&Tag ChromHMM
aCRE_cell=read.delim(paste0(peak_folder,"all_aCRE.TSS.access.markalone_cuttag.tsv.gz"), header=T, stringsAsFactors = F)
TEAD4final=left_join(TEAD4final, aCRE_cell[,c(1,8:10,12:14)], by="peakID", copy=F)
a1=TEAD4final[which(!is.na(TEAD4final$iPSC_CT)),]%>%group_by(iPSC_CT)%>%dplyr::summarise(chip=n())%>%dplyr::mutate(chip_percent=chip/sum(chip))
b1=aCRE_cell[which(!is.na(aCRE_cell$iPSC_CT)),]%>%group_by(iPSC_CT)%>%dplyr::summarise(bg=n())%>%dplyr::mutate(bg_percent=bg/sum(bg))

ONECUTfinal=left_join(ONECUTfinal, aCRE_cell[,c(1,8:10,12:14)], by="peakID", copy=F)
a2=ONECUTfinal[which(!is.na(ONECUTfinal$Neuron_CT)),]%>%group_by(Neuron_CT)%>%dplyr::summarise(chip=n())%>%dplyr::mutate(chip_percent=chip/sum(chip))
b2=aCRE_cell[which(!is.na(aCRE_cell$Neuron_CT)),]%>%group_by(Neuron_CT)%>%dplyr::summarise(bg=n())%>%dplyr::mutate(bg_percent=bg/sum(bg))
ab1=left_join(a1,b1,by="iPSC_CT")
ab1$rate=ab1$chip_percent/ab1$bg_percent
ab2=left_join(a2,b2,by="Neuron_CT")
ab2$rate=ab2$chip_percent/ab2$bg_percent

# a3=TEAD4final[which(TEAD4final$tCRE=="Yes" & !is.na(TEAD4final$iPSC_CT)),]%>%group_by(iPSC_CT)%>%dplyr::summarise(chip=n())%>%dplyr::mutate(chip_percent=chip/sum(chip))
# b3=aCRE_cell[which(aCRE_cell$tCRE=="Yes" & !is.na(aCRE_cell$iPSC_CT)),]%>%group_by(iPSC_CT)%>%dplyr::summarise(bg=n())%>%dplyr::mutate(bg_percent=bg/sum(bg))
# a4=ONECUTfinal[which(ONECUTfinal$tCRE=="Yes" & !is.na(ONECUTfinal$Neuron_CT)),]%>%group_by(Neuron_CT)%>%dplyr::summarise(chip=n())%>%dplyr::mutate(chip_percent=chip/sum(chip))
# b4=aCRE_cell[which(aCRE_cell$tCRE=="Yes" & !is.na(aCRE_cell$Neuron_CT)),]%>%group_by(Neuron_CT)%>%dplyr::summarise(bg=n())%>%dplyr::mutate(bg_percent=bg/sum(bg))
# ab3=left_join(a3,b3,by="iPSC_CT")
# ab3$rate=ab3$chip_percent/ab3$bg_percent
# ab4=left_join(a4,b4,by="Neuron_CT")
# ab4$rate=ab4$chip_percent/ab4$bg_percent
# ab1$group="All"
# ab2$group="All"
# ab3$group="Transcribed"
# ab4$group="Transcribed"
# ab1=rbind(ab1,ab3)
# ab2=rbind(ab2,ab4)

write.table(TEAD4final,gzfile(paste0(ChIP_folder,"final_data/TEAD4_ChIP_aCRE_final.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
#write.table(ab1,paste0(ChIP_folder,"final_data/TEAD4_ChIP_aCRE_summary.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)

write.table(ONECUTfinal,gzfile(paste0(ChIP_folder,"final_data/ONECUT_ChIP_aCRE_final.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
#write.table(ab2,paste0(ChIP_folder,"final_data/ONECUT_ChIP_aCRE_summary.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)


#===============================================================================
# add motif and tobias result
# consider TEAD[1-4] and ONECUT[1-3]
all_FP=read.delim(paste0(FP_output,"motif_400_CRE463966.matrix.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
all_FP_TEAD=all_FP[,c(1,grep("TEAD",colnames(all_FP)))]
all_FP_ONECUT=all_FP[,c(1,grep("ONECUT",colnames(all_FP)))]
all_FP_TEAD_ID=all_FP_TEAD$peakID[which(rowSums(all_FP_TEAD[,c(2:5)]=="Yes")>0)]
all_FP_ONECUT_ID=all_FP_ONECUT$peakID[which(rowSums(all_FP_ONECUT[,c(2:4)]=="Yes")>0)]

files=list.files(path=FP_output, pattern="FP_")
files_TEAD=files[grep("FP_bound_iPS",files)]
TEAD_FP_bound_ID <- character(0)
for(i in 1:length(files_TEAD)){
  TEAD_FP_bound=read.delim(paste0(FP_output,files_TEAD[i]),header=T)
  TEAD_FP_bound=TEAD_FP_bound[,c(1,grep("TEAD",colnames(TEAD_FP_bound)))]
  TEAD_FP_bound_ID1=TEAD_FP_bound$peakID[which(rowSums(!is.na(TEAD_FP_bound[,c(2:5)]))>0)]
  TEAD_FP_bound_ID=union(TEAD_FP_bound_ID,TEAD_FP_bound_ID1)
}
files_ONECUT=files[grep("FP_bound_neuron",files)]
ONECUT_FP_bound_ID <- character(0)
for(i in 1:length(files_ONECUT)){
  ONECUT_FP_bound=read.delim(paste0(FP_output,files_ONECUT[i]),header=T)
  ONECUT_FP_bound=ONECUT_FP_bound[,c(1,grep("ONECUT",colnames(ONECUT_FP_bound)))]
  ONECUT_FP_bound_ID1=ONECUT_FP_bound$peakID[which(rowSums(!is.na(ONECUT_FP_bound[,c(2:4)]))>0)]
  ONECUT_FP_bound_ID=union(ONECUT_FP_bound_ID,ONECUT_FP_bound_ID1)
}
aCREe=read.delim(paste0(peak_folder,"peaks.merged.all.markalone_k16.subclass.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
aCREe$chr=sapply(strsplit(aCREe$V4,"-"),"[",1)
aCREe=aCREe[which(nchar(aCREe$chr)<=5),] #only main chr

aCREe$TEAD_motif="No"
aCREe$TEAD_motif[which(aCREe$V4 %in% all_FP_TEAD_ID)]="Yes"
aCREe$ONECUT_motif="No"
aCREe$ONECUT_motif[which(aCREe$V4 %in% all_FP_ONECUT_ID)]="Yes"
aCREe$TEAD_FootPrint="No"
aCREe$TEAD_FootPrint[which(aCREe$V4 %in% TEAD_FP_bound_ID)]="Yes"
aCREe$ONECUT_FootPrint="No"
aCREe$ONECUT_FootPrint[which(aCREe$V4 %in% ONECUT_FP_bound_ID)]="Yes"
write.table(aCREe,gzfile(paste0(path_fig4_data,"peaks.merged.all.markalone_k16.subclass_TEAD_ONECUT.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)

#chip
TEAD_chip=read.delim(paste0(ChIP_folder,"final_data/TEAD4_ChIP_aCRE_final.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
TEAD_chip$TEAD_motif="No"
TEAD_chip$TEAD_motif[which(TEAD_chip$peakID %in% all_FP_TEAD_ID)]="Yes"
TEAD_chip$TEAD_FootPrint="No"
TEAD_chip$TEAD_FootPrint[which(TEAD_chip$peakID %in% TEAD_FP_bound_ID)]="Yes"
write.table(TEAD_chip,gzfile(paste0(ChIP_folder,"final_data/TEAD4_ChIP_aCRE_final.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

ONECUT_chip=read.delim(paste0(ChIP_folder,"final_data/ONECUT_ChIP_aCRE_final.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
ONECUT_chip$ONECUT_motif="No"
ONECUT_chip$ONECUT_motif[which(ONECUT_chip$peakID %in% all_FP_ONECUT_ID)]="Yes"
ONECUT_chip$ONECUT_FootPrint="No"
ONECUT_chip$ONECUT_FootPrint[which(ONECUT_chip$peakID %in% ONECUT_FP_bound_ID)]="Yes"
write.table(ONECUT_chip,gzfile(paste0(ChIP_folder,"final_data/ONECUT_ChIP_aCRE_final.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

aCREe$TEAD_ChIP="No"
aCREe$TEAD_ChIP[which(aCREe$V4 %in% TEAD_chip$peakID)]="Yes"
aCREe$ONECUT_ChIP="No"
aCREe$ONECUT_ChIP[which(aCREe$V4 %in% ONECUT_chip$peakID)]="Yes"
write.table(aCREe,gzfile(paste0(path_fig4_data,"peaks.merged.all.markalone_k16.subclass_TEAD_ONECUT.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================











#Code not used
#===============================
# use jasper 2024 now.
# the setting is the same as tobias, no need to re-do
library(motifmatchr)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SummarizedExperiment)
set.seed(2017)
library(rtracklayer)
library(GenomicRanges)

pfm0 <- readJASPARMatrix(paste0(primary_folder,"Data_and_code/ChromVar/jaspar2024/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.adjust_TEAD_ONECUT.txt"), matrixClass = "PFM")
bed0 <- import(paste0(primary_folder,"Data_and_code/scATAC/merged_peak/filtered_peaks.merged.all.bed.gz"))
seqlevelsStyle(bed0) <- "UCSC"
rse <- SummarizedExperiment(rowRanges = bed0)
genome(rse) <- "hg38"
motif_ix <- matchMotifs(pfm0, rse, genome = BSgenome.Hsapiens.UCSC.hg38, p.cutoff = 1e-04)
mm <- as.data.frame(as.matrix(motifMatches(motif_ix)))
gr <- rowRanges(rse)
mm_df <- cbind(
  data.frame(ID = mcols(gr)$name),
  as.matrix(motifMatches(motif_ix)))
motif_id=read.delim(paste0(primary_folder,"Data_and_code/ChromVar/jaspar2024/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.ID_TEAD_ONECUT.tsv"), header=F, stringsAsFactors = F)
colnames(mm_df)[2:8]=motif_id$V2
mm_df1=reshape2::melt(mm_df, id=1)
mm_df1=mm_df1[which(mm_df1$value == "TRUE"),]

TEAD_aCRE=mm_df1[grep("TEAD",mm_df1$variable),]%>%group_by(ID)%>%dplyr::summarise(motifName=paste(variable, collapse=";"))
ONECUT_aCRE=mm_df1[grep("ONECUT",mm_df1$variable),]%>%group_by(ID)%>%dplyr::summarise(motifName=paste(variable, collapse=";"))


#===============================================================================
# define motif scores on several motifs by homer
# use p.cutoff = 1e-04 for all motif, not use this version
setwd(paste0(ChIP_folder,"final_data/TEAD4/"))
system(paste0("/home/yip/homer/bin/findMotifsGenome.pl TEAD4_ChIP_aCRE.bed hg38 ./motiffind -bg ",primary_folder,"Data_and_code/scATAC/merged_peak/filtered_peaks.merged.all.bed.gz -find ./motiffind/jTEAD4_p_top.motif > ./motiffind/homer_TEAD4_motiffind.txt"))
setwd(paste0(ChIP_folder,"final_data/ONECUT2/"))
system(paste0("/home/yip/homer/bin/findMotifsGenome.pl ONECUT_ChIP_aCRE.bed hg38 ./motiffind -bg ",primary_folder,"Data_and_code/scATAC/merged_peak/filtered_peaks.merged.all.bed.gz -find ./motiffind/jONECUT_p_top.motif > ./motiffind/homer_ONECUT2_motiffind.txt"))

setwd(paste0(ChIP_folder,"final_data/TEAD4/motiffind/"))
TEAD4i2=read.delim("homer_TEAD4_motiffind.txt", header=T, stringsAsFactors = F, check.names = T)
TEAD4i2=TEAD4i2%>%group_by(PositionID,Motif.Name)%>%dplyr::slice(which.max(MotifScore))
TEAD4i3=TEAD4i2[,c(1,4,6)]%>%pivot_wider(names_from = Motif.Name, values_from = MotifScore)
TEAD4final=read.delim(paste0(ChIP_folder,"final_data/TEAD4_ChIP_aCRE_final.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
TEAD4final=left_join(TEAD4final,TEAD4i3,by=c("peakID"="PositionID"))
colnames(TEAD4final)[4]="ChromVar_hit"
write.table(TEAD4final,gzfile(paste0(ChIP_folder,"final_data/TEAD4_ChIP_aCRE_final.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
#
setwd(paste0(ChIP_folder,"final_data/ONECUT2/motiffind/"))
ONECUTi2=read.delim("homer_ONECUT2_motiffind.txt", header=T, stringsAsFactors = F, check.names = T)
ONECUTi2=ONECUTi2%>%group_by(PositionID,Motif.Name)%>%dplyr::slice(which.max(MotifScore))
ONECUTi3=ONECUTi2[,c(1,4,6)]%>%pivot_wider(names_from = Motif.Name, values_from = MotifScore)
ONECUTi3=ONECUTi3[,c(1,2,4,7,3,5,6)]
ONECUTfinal=read.delim(paste0(ChIP_folder,"final_data/ONECUT_ChIP_aCRE_final.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
ONECUTfinal=left_join(ONECUTfinal,ONECUTi3,by=c("peakID"="PositionID"))
colnames(ONECUTfinal)[4]="ChromVar_hit"
write.table(ONECUTfinal,gzfile(paste0(ChIP_folder,"final_data/ONECUT_ChIP_aCRE_final.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)


