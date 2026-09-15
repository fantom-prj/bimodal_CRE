library(dplyr) 
library(magrittr)
library(edgeR)
library(knitr)
library(ggplot2)
library(stringr)
library(ggthemes)
library(ggrepel)
library(reshape2)
library(ggsignif)
library(gridExtra)
library(grid)
library(data.table)
library(tidyr)
library(shadowtext)
library(ggrastr)
library(Seurat)

###color#####
library("ggsci")
mypal = pal_npg("nrc", alpha = 1)(10)
mypal2=pal_bmj()(9)
mypal
library("scales")
show_col(mypal)
show_col(mypal2)

#####################
primary_folder="/analysisdata/fantom6/Interactome/single_cell_wallace/Figure/"
FP_folder=paste0(primary_folder,"Data_and_code/Fig4.Foot_Printing/")
FP_output=paste0(FP_folder,"output/")
primary_FP_path="/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/subset_BAM/BINDetect_compare_output/"


setwd(FP_folder)

#===============================================================================
# generate cell cluster specific bam file
#===============================================================================
# subset the scATAC bam for each cell cluster
pca=read.delim(paste0(primary_folder,"Data_and_code/coembed/scDART_n_SEACell/scDART_SEACell_for_all_clusters.tsv"), header=T, stringsAsFactors = F, check.names = F)

# extract cell barcode from each cell cluster
path12="/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/subset_BAM/"
pca=pca[which(pca$cluster != "neuron_progenitor"),]
pca1=pca[which(pca$dataType == "ATAC"),]
pca1$cell=sapply(strsplit(pca1$V1,"_"),"[",1)
pca1$CB=sapply(strsplit(pca1$V1,"_"),"[",2)
cluster=unique(pca1$cluster)
for (i in 1: length(cluster)){
  write.table(pca1[which(pca1$cell == "iPS" & pca1$cluster == cluster[i]),16], paste0(path12,"iPS_",cluster[i],"_CB.csv"), sep=" ", quote=F, row.names=F, col.names=F)
  write.table(pca1[which(pca1$cell == "NSC" & pca1$cluster == cluster[i]),16], paste0(path12,"NSC_",cluster[i],"_CB.csv"), sep=" ", quote=F, row.names=F, col.names=F) 
  write.table(pca1[which(pca1$cell == "NRN" & pca1$cluster == cluster[i]),16], paste0(path12,"NRN_",cluster[i],"_CB.csv"), sep=" ", quote=F, row.names=F, col.names=F) 
}

system("sh subset.bam.sh")
system("sh merge.bam.sh")
# bam files after subset are not provided,
# original bam files please refer to DRA019608 (DRR618513- DRR618515)

#=====================
# run tobias
# ATAC correct
system("sh 20231229tobias.sh")
# FootprintScores
system("sh 20231230tobias_step2.sh")
# BINDetect
system("sh 20240101tobias_step3.sh")

# raw output of tobias footprinting is not provided

#===============================================================================
# summarize tobias results
# extract 400 TF motifs with expression support
together2=read.delim(paste0(primary_folder,"Fig3/data/auto_correlation_to_p_en.z.summarise_noFP.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

need=unique(paste0(together2$motifName,"_",together2$motifID))
need=gsub("::","",need)

# all CRE matrix
scACREv3=read.delim(paste0(primary_folder,"Data_and_code/scATAC/merged_peak/all_aCRE.PE.final.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

scACREv4_motif=scACREv3[,c(1:5)]
scACREv4_FP10p=scACREv4_motif
scACREv4_FP02=scACREv4_motif

for (i in 1:length(need)){
  data=fread(paste0(primary_FP_path,need[i],"/beds/",need[i],"_all.bed"), header=F)
  cutoffs=sapply(data[, c(13:22)], quantile, probs = 0.9, na.rm = TRUE)
  data=data%>%mutate(cutoff10percent = ifelse(rowSums(sweep(data[, c(13:22)], 2, cutoffs, `>=`)) > 0,"Yes","No"))
  data=data%>%mutate(cutoff02 = ifelse(rowSums(sweep(data[, c(13:22)], 2, 0.2, `>=`)) > 0,"Yes","No"))
  scACREv4_motif[[need[i]]] = ifelse(scACREv4_motif$peakID %in% data$V10,"Yes","No")
  scACREv4_FP10p[[need[i]]] = ifelse(scACREv4_FP10p$peakID %in% data$V10[which(data$cutoff10percent == "Yes")],"Yes","No")
  scACREv4_FP02[[need[i]]] = ifelse(scACREv4_FP02$peakID %in% data$V10[which(data$cutoff02 == "Yes")],"Yes","No")
}

scACREv4_motif=scACREv4_motif%>%mutate(anyMotif = ifelse(rowSums(scACREv4_motif[,6:405] == "Yes") > 0, "Yes","No"))
length(which(scACREv4_motif$anyMotif == "No")) #100 -> those not included in main chromosome
scACREv4_FP10p=scACREv4_FP10p%>%mutate(anyFP10p = ifelse(rowSums(scACREv4_FP10p[,6:405] == "Yes") > 0, "Yes","No"))
length(which(scACREv4_FP10p$anyFP10p == "No")) #31806 /463966 = 6.9%
scACREv4_FP02=scACREv4_FP02%>%mutate(anyFP02 = ifelse(rowSums(scACREv4_FP02[,6:405] == "Yes") > 0, "Yes","No"))
length(which(scACREv4_FP02$anyFP02 == "No")) #121138 /463966 = 26.1%

setwd(FP_output)
write.table(scACREv4_motif,gzfile("motif_400_CRE463966.matrix.tsv.gz"),col.names=T, row.names=F, sep="\t", quote=F)
write.table(scACREv4_FP10p,gzfile("Footprint_10percent_cutoff_400_CRE463966.matrix.tsv.gz"),col.names=T, row.names=F, sep="\t", quote=F)
write.table(scACREv4_FP02,gzfile("Footprint_0.2_cutoff_400_CRE463966.matrix.tsv.gz"),col.names=T, row.names=F, sep="\t", quote=F)
rm(scACREv4_motif,scACREv4_FP10p,scACREv4_FP02)

#===============================================================================
# extract "bound" TF from different cluster
scACREv3=scACREv3[,c(1:5)]
scACREv4_bound=scACREv3

files=list.files(path=paste0(primary_FP_path,"AhrArnt_MA0006.2/beds"), pattern="_bound.bed")
cluster=gsub("AhrArnt_MA0006.2_out_","",files)
cluster=gsub("_CB_footprints_bound.bed","",cluster)

for (j in 1:length(cluster)){
for (i in 1:length(need)){
  data=fread(paste0(primary_FP_path,need[i],"/beds/",need[i],"_out_",cluster[j],"_CB_footprints_bound.bed"), select=c(10,13), header=F)
  data=data%>%group_by(V10)%>%slice_max(V13)
  data=data[!duplicated(data$V10),]
  colnames(data)[2]=need[i]
  scACREv4_bound=left_join(scACREv4_bound, data, by=c("peakID"="V10"),copy=F)}
write.table(scACREv4_bound, gzfile(paste0(FP_output,"FP_bound_",cluster[j],"_400_CRE463966.matrix.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)
scACREv4_bound=scACREv3}

