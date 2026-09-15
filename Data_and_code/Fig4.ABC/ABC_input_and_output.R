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
library(ggrastr)

#===============================================================================
primary_folder="/analysisdata/fantom6/Interactome/single_cell_wallace/Figure/"
path_fig4_data=paste0(primary_folder,"fig4/data/")
ABC_folder=paste0(primary_folder,"Data_and_code/ABC/")

#===============================================================================
# incorporate merged CRE region as enhancer/promoter region
peak_a=read.delim(paste0(primary_folder,"Data_and_code/scATAC/merged_peak/peaks.merged.all.main_chr.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
summary(peak_a$V3-peak_a$V2)

peak_a <- peak_a %>%
  mutate(chr_num = as.numeric(str_extract(V1, "\\d+")),
         chr_num = ifelse(is.na(chr_num), 99, chr_num)) %>% 
  arrange(chr_num, V1, V2) %>%
  select(-chr_num)

write.table(peak_a[,c(1:3)],"/home/yip/tool2026/ABC-Enhancer-Gene-Prediction-main/result_sc/iPSC/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed", row.names=F, col.names=F, sep="\t", quote=F)
write.table(peak_a[,c(1:3)],"/home/yip/tool2026/ABC-Enhancer-Gene-Prediction-main/result_sc/NSC/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed", row.names=F, col.names=F, sep="\t", quote=F)
write.table(peak_a[,c(1:3)],"/home/yip/tool2026/ABC-Enhancer-Gene-Prediction-main/result_sc/Neuron/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed", row.names=F, col.names=F, sep="\t", quote=F)
# put them into ABC pipeline

#===============================================================================
# prepare input: scATAC to tagalign
# sc BAM file from DRA019608 (DRR618513- DRR618515)
# intermediate files too big, not provided.
setwd(paste0(ABC_folder,"input"))
system("sbatch bam_to_tagalign.sh")

#=====
#keep only main chr and resort according to ABC_model requirement
setwd(paste0(ABC_folder,"input/tagalign"))
files=list.files()
files=files[-grep("tbi", files)]

for (i in 1:3){
data=fread(files[i], header=F)
data=data[which(nchar(data$V1)<6),]

data <- data %>%
  mutate(chr_num = as.numeric(str_extract(V1, "\\d+")),
         chr_num = ifelse(is.na(chr_num), 99, chr_num)) %>% 
    arrange(chr_num, V1, V2) %>%
  select(-chr_num)

write.table(data,files[i],row.names=F, col.names=F, sep="\t", quote=F)
}

#===============================================================================
# prepare Hi-C data to bedpe
# HiC data from DRA019572 (DRR614942- DRR614953)
# intermediate files too big, not provided.
setwd(paste0(ABC_folder,"input"))
system("sbatch HiCPE.sbatch")

# split Hi-C files into chromasome base
setwd(paste0(ABC_folder,"input"))
system("sh split_hic_bedpe.sh")

#===============================================================================
# H3K27ac BAM files used as they are
# CUT&Tag data from DRA019568 (DRR614873- DRR614905)

#===============================================================================
# run ABC model separately on iPSC, NSC and Neuron
# run according to the "config.yaml" and "config_biosamples_NeuronSeries_sc.tsv" copied from the config folder

# essential to touch
setwd("/home/yip/tool2026/ABC-Enhancer-Gene-Prediction-main")
system("touch result_sc/iPSC/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed")
system("touch result_sc/NSC/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed")
system("touch result_sc/Neuron/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed")

# -> dry run to confirm successful incorporation

snakemake -n -p \
result_sc/iPSC/Predictions/EnhancerPredictionsAllPutative.tsv.gz \
result_sc/NSC/Predictions/EnhancerPredictionsAllPutative.tsv.gz \
result_sc/Neuron/Predictions/EnhancerPredictionsAllPutative.tsv.gz

#====
#real run
setwd(ABC_folder)
system("sbatch ABC_sc.sh")

#===============================================================================
# gather ABC results
setwd(ABC_folder)
peak_a=read.delim(paste0(primary_folder,"Data_and_code/scATAC/merged_peak/peaks.merged.all.main_chr.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
peak_a$V4=gsub("-","_",peak_a$V4)
colnames(peak_a)[4]="peakID"

path1="/home/yip/tool2026/ABC-Enhancer-Gene-Prediction-main/result_sc/"
files=list.files(path=path1, pattern="EnhancerPredictionsFull_threshold0.02_self_promoter.tsv", recursive =T)
cells=sapply(strsplit(files,"\\/"),"[",1)

data0=data.frame()
for (i in 1:3){
data=fread(paste0(path1,files[i]), header=T)
data=left_join(data, peak_a[,c(1:4)], by=c("chr"="V1", "start"="V2", "end"="V3"),copy=F)
data=data[which(data$class != "promoter"),] #remove self
data0=rbind(data0, data)}
write.table(data0,gzfile("EnhancerPredictionsFull_threshold0.02_3cell_eRNAalone.tsv.gz"), row.names=F, col.names=T, sep="\t", quote=F)

#=====================
#get the overall matrix
setwd(ABC_folder)
path1="/home/yip/tool2026/ABC-Enhancer-Gene-Prediction-main/result_sc/"
files=list.files(path=path1, pattern="EnhancerPredictionsAllPutative.tsv", recursive =T)
cells=sapply(strsplit(files,"\\/"),"[",1)
data0=data.frame()
for (i in 1:3){
data=fread(paste0(path1,files[i]), header=T, select=c(1,2,3,11,28))
data$celltype=cells[i]
data0=rbind(data0, data)}

data0=left_join(data0, peak_a[,c(1:4)], by=c("chr"="V1", "start"="V2", "end"="V3"),copy=F)
data0$pairID=paste0(data0$peakID,"::",data0$TargetGene)
setDT(data0)
data0a <- dcast(data0[, .(id = pairID, key = celltype, value = ABC.Score)], id ~ key, value.var = "value")

data0a$peakID=sapply(strsplit(as.character(data0a$id),"::"),"[",1)
data0a$TargetGene=sapply(strsplit(as.character(data0a$id),"::"),"[",2)
data0a=data0a[which(data0a$NSC+data0a$Neuron+data0a$iPSC!=0),]
data0a=data0a[,c(5,6,4,2,3)]
write.table(data0a,gzfile("ABCscore_matrix_3cells.tsv.gz"), row.names=F, col.names=T, sep="\t", quote=F)

#===============================================================================
# add target gene expression level
data0=read.delim(paste0(ABC_folder,"EnhancerPredictionsFull_threshold0.02_3cell_eRNAalone.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

avg_expr=read.delim(paste0(primary_folder,"Data_and_code/sc5nRNA/output/gene_normalizedExp_perSample.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
avg_expr$geneName=rownames(avg_expr)
colnames(avg_expr)[1]="iPSC"
avg_exprm=reshape2::melt(avg_expr, id=4)
data0=left_join(data0, avg_exprm, by=c("CellType"="variable","TargetGene"="geneName"))
colnames(data0)[34]="target_Exp"

CPM= read.delim(paste0(primary_folder,"Data_and_code/sc5nRNA/output/gene_RLE_CPM_perSample.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
CPM$geneName=rownames(CPM)
colnames(CPM)[1]="iPSC"
CPMm=reshape2::melt(CPM, id=4)
data0=left_join(data0, CPMm, by=c("CellType"="variable","TargetGene"="geneName"))
colnames(data0)[35]="target_CPM"

write.table(data0,gzfile(paste0(ABC_folder,"EnhancerPredictionsFull_threshold0.02_3cell_eRNAalone.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

#===============================================================================
# supplementary table 
# data0=read.delim("EnhancerPredictionsFull_threshold0.02_3cell_eRNAalone.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
# data0a=spread(data0[,c("V4","TargetGene","ex5cluster_class","CpGTATA","n5_string","cell","ABC.Score")], key=6, value=7)
# colnames(data0a)[1]="merged_tCRE_ID"
# data0a=data0a[which(data0a$CpGTATA != "Others"),]

#===============================================================================




