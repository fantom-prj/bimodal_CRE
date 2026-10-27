# Data manipulation
library(dplyr)
library(magrittr)
library(reshape2)

# Plot
library(ggplot2)
library(ggpubr)
library(rstatix)
library(patchwork)

# File I/O
library(readxl)

# Seurat
#library(SeuratDisk)
library(Seurat)

# Signac
set.seed(2022)
library(Signac)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(GenomicRanges)
library(future)

###color#####
library("ggsci")
mypal = pal_npg("nrc", alpha = 1)(9)
mypal
library("scales")
show_col(mypal)

#####################
#primary_folder=[primary_folder]
primary_folder="/analysisdata/fantom6/Interactome/single_cell_wallace/Figure/"
RNA_folder=paste0(primary_folder,"Data_and_code/DATA/")
ATAC_folder=paste0(primary_folder,"Data_and_code/scATAC/")
path_scDART=paste0(primary_folder,"Data_and_code/Fig1.coembed/scDART/")
path_SEACell=paste0(primary_folder,"Data_and_code/Fig1.coembed/SEACell/")
path_fig1_data=paste0(primary_folder,"Fig1/data/")
resource_folder=paste0(primary_folder,"Data_and_code/resources/")
peak_folder=paste0(primary_folder,"Data_and_code/scATAC/merged_peak/")

#===============================================================================
#Put back the PCA and latent into each co-embedded cluster, prepare latent files for python SEACells

All.rna = readRDS(paste0(primary_folder,"Data_and_code/DATA/20220331_All.rna.flt.rds"))
Idents(All.rna) = "celltype_manual"
identity=unique(Idents(All.rna))

###
for (i in 1: length(identity)){
  setwd(path_SEACell)
  system(paste0("mkdir ", path_SEACell, identity[i]))
  
  path1=paste0(path_scDART,identity[i])
  path2=paste0(path_SEACell,identity[i])
  
  setwd(path1)
  
  ATAC.latent=read.csv("ATAC.latent.csv")
  RNA.latent=read.csv("RNA.latent.csv")
  ATAC.PCA=read.csv("ATAC.PCA.csv")
  RNA.PCA=read.csv("RNA.PCA.csv")
  ATAC.label=fread("ATAC.count.csv",header=F, select=1, skip=1)
  RNA.label=fread("RNA.count.csv",header=F, select=1, skip=1)
  
  ATAC.latent=data.frame(cbind(ATAC.label,ATAC.latent[,c(2:9)]))
  ATAC.PCA=data.frame(cbind(ATAC.label,ATAC.PCA[,c(2:9)]))
  RNA.latent=data.frame(cbind(RNA.label,RNA.latent[,c(2:9)]))
  RNA.PCA=data.frame(cbind(RNA.label,RNA.PCA[,c(2:9)]))
  
  rownames(ATAC.latent)=ATAC.latent$V1
  rownames(ATAC.PCA)=ATAC.PCA$V1
  rownames(RNA.latent)=RNA.latent$V1
  rownames(RNA.PCA)=RNA.PCA$V1
  
  ATAC.latent=ATAC.latent[,-1]
  ATAC.PCA=ATAC.PCA[,-1]
  RNA.latent=RNA.latent[,-1]
  RNA.PCA=RNA.PCA[,-1]
  latent=rbind(ATAC.latent,RNA.latent)
  colnames(latent)=c("D1","D2","D3","D4","D5","D6","D7","D8")
  PCA=rbind(ATAC.PCA,RNA.PCA)
  colnames(PCA)=c("P1","P2","P3","P4","P5","P6","P7","P8")
  
  setwd(path2)
  write.csv(PCA, "scDART_PCA.csv", row.names=T)
  write.csv(latent, "scDART_latent.csv", row.names=T)}

#===============================================================================
# schwann cluster is too small, combine progenitor and schwann into one cluster

All.rna = readRDS(paste0(primary_folder,"Data_and_code/DATA/20220331_All.rna.flt.rds"))
Idents(All.rna) = "celltype_manual"

identity="neuron_progenitor_n_schwann"
###
setwd(path_SEACell)
system(paste0("mkdir ", path_SEACell, identity[1]))

path1=paste0(path_scDART,identity[1])
path2=paste0(path_SEACell,identity[1])

setwd(path1)

ATAC.latent=read.csv("ATAC.latent.csv")
RNA.latent=read.csv("RNA.latent.csv")
ATAC.PCA=read.csv("ATAC.PCA.csv")
RNA.PCA=read.csv("RNA.PCA.csv")
ATAC.label=fread("ATAC.count.csv",header=F, select=1, skip=1)
RNA.label=fread("RNA.count.csv",header=F, select=1, skip=1)

ATAC.latent=data.frame(cbind(ATAC.label,ATAC.latent[,c(2:9)]))
ATAC.PCA=data.frame(cbind(ATAC.label,ATAC.PCA[,c(2:9)]))
RNA.latent=data.frame(cbind(RNA.label,RNA.latent[,c(2:9)]))
RNA.PCA=data.frame(cbind(RNA.label,RNA.PCA[,c(2:9)]))

rownames(ATAC.latent)=ATAC.latent$V1
rownames(ATAC.PCA)=ATAC.PCA$V1
rownames(RNA.latent)=RNA.latent$V1
rownames(RNA.PCA)=RNA.PCA$V1

ATAC.latent=ATAC.latent[,-1]
ATAC.PCA=ATAC.PCA[,-1]
RNA.latent=RNA.latent[,-1]
RNA.PCA=RNA.PCA[,-1]
latent=rbind(ATAC.latent,RNA.latent)
colnames(latent)=c("D1","D2","D3","D4","D5","D6","D7","D8")
PCA=rbind(ATAC.PCA,RNA.PCA)
colnames(PCA)=c("P1","P2","P3","P4","P5","P6","P7","P8")

setwd(path2)
write.csv(PCA, "scDART_PCA.csv", row.names=T)
write.csv(latent, "scDART_latent.csv", row.names=T)

#===============================================================================
#run SEACell
#python [cluster_name]_SEACell.py in "Fig1.coembed"

#===============================================================================
#Parse SEACell result
files=list.files(pattern="SEACell_assignment.tsv" , path=path_SEACell, recursive =T)
identity = sapply(strsplit(files,"/"),"[",1)
PCA0=data.frame()
latent0=data.frame()
for (i in 1: length(identity)){
  path2=paste0(path_SEACell,identity[i])
  setwd(path2)
  seacell=read.delim("SEACell_assignment.tsv", header=T, stringsAsFactors =F)
  seacell$metacell=paste0(identity[i],"_",seacell$SEACell)
  seacell$metacell=gsub("SEACell-","meta",seacell$metacell)
  
  PCA=read.csv("scDART_PCA.csv")
  latent=read.csv("scDART_latent.csv")
  
  PCA$dataType="ATAC"
  PCA$dataType[which(nchar(sapply(strsplit(PCA$X,"_"),"[",1))>5)]="RNA"
  latent$dataType="ATAC"
  latent$dataType[which(nchar(sapply(strsplit(latent$X,"_"),"[",1))>5)]="RNA"
  
  PCA=left_join(PCA, seacell[,c(1,3)],by=c("X"="cell_barcode"), copy=F)
  latent=left_join(latent, seacell[,c(1,3)],by=c("X"="cell_barcode"), copy=F)
  colnames(PCA)[11]="SEACell"
  colnames(latent)[11]="SEACell"
  PCA$cluster=identity[i]
  latent$cluster=identity[i]
  latent0=rbind(latent0,latent)
  PCA0=rbind(PCA0, PCA)}

#==
library(uwot)

latent0 <- latent0 %>%
  group_by(cluster) %>%
  group_modify(~{
    x = as.matrix(.x[, grep("^D", colnames(.x))])
    u = uwot::umap(x, n_neighbors = min(15, nrow(x) - 1), min_dist = 0.3, n_components = 2, metric = "euclidean")
    .x$UMAP1 = u[,1]
    .x$UMAP2 = u[,2]
    .x}) %>% ungroup()

PCA0 <- PCA0 %>%
  group_by(cluster) %>%
  group_modify(~{
    x = as.matrix(.x[, grep("^P", colnames(.x))])
    u = uwot::umap(x, n_neighbors = min(15, nrow(x) - 1), min_dist = 0.3, n_components = 2, metric = "euclidean")
    .x$UMAP1 = u[,1]
    .x$UMAP2 = u[,2]
    .x}) %>% ungroup()

write.table(latent0,gzfile(paste0(path_SEACell,"scDART_SEACell_for_all_clusters.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)
write.table(PCA0,gzfile(paste0(path_SEACell,"scDART_SEACell_PCA_for_all_clusters.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

write.table(latent0,gzfile(paste0(path_fig1_data,"scDART_SEACell_for_all_clusters.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)
# -> for fig 1c
# -> for fig ex1b-c
#===============================================================================

