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
###color#####
library("ggsci")
mypal = pal_npg("nrc", alpha = 1)(9)
mypal
library("scales")
show_col(mypal)

# Plot
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

#
#primary_folder=[primary_folder]
primary_folder="/analysisdata/fantom6/Interactome/single_cell_wallace/Figure/"
metacell=paste0(primary_folder,"Data_and_code/Fig1.metacell/")
scDART=paste0(primary_folder,"Data_and_code/Fig1.coembed/scDART_n_SEACell/")
bigDATA=paste0(primary_folder,"Data_and_code/DATA/")
path_fig1_data=paste0(primary_folder,"Fig1/data/")

setwd(metacell)
#=============================================
# before putting together as metacells
#=============================================
# marker genes - per cluster
All.rna.flt=readRDS(paste0(bigDATA,"20220331_All.rna.flt.rds"))
meta0=data.frame(All.rna.flt@meta.data)
meta0$cellID=rownames(meta0)
meta0=left_join(meta0,cluster_info, by=c("celltype_manual"="cluster"),copy=F)
meta0$cluster1[which(is.na(meta0$cluster1))]="Neuron_0"
meta0$label[which(is.na(meta0$label))]="NPC_like"
write.table(meta0,gzfile(paste0(metacell,"output/scRNA_n24434.meta.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

need=unique(meta0$celltype_manual)
cluster_markers = data.frame()
for (i in 1: length(need)){
  cluster_markers1 <- FindMarkers(All.rna.flt, ident.1 = need[i],logfc.threshold = 1,test.use = "wilcox")
  cluster_markers1$celltype_manual=need[i]
  cluster_markers1$geneID=rownames(cluster_markers1)
  cluster_markers=rbind(cluster_markers,cluster_markers1)}
cluster_markers1=cluster_markers[which(cluster_markers$avg_log2FC > 1 & cluster_markers$p_val_adj < 0.01),]
cluster_markers1=left_join(cluster_markers1,cluster_info[,c(1,3)],by=c("celltype_manual"="cluster"),copy=F)
cluster_markers1$label[is.na(cluster_markers1$label)]="NPC_like"
cluster_markers1=cluster_markers1[,c(1:5,8,7)]
write.table(cluster_markers1,gzfile(paste0(metacell,"output/cluster.marker.gene.tsv.gz")),row.names=F, col.names=T, sep="\t", quote=F)

#====
#marker genes for sample type (iPSC, NSC, Neuron)
Idents(All.rna.flt) <- "sampleType"
sample_markers <- FindAllMarkers(object = All.rna.flt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
write.table(sample_markers,gzfile(paste0(metacell,"output/sample.marker.gene.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)
sample_markers1=sample_markers[which(sample_markers$avg_log2FC > 1 & sample_markers$p_val_adj < 0.01),]

#====
#marker genes for mature neuron (both 1 and 2)
Idents(All.rna.flt) <- "celltype_manual"
Mneuron_markers <- FindMarkers(object = All.rna.flt, ident.1 = c("neuron_mature", "neuron_mature_2"), only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Mneuron_markers$cluster="mature_neuron"
Mneuron_markers$gene=rownames(Mneuron_markers)
write.table(Mneuron_markers,gzfile(paste0(metacell,"output/Mature_neuron_2cluster.marker.gene.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)
Mneuron_markers1=Mneuron_markers[which(Mneuron_markers$avg_log2FC > 1 & Mneuron_markers$p_val_adj < 0.01),]

#====================================
# select representative genes
clusters_of_interest=c("iPSC","iPSC","iPSC","iPSC","NSC","NSC","NSC","NPC_like","Immature_neuron","Immature_neuron","Mature_neuron_1","Mature_neuron_1","Mature_neuron_2","Mature_neuron_2","OPC_like")
genes_of_interest=c("POU5F1","EPCAM","UBE2C","HIST1H1E","WLS","TAGLN","CDKN1A","NEUROD1","NEFM","NHLH1","STMN2","RTN1","ONECUT2","HOXB9","DCN")
interest=data.frame(cbind(clusters_of_interest,genes_of_interest))
interest$header=paste0(interest$clusters_of_interest,":\n",interest$genes_of_interest)

need=which(rownames(All.rna.flt) %in% genes_of_interest)
expression_matrix <- GetAssayData(All.rna.flt, slot = "data")
expression_subset <- expression_matrix[need,]
expression_df <- as.data.frame(t(as.matrix(expression_subset))) 
expression_df$cellID=rownames(expression_df)
umap3=read.delim(paste0(metacell,"output/harmony.scRNA.umap.with.metacell.tsv.gz"), header=T, stringsAsFactors =F, check.names = F)
umap3=left_join(umap3,cluster_info, by="cluster", copy=F)
umap3$label=factor(umap3$label, levels=label)
umap3$tech="UMAP of sc5nRNA-seq"
colnames(umap3)[c(4,5)]=c("UMAP1", "UMAP2")
umap3=left_join(umap3[which(umap3$group == "single cell"),], expression_df, by="cellID", copy=F)
umap3=left_join(umap3, meta0[,c(1:35)], by="cellID", copy=F)
write.table(umap3,gzfile(paste0(metacell,"output/scRNA_n24434.with_UMAP_meta.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(umap3,gzfile(paste0(path_fig1_data,"scRNA_n24434.with_UMAP_meta.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#=====
#repeat for ATAC activity
all.Atac=readRDS(paste0(bigDATA,"20220423_All.atac.rds"))
Ameta0=data.frame(all.Atac@meta.data)
Ameta0$cellID=rownames(Ameta0)
DefaultAssay(all.Atac) = "GeneActivity"
Idents(all.Atac) = "celltype_manual_LT"

need=unique(Ameta0$celltype_manual_LT)
cluster_markers = data.frame()
for (i in 1: length(need)){
  cluster_markers1 <- FindMarkers(all.Atac, ident.1 = need[i],logfc.threshold = 0.25,test.use = "wilcox")
  cluster_markers1$celltype_manual_LT=need[i]
  cluster_markers1$geneID=rownames(cluster_markers1)
  cluster_markers=rbind(cluster_markers,cluster_markers1)}
cluster_markers1=cluster_markers[which(cluster_markers$avg_log2FC > 0.25 & cluster_markers$p_val_adj < 0.01),]
write.table(cluster_markers1,gzfile(paste0(metacell,"output/ATAC.cluster.marker.gene.tsv.gz")),row.names=F, col.names=T, sep="\t", quote=F)

expression_matrix <- GetAssayData(all.Atac, assay="GeneActivity", slot = "data")
need=which(rownames(all.Atac) %in% genes_of_interest)
expression_subset <- expression_matrix[need,]
expression_df <- as.data.frame(t(as.matrix(expression_subset))) 
expression_df$cellID=rownames(expression_df)
Aumap3=read.delim(paste0(metacell,"output/harmony.snATAC.umap.with.metacell.tsv.gz"), header=T, stringsAsFactors =F, check.names = F)
Aumap3=left_join(Aumap3,cluster_info, by="cluster", copy=F)
Aumap3$label=factor(Aumap3$label, levels=label)
Aumap3$tech="UMAP of snATAC-seq"
colnames(Aumap3)[c(4,5)]=c("UMAP1", "UMAP2")
Aumap3=left_join(Aumap3[which(Aumap3$group == "single cell"),], expression_df, by="cellID", copy=F)
Aumap3=left_join(Aumap3, Ameta0, by="cellID", copy=F)
write.table(Aumap3,gzfile(paste0(metacell,"output/snATAC_n15048.with_UMAP_meta.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(Aumap3,gzfile(paste0(path_fig1_data,"snATAC_n15048.with_UMAP_meta.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#=====
# repeat for ATAC imputed RNA expression
DefaultAssay(all.Atac) = "RNA"
Idents(all.Atac) <- "celltype_manual_LT"
need=unique(Ameta0$celltype_manual_LT)
cluster_markers = data.frame()

for (i in 1: length(need)){
  cluster_markers1 <- FindMarkers(all.Atac, ident.1 = need[i],logfc.threshold = 0.5,test.use = "wilcox")
  cluster_markers1$celltype_manual_LT=need[i]
  cluster_markers1$geneID=rownames(cluster_markers1)
  cluster_markers=rbind(cluster_markers,cluster_markers1)}
cluster_markers1=cluster_markers[which(cluster_markers$avg_log2FC > 0.5 & cluster_markers$p_val_adj < 0.01),]
write.table(cluster_markers1,gzfile(paste0(metacell,"output/ATAC.cluster.marker.gene_imputed_RNA.tsv.gz")),row.names=F, col.names=T, sep="\t", quote=F)

expression_matrix <- GetAssayData(all.Atac, assay="RNA", slot = "data")
need=which(rownames(expression_matrix) %in% genes_of_interest)
expression_subset <- expression_matrix[need,]
expression_df <- as.data.frame(t(as.matrix(expression_subset))) 
expression_df$cellID=rownames(expression_df)

Aumap3=read.delim(paste0(metacell,"output/harmony.snATAC.umap.with.metacell.tsv.gz"), header=T, stringsAsFactors =F, check.names = F)
Aumap3=left_join(Aumap3,cluster_info, by="cluster", copy=F)
Aumap3$label=factor(Aumap3$label, levels=label)
Aumap3$tech="UMAP of snATAC-seq"
colnames(Aumap3)[c(4,5)]=c("UMAP1", "UMAP2")
Aumap3=left_join(Aumap3[which(Aumap3$group == "single cell"),], expression_df, by="cellID", copy=F)
Aumap3=left_join(Aumap3, Ameta0, by="cellID", copy=F)
write.table(Aumap3,gzfile(paste0(metacell,"output/snATAC_n15048.with_UMAP_meta_imputed_RNA.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(Aumap3,gzfile(paste0(path_fig1_data,"snATAC_n15048.with_UMAP_meta_imputed_RNA.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

# -> for fig 1e & ex1a

#===================
pca=read.delim(paste0(scDART,"scDART_SEACell_for_all_clusters.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
pca=pca[which(pca$cluster != "neuron_progenitor"),]
pca1=pca[,c(1:9)]
umap1=umap(pca1)
pca1=cbind(pca1,umap1)
pca1=left_join(pca1,pca[,c(1,10:12)],by="V1",copy=F)
colnames(pca1)[c(10,11)]=c("UMAP1","UMAP2")
write.table(pca1, gzfile(paste0(scDART,"scDART_SEACell_for_all_clusters.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
#count sum for each metacell
pca1=read.delim(paste0(scDART,"scDART_SEACell_for_all_clusters.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
pca2=pca1%>%group_by(SEACell,dataType)%>%dplyr::summarise(count=n())
pca2=spread(pca2,key=2,value=3)
pca2$number=1:nrow(pca2)
pca1=left_join(pca1,pca2[,c(1,4)],by="SEACell",copy=F)

write.table(pca1,gzfile(paste0(metacell,"output/metadata.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)
write.table(pca2,gzfile(paste0(metacell,"output/metadata2.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

write.table(pca1,gzfile(paste0(path_fig1_data,"metadata.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)
# -> for fig 1c

#===============================================================================
#=============================================
# start putting together as metacells
#=============================================
#===============================================================================
# prepare RNA-seq sc matrix
path2="/analysisdata/fantom6/Interactome/sc5cellranger_Kouno/hg38/"
files1=list.files(path=path2, pattern="barcodes.tsv.gz", recursive=T)
files1=files1[grep("filtered_feature_bc_matrix",files1)]
files2=list.files(path=path2, pattern="features.tsv.gz", recursive=T)
files2=files2[grep("filtered_feature_bc_matrix",files2)]
files3=list.files(path=path2, pattern="matrix.mtx.gz", recursive=T)
files3=files3[grep("filtered_feature_bc_matrix",files3)]

i=1
barcodes=fread(paste0(path2,files1[i]), stringsAsFactors = F, header=F)
barcodes$V2=1:nrow(barcodes)
barcodes$V1=paste0(barcodes$V1,"_",i)
barcodes=left_join(barcodes, pca1[,c(1,15)], by="V1", copy=F)
barcodes2=barcodes[which(!is.na(barcodes$number)),]

m1=read.delim(paste0(path2,files3[i]), stringsAsFactors = F, header=F, skip=3, sep=" ")
m1=left_join(m1,barcodes2[,c(2,3)], by="V2",copy=F)
m2=m1[which(!is.na(m1$number)),]
m3=m2%>%group_by(V1,number)%>%dplyr::summarise(V3=sum(V3))

for (i in 2:length(files1)){
  barcodes=fread(paste0(path2,files1[i]), stringsAsFactors = F, header=F)
  barcodes$V2=1:nrow(barcodes)
  barcodes$V1=paste0(barcodes$V1,"_",i)
  barcodes=left_join(barcodes, pca1[,c(1,15)], by="V1", copy=F)
  barcodes2=barcodes[which(!is.na(barcodes$number)),]
  
  m1=read.delim(paste0(path2,files3[i]), stringsAsFactors = F, header=F, skip=3, sep=" ")
  m1=left_join(m1,barcodes2[,c(2,3)], by="V2",copy=F)
  m2=m1[which(!is.na(m1$number)),]
  im3=m2%>%group_by(V1,number)%>%dplyr::summarise(V3=sum(V3))
  m3=rbind(m3,im3)}

setwd(paste0(metacell,"scRNA"))
m4=m3%>%group_by(V1,number)%>%dplyr::summarise(V3=sum(V3))
write.table(m4,"matrix.mtx",col.names=F,row.names=F, sep=" ",quote=F)
m=read.delim(paste0(path2,files3[1]), stringsAsFactors = F, header=F, nrows=3)
m$V1[3]=paste0("36601 ",nrow(pca2), " ", nrow(m4))

write.table(m,"matrix.mtx.head",col.names=F,row.names=F, sep=" ",quote=F)
genes=read.delim(paste0(path2,files2[1]), stringsAsFactors = F, header=F)
write.table(genes,"genes.tsv",col.names=F,row.names=F, sep="\t",quote=F)
write.table(pca2[,1],"barcodes.tsv", col.names=F, row.names=F, sep="\t", quote=F)

#===============================================================================
# prepare TSS sc matrix
SCAFE=paste0(primary_folder,"/Data_and_code/SCAFE/sc5end/count/aCRE_original/")
files1=list.files(path=SCAFE, pattern="barcodes.tsv", recursive=T)
files2=list.files(path=SCAFE, pattern="genes.tsv", recursive=T)
files3=list.files(path=SCAFE, pattern="matrix.mtx", recursive=T)
files.names=sapply(strsplit(files1,"\\/"),"[",1)

i=1
barcodes=fread(paste0(SCAFE,files1[i]), stringsAsFactors = F, header=F)
barcodes$V2=1:nrow(barcodes)
barcodes$V1=paste0(barcodes$V1,"_",i)
barcodes=left_join(barcodes, pca1[,c(1,15)], by="V1", copy=F)
barcodes2=barcodes[which(!is.na(barcodes$number)),]

m1=read.delim(paste0(SCAFE,files3[i]), stringsAsFactors = F, header=F, skip=3, sep=" ")
m1=left_join(m1,barcodes2[,c(2,3)], by="V2",copy=F)
m2=m1[which(!is.na(m1$number)),]
m3=m2%>%group_by(V1,number)%>%dplyr::summarise(V3=sum(V3))

for (i in 2:length(files1)){
  barcodes=fread(paste0(SCAFE,files1[i]), stringsAsFactors = F, header=F)
  barcodes$V2=1:nrow(barcodes)
  barcodes$V1=paste0(barcodes$V1,"_",i)
  barcodes=left_join(barcodes, pca1[,c(1,15)], by="V1", copy=F)
  barcodes2=barcodes[which(!is.na(barcodes$number)),]
  
  m1=read.delim(paste0(SCAFE,files3[i]), stringsAsFactors = F, header=F, skip=3, sep=" ")
  m1=left_join(m1,barcodes2[,c(2,3)], by="V2",copy=F)
  m2=m1[which(!is.na(m1$number)),]
  im3=m2%>%group_by(V1,number)%>%dplyr::summarise(V3=sum(V3))
  m3=rbind(m3,im3)}

setwd(paste0(metacell,"tCRE_ori_ATAC_region"))
m4=m3%>%group_by(V1,number)%>%dplyr::summarise(V3=sum(V3))
write.table(m4,"matrix.mtx",col.names=F,row.names=F, sep=" ",quote=F)

m=read.delim(paste0(SCAFE,files3[1]), stringsAsFactors = F, header=F, nrows=3)
m$V1[3]=paste0("463966 ",nrow(pca2), " ", nrow(m4))
write.table(m,"matrix.mtx.head",col.names=F,row.names=F, sep=" ",quote=F)
genes=read.delim(paste0(SCAFE,files2[1]), stringsAsFactors = F, header=F)
write.table(genes,"genes.tsv",col.names=F,row.names=F, sep="\t",quote=F)
write.table(pca2[,1],"barcodes.tsv", col.names=F, row.names=F, sep="\t", quote=F)

#===============================================================================
# ATAC sc matrix preparation
All.atac.raw <- readRDS(paste0(bigDATA,"All.atac.raw.3samples.rds"))
data = as.data.frame(summary(All.atac.raw@assays[["peaks"]]@counts))
peaks = data.frame(cbind(rownames(All.atac.raw),rownames(All.atac.raw)))
barcode = data.frame(colnames(All.atac.raw))
barcode$V2=1:nrow(barcode)
colnames(barcode)[1]="V1"
barcode=left_join(barcode, pca1[,c(1,15)], by="V1", copy=F)
barcode2=barcode[which(!is.na(barcode$number)),]

setwd(paste0(metacell,"scATAC"))
data=left_join(data,barcode2[,c(2,3)], by=c("j"="V2"),copy=F)
data2=data[which(!is.na(data$number)),]
data3=data2%>%group_by(i,number)%>%dplyr::summarise(x=sum(x))
write.table(data3,"matrix.mtx",col.names=F,row.names=F, sep=" ",quote=F)

m=read.delim(paste0(path2,files3[1]), stringsAsFactors = F, header=F, nrows=3)
m$V1[3]=paste0(nrow(peaks)," ",nrow(pca2), " ", nrow(data3))
write.table(m,"matrix.mtx.head",col.names=F,row.names=F, sep=" ",quote=F)
write.table(peaks,"genes.tsv",col.names=F,row.names=F, sep="\t",quote=F)
write.table(pca2[,1],"barcodes.tsv", col.names=F, row.names=F, sep="\t", quote=F)

#=============================================================
# put together ATAC and TSS signal per CRE
combined <- Read10X(data.dir = paste0(metacell,"scRNA"))
combined <- CreateSeuratObject(counts = combined, project = "meta_cell", min.cells = 0, min.features = 0)

#add tCRE
tCRE <- Read10X(data.dir = paste0(metacell,"tCRE_ori_ATAC_region"))
tCRE.assay <- CreateChromatinAssay(counts = tCRE, sep = c("-", "-"))
combined[["tCRE"]] <- tCRE.assay

#add aCRE
aCRE <- Read10X(data.dir = paste0(metacell,"scATAC"))
aCRE.assay <- CreateChromatinAssay(counts = aCRE, sep = c("-", "-"))
combined[["aCRE"]] <- aCRE.assay

#add metadata
metadata=data.frame(combined@meta.data)
md1=read.delim(paste0(metacell,"output/metadata2.tsv"), header=T, stringsAsFactors = F, check.names = F)
pca1=read.delim(paste0(metacell,"output/metadata.tsv"),header=T, stringsAsFactors = F, check.names = F)
md1=left_join(md1,unique(pca1[,c(13,14)]),by="SEACell",copy=F)
pca1$sampleID="iPS"
pca1$sampleID[grep("NSC_",pca1$V1)]="NSC"
pca1$sampleID[grep("_5",pca1$V1)]="NSC"
pca1$sampleID[grep("_6",pca1$V1)]="NSC"
pca1$sampleID[grep("NRN_",pca1$V1)]="Neuron"
pca1$sampleID[grep("_3",pca1$V1)]="Neuron"
pca1$sampleID[grep("_4",pca1$V1)]="Neuron"
pca2=pca1%>%group_by(SEACell,sampleID)%>%dplyr::summarise(count=n())
pca2a=pca2%>%group_by(SEACell)%>%dplyr::slice_max(count)
colnames(pca2a)[2]="major_sampleID"
pca3=spread(pca2,key=2, value=3)
pca3[is.na(pca3)]=0
md1=left_join(md1,pca3,by="SEACell",copy=F)
md1=left_join(md1,pca2a[,c(1,2)],by="SEACell",copy=F)

colnames(md1)[c(2:4)]=c("ATAC_cell","RNA_cell","Total_cell")
md1$Total_cell=md1$ATAC_cell+md1$RNA_cell
metadata=cbind(metadata,md1)
combined@meta.data=metadata

DefaultAssay(combined) <- "RNA"
#basic
Sys.time()
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined)
combined <- RunUMAP(combined, dims = 1:30)
combined <- FindNeighbors(combined, dims = 1:30)
combined <- FindClusters(combined,  algorithm = 1 )
combined$log1p_Louvain_res.0.8 <-  combined$seurat_clusters
combined$RNA_snn_res.0.8 <- NULL
combined <- FindClusters(combined, algorithm = 4)
combined$log1p_Leiden_res.0.8 <-  combined$seurat_clusters
combined$RNA_snn_res.0.8 <- NULL
Sys.time() 

#cell cycle scoring
data("cc.genes")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

combined <- CellCycleScoring(combined, g2m.features = g2m.genes, s.features = s.genes)
combined$Phase <- factor(combined$Phase)

#===============================================================================
# select 391 metacells
combined=subset(x = combined, subset = ATAC_cell >= 3 & RNA_cell >=3)
saveRDS(combined, paste0(bigDATA,"metacell.combined.feature.Rds"))

#=================================================================
#define which CREs are transcribed, read > 0 in any sample
combined=readRDS(paste0(bigDATA,"metacell.combined.feature.Rds"))
data1=data.frame(rowSums(data.frame(combined@assays[["tCRE"]]@counts)))
data1$peakID=rownames(data1)
colnames(data1)[1]="tCRE_count"
write.table(data1[,c(2,1)],gzfile(paste0(metacell,"output/tCRE_count_in_all_CRE.tsv.gz")), col.names=T, row.names=F, sep="\t",quote=F)

#===============================================================================
#manually curate TF function
setwd(paste0(metacell,"output/TF_motif/"))
motif.df=read.delim(paste0(primary_folder,"Data_and_code/ChromVar/jaspar2024/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.ID.tsv"), header=F, stringsAsFactors = F)

motif.df1=motif.df
motif.df1$motifName1=sapply(strsplit(motif.df1$V2,"::"),"[",1)
motif.df1$motifName2=sapply(strsplit(motif.df1$V2,"::"),"[",2)
motif.df2=melt(motif.df1, id=c(1,2))
motif.df2=motif.df2[which(!is.na(motif.df2$value)),]
motif.df2$value=sapply(strsplit(motif.df2$value,"\\."),"[",1)
write.table(unique(motif.df2[,c(4)]),"all.motif.name.query.tsv", row.names=F, col.names=F, sep="\t", quote=F)

name.check=read.csv("hgnc-symbol-check_2024.csv", header=T, skip=1)
name.check=name.check[which(name.check$Match.type == "Approved symbol"),]
motif.df2=left_join(motif.df2,unique(name.check), by=c("value"="Input"), copy=F)
motif.df2$Approved.symbol[which(is.na(motif.df2$Approved.symbol))]=motif.df2$value[which(is.na(motif.df2$Approved.symbol))]
motif.df2$Approved.symbol=toupper(motif.df2$Approved.symbol)
write.table(unique(motif.df2$Approved.symbol),"all_TF_HGNC.tsv", col.names = F, row.names=F, sep="\t", quote=F)

# id_mapping_tools from uniprot
function1=read.delim("idmapping_2024_04_25.tsv", header=T, stringsAsFactors = F, check.names = F)
length(which(function1$`Function [CC]` != ""))
function1=function1[which(function1$`Function [CC]` != ""),]
function2=function1%>%group_by(From)%>%dplyr::slice_max(Annotation)
function2=function2[which(function2$Reviewed == "reviewed"),]
function2$TF_activity=""
function2$TF_activity[grep("repress",function2$`Function [CC]`)]="repressor"
write.table(function2, "idmapping_2024_04_25.tsv", col.names=T, row.names=F, sep="\t", quote=F)

# after manual curation
function2=read.delim("idmapping_2024_04_25.tsv",header=T, stringsAsFactors = F)
function2$curated="no"
function2$curated[c(1:255,570)]="yes"
function2$TF_activity2=""
function2$TF_activity2[grep("activate",function2$Function..CC.)]="activator"
function2$TF_activity2[grep("activated",function2$Function..CC.)]=""
function2$TF_activity2[grep("activator",function2$Function..CC.)]="activator"
function2$TF_activity[which(function2$curated=="no")]=function2$TF_activity2[which(function2$curated=="no")]
function2$TF_activity[which(function2$TF_activity == "unknown" | function2$TF_activity == "uknown" | function2$TF_activity=="")]="unknown"
function2%>%group_by(TF_activity)%>%dplyr::summarise(count=n())
write.table(function2, "idmapping_2024_04_25.tsv", col.names=T, row.names=F, sep="\t", quote=F)
function4=function2%>%group_by(From)%>%dplyr::slice_min(TF_activity)
function4=function4[-which(function4$Entry=="Q14494" & function4$From == "NRF1"),]
motif.df2=left_join(motif.df2, function4, by=c("Approved.symbol"="From"), copy=F)

expression=data.frame(t(combined@assays[["RNA"]]@data))
expression1=expression[,which(colnames(expression)%in% unique(motif.df2$Approved.symbol))]
expression1$metacell=rownames(expression1)
expression2=melt(expression1, id=736)
write.table(expression2,gzfile("TF.expression.per.meta.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)
expression3=expression2%>%group_by(variable)%>%dplyr::summarise(max=max(value),sum=sum(value),mean=mean(value),count=sum(value!=0))

motif.df4=left_join(motif.df2, expression3, by=c("Approved.symbol"="variable"), copy=F)
write.table(motif.df4,"motif.info.with_TF.expression.tsv", col.names=T, row.names=F, sep="\t", quote=F)
motif.df4rep=motif.df4[grep("repress",motif.df4$TF_activity),]
motif.df4rep=motif.df4rep[which(motif.df4rep$max >=0.2),]
motif.df4rep=motif.df4rep%>%group_by(V1)%>%dplyr::mutate(nTF=n())
motifout=motif.df4rep$V1[intersect(grep(":",motif.df4rep$V2), which(motif.df4rep$nTF ==1))]
motif.df4rep=motif.df4rep[-which(motif.df4rep$V1 %in% motifout),]
length(unique(motif.df4rep$V1)) #145
write.table(motif.df4rep,"motif.info.with_TF.expression.repress.max02.tsv", col.names=T, row.names=F, sep="\t", quote=F)

motif.df1=left_join(motif.df1,unique(motif.df2[,c(4,6)]), by=c("motifName1"="value"))
motif.df1=left_join(motif.df1,unique(motif.df2[,c(4,6)]), by=c("motifName2"="value"), suffix=c(".TF1",".TF2"))

motif.df5=left_join(motif.df1, expression2, by=c("Approved.symbol.TF1"="variable"), copy=F)
motif.df5=left_join(motif.df5, expression2, by=c("Approved.symbol.TF2"="variable", "metacell"="metacell"), copy=F, suffix=c(".TF1",".TF2"))
colnames(motif.df5)[c(1:2)]=c("motifID","motifName")
write.table(motif.df5,gzfile("motif.base.TF.expression.per.meta.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)
# -> expression of TFs from all 879 motifs

#===============================================================================
#filter motif for expression, max 0.2 for 3 metacells --> 400 motif
motif.df5=read.delim("motif.base.TF.expression.per.meta.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
motif.df6a=motif.df5[which(is.na(motif.df5$motifName2)),]
motif.df6b=motif.df5[which(!is.na(motif.df5$motifName2)),]
motif.df6a_summary=motif.df6a[which(motif.df6a$value.TF1 >= 0.2),]%>%group_by(motifID)%>%dplyr::summarise(count=n())
need1=motif.df6a_summary$motifID[which(motif.df6a_summary$count>=3)] #377
motif.df6b_summary=motif.df6b[which(motif.df6b$value.TF1 >= 0.2 & motif.df6b$value.TF2 >=0.2),]%>%group_by(motifID)%>%dplyr::summarise(count=n())
need2=motif.df6b_summary$motifID[which(motif.df6b_summary$count>=3)] #23
length(union(need1,need2)) #400
motif.df6=motif.df5[which(motif.df5$motifID %in% union(need1,need2)),]
write.table(motif.df6,gzfile("motif.base.TF.expression.per.max02.cell3.meta.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)

# motif with TF expression and literature base function is stored in paste0(metacell,"output/TF_motif/")

#=====================================================================================================
#take the single cell pseudotime and take the mean for metacell
pseudotime=read.delim(paste0(primary_folder,"Data_and_code/sc5nRNA/output/pseudotime_table_scRNA.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
pseudotime=left_join(pseudotime[,c(1,37:41)],pca1[,c(1,13,14)],by=c("cellID" = "V1"),copy=F)
pseudotime1=pseudotime%>%group_by(SEACell)%>%dplyr::summarise(pseudotime=mean(pseudotime))
pseudotime1=pseudotime1[order(pseudotime1$pseudotime),]
pseudotime1$T1=1:nrow(pseudotime1)
write.table(pseudotime1, gzfile(paste0(metacell,"output/pseudotime.scRNA.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

#=====================================================================================================
# take the scRNA UMAP and transfer to metacell
cluster1=c("iPS_0_undifferentiaed","iPS_1_neuroepithelial","NSC_0_iPS_like", "NSC_1_differentiating", "NSC_2_astrocyte_like", "neuron_progenitor_n_schwann", "neuron_immature", "neuron_mature","neuron_mature_2", "OPC_like")
All.RNA <- readRDS("/analysisdata/fantom6/Interactome/scRNA_scATAC_analysis_JC/20220331_All.rna.flt.rds")
umap1=data.frame(All.RNA@reductions[["umap_harmony"]]@cell.embeddings)
umap1$cellID=rownames(umap1)
umap1=left_join(umap1,pca1[,c(1,13,14)],by=c("cellID"="V1"),copy=F)
umap2=umap1%>%group_by(SEACell,cluster)%>%dplyr::summarise(umapharmony_1=mean(umapharmony_1),umapharmony_2=mean(umapharmony_2))
umap2$group="meta cell"
umap1$group="single cell"
colnames(umap2)[1]="cellID"
umap3=rbind(umap1[,c(3,5,6,1,2)],umap2[,c(1,2,5,3,4)])
write.table(umap3,gzfile(paste0(metacell,"output/harmony.scRNA.umap.with.metacell.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

write.table(umap3,gzfile(paste0(path_fig1_data,"harmony.scRNA.umap.with.metacell.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
# ->for fig 1b

# take the snATAC UMAP and transfer to metacell
All.atac <- readRDS("/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/All.atac.rds")
Aumap1=data.frame(All.atac@reductions[["umap_integrated"]]@cell.embeddings)
Aumap1$cellID=rownames(Aumap1)
Aumap1=left_join(Aumap1,pca1[,c(1,13,14)],by=c("cellID"="V1"),copy=F)
Aumap2=Aumap1%>%group_by(SEACell,cluster)%>%dplyr::summarise(UMAP_1=mean(UMAP_1),UMAP_2=mean(UMAP_2))
Aumap2$group="meta cell"
Aumap1$group="single cell"
colnames(Aumap2)[1]="cellID"
Aumap3=rbind(Aumap1[,c(3,5,6,1,2)],Aumap2[,c(1,2,5,3,4)])
Aumap3$cluster=factor(Aumap3$cluster, levels=cluster1)
write.table(Aumap3,gzfile(paste0(metacell,"output/harmony.snATAC.umap.with.metacell.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

write.table(Aumap3,gzfile(paste0(path_fig1_data,"harmony.snATAC.umap.with.metacell.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
# ->for fig 1b

#=============================
#get normalized count for aCRE and tCRE
combined=readRDS(paste0(bigDATA,"metacell.combined.feature.Rds"))
combined <- NormalizeData(
  object = combined,
  assay = "aCRE", 
  normalization.method = "LogNormalize", 
  scale.factor = 10000)
combined <- NormalizeData(
  object = combined,
  assay = "tCRE", 
  normalization.method = "LogNormalize", 
  scale.factor = 10000)
raw_counts <- GetAssayData(combined, assay = "aCRE", slot = "counts")
normalized_counts <- GetAssayData(combined, assay = "aCRE", slot = "data")
all.equal(raw_counts, normalized_counts)
saveRDS(combined, paste0(bigDATA,"metacell.combined.feature.Rds"))

#===============================================================================
# expression level for each cell cluster
combined=readRDS(paste0(bigDATA,"metacell.combined.feature.Rds"))
avg_exp = AverageExpression(combined, group.by = "cluster",  assays = "RNA")
exp=data.frame(avg_exp$RNA)
exp$geneName=rownames(exp)
gtf=read.delim(paste0(primary_folder,"Data_and_code/resources/cellranger_genes.gtf.gz"), skip=5, header=F, stringsAsFactors = F)
gtf=gtf[which(gtf$V3 == "gene"),]
gtf$geneID=sapply(strsplit(gtf$V9,"; "),"[",1)
gtf$geneName=sapply(strsplit(gtf$V9,"; "),"[",4)
gtf$geneID=gsub("gene_id ","",gtf$geneID)
gtf$geneName=gsub("gene_name ","",gtf$geneName)
exp=left_join(exp,gtf[,c(10:11)],by="geneName",copy=F)
write.table(exp, gzfile(paste0(metacell,"output/RNA_expression_per_cluster.matrix.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

ATACm = as.data.frame(AverageExpression(object = combined, group.by = "cluster", assays = "aCRE", slot = "data")$aCRE)
ATACm$peakID=rownames(ATACm)
TSSm = as.data.frame(AverageExpression(object = combined, group.by = "cluster", assays = "tCRE", slot = "data")$tCRE)
TSSm$peakID=rownames(TSSm)
write.table(ATACm, gzfile(paste0(metacell,"output/CRE_ATAC_per_cluster.matrix.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(TSSm, gzfile(paste0(metacell,"output/CRE_TSS_per_cluster.matrix.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================


