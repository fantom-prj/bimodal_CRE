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
path_SEACell=paste0(primary_folder,"Data_and_code/Fig1.coembed/SEACell/")
bigDATA=paste0(primary_folder,"Data_and_code/DATA/")
path_fig1_data=paste0(primary_folder,"Fig1/data/")
path_fig2_data=paste0(primary_folder,"Fig2/data/")

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

#===============================================================================
#count sum for each metacell
latent0=read.delim(paste0(path_SEACell,"scDART_SEACell_for_all_clusters.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
latent1=latent0%>%group_by(SEACell,dataType)%>%dplyr::summarise(count=n())
latent1=spread(latent1,key=2,value=3)
latent1$ATAC[is.na(latent1$ATAC)]=0
latent1$number=1:nrow(latent1)
latent0=left_join(latent0,latent1[,c(1,4)],by="SEACell",copy=F)

#keep all the 395 metacells first, filter later
write.table(latent0,gzfile(paste0(metacell,"output/metadata.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)
write.table(latent1,gzfile(paste0(metacell,"output/metadata2.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
#=============================================
# start putting together as metacells
#=============================================
#===============================================================================
# prepare RNA-seq sc matrix -> meta matrix
path2=paste0(primary_folder,"Data_and_code/sc5nRNA/cellranger/")
files1=list.files(path=path2, pattern="barcodes.tsv.gz", recursive=T)
files1=files1[grep("filtered_feature_bc_matrix",files1)]
files2=list.files(path=path2, pattern="features.tsv.gz", recursive=T)
files2=files2[grep("filtered_feature_bc_matrix",files2)]
files3=list.files(path=path2, pattern="matrix.mtx.gz", recursive=T)
files3=files3[grep("filtered_feature_bc_matrix",files3)]

m3=data.frame()
for (i in 1:length(files1)){
  barcodes=fread(paste0(path2,files1[i]), stringsAsFactors = F, header=F)
  barcodes$V2=1:nrow(barcodes)
  barcodes$V1=paste0(barcodes$V1,"_",i)
  barcodes=left_join(barcodes, latent0[,c(2,15)], by=c("V1"="X"), copy=F)
  barcodes2=barcodes[which(!is.na(barcodes$number)),]
  
  m1=read.delim(paste0(path2,files3[i]), stringsAsFactors = F, header=F, skip=3, sep=" ")
  m1=left_join(m1,barcodes2[,c(2,3)], by="V2",copy=F)
  m2=m1[which(!is.na(m1$number)),]
  im3=m2%>%group_by(V1,number)%>%dplyr::summarise(V3=sum(V3))
  m3=rbind(m3,im3)}

m4=m3%>%group_by(V1,number)%>%dplyr::summarise(V3=sum(V3))
write.table(m4,paste0(metacell,"scRNA/matrix.mtx"),col.names=F,row.names=F, sep=" ",quote=F)

m=read.delim(paste0(path2,files3[1]), stringsAsFactors = F, header=F, nrows=3)
m$V1[3]=paste0("36601 ",nrow(latent1), " ", nrow(m4))
write.table(m,paste0(metacell,"scRNA/matrix.mtx.head"),col.names=F,row.names=F, sep=" ",quote=F)

genes=read.delim(paste0(path2,files2[1]), stringsAsFactors = F, header=F)
write.table(genes,paste0(metacell,"scRNA/genes.tsv"),col.names=F,row.names=F, sep="\t",quote=F)
write.table(latent1[,1],paste0(metacell,"scRNA/barcodes.tsv"), col.names=F, row.names=F, sep="\t", quote=F)

#===============================================================================
# prepare TSS sc matrix -> metacell matrix
SCAFE=paste0(primary_folder,"/Data_and_code/SCAFE/sc5end/count/aCRE_original/")
files1=list.files(path=SCAFE, pattern="barcodes.tsv", recursive=T)
files2=list.files(path=SCAFE, pattern="genes.tsv", recursive=T)
files3=list.files(path=SCAFE, pattern="matrix.mtx", recursive=T)
files.names=sapply(strsplit(files1,"\\/"),"[",1)

m3=data.frame()
for (i in 1:length(files1)){
  barcodes=fread(paste0(SCAFE,files1[i]), stringsAsFactors = F, header=F)
  barcodes$V2=1:nrow(barcodes)
  barcodes$V1=paste0(barcodes$V1,"_",i)
  barcodes=left_join(barcodes, latent0[,c(2,15)], by=c("V1"="X"), copy=F)
  barcodes2=barcodes[which(!is.na(barcodes$number)),]
  
  m1=read.delim(paste0(SCAFE,files3[i]), stringsAsFactors = F, header=F, skip=3, sep=" ")
  m1=left_join(m1,barcodes2[,c(2,3)], by="V2",copy=F)
  m2=m1[which(!is.na(m1$number)),]
  im3=m2%>%group_by(V1,number)%>%dplyr::summarise(V3=sum(V3))
  m3=rbind(m3,im3)}

setwd(paste0(metacell,"tCRE_ori_ATAC_region"))
m4=m3%>%group_by(V1,number)%>%dplyr::summarise(V3=sum(V3))
write.table(m4,paste0(metacell,"tCRE_ori_ATAC_region/matrix.mtx"),col.names=F,row.names=F, sep=" ",quote=F)

m=read.delim(paste0(SCAFE,files3[1]), stringsAsFactors = F, header=F, nrows=3)
m$V1[3]=paste0("463966 ",nrow(latent1), " ", nrow(m4))
write.table(m,paste0(metacell,"tCRE_ori_ATAC_region/matrix.mtx.head"),col.names=F,row.names=F, sep=" ",quote=F)
genes=read.delim(paste0(SCAFE,files2[1]), stringsAsFactors = F, header=F)
write.table(genes,paste0(metacell,"tCRE_ori_ATAC_region/genes.tsv"),col.names=F,row.names=F, sep="\t",quote=F)
write.table(latent1[,1],paste0(metacell,"tCRE_ori_ATAC_region/barcodes.tsv"), col.names=F, row.names=F, sep="\t", quote=F)

#===============================================================================
# ATAC sc matrix preparation -> metacell matrix
All.atac.raw <- readRDS(paste0(bigDATA,"All.atac.raw.3samples.rds"))
data = as.data.frame(summary(All.atac.raw@assays[["peaks"]]@counts))
peaks = data.frame(cbind(rownames(All.atac.raw),rownames(All.atac.raw)))
barcode = data.frame(colnames(All.atac.raw))
barcode$V2=1:nrow(barcode)
colnames(barcode)[1]="V1"
barcode=left_join(barcode, latent0[,c(2,15)], by=c("V1"="X"), copy=F)
barcode2=barcode[which(!is.na(barcode$number)),]

data=left_join(data,barcode2[,c(2,3)], by=c("j"="V2"),copy=F)
data2=data[which(!is.na(data$number)),]
data3=data2%>%group_by(i,number)%>%dplyr::summarise(x=sum(x))
write.table(data3,paste0(metacell,"scATAC/matrix.mtx"),col.names=F,row.names=F, sep=" ",quote=F)

m$V1[3]=paste0(nrow(peaks)," ",nrow(latent1), " ", nrow(data3))
write.table(m,paste0(metacell,"scATAC/matrix.mtx.head"),col.names=F,row.names=F, sep=" ",quote=F)
write.table(peaks,paste0(metacell,"scATAC/genes.tsv"),col.names=F,row.names=F, sep="\t",quote=F)
write.table(latent1[,1],paste0(metacell,"scATAC/barcodes.tsv"), col.names=F, row.names=F, sep="\t", quote=F)


#===============================================================================
# put together ATAC and TSS signal per CRE at metacell level
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
md1=read.delim(paste0(metacell,"output/metadata2.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
latent0=read.delim(paste0(metacell,"output/metadata.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
md1=left_join(md1,unique(latent0[,c(12,1)]),by="SEACell",copy=F)
latent0$sampleID="iPS"
latent0$sampleID[grep("NSC_",latent0$X)]="NSC"
latent0$sampleID[grep("_5",latent0$X)]="NSC"
latent0$sampleID[grep("_6",latent0$X)]="NSC"
latent0$sampleID[grep("NRN_",latent0$X)]="Neuron"
latent0$sampleID[grep("_3",latent0$X)]="Neuron"
latent0$sampleID[grep("_4",latent0$X)]="Neuron"
latent1=latent0%>%group_by(SEACell,sampleID)%>%dplyr::summarise(count=n())
latent1a=latent1%>%group_by(SEACell)%>%dplyr::slice_max(count)
colnames(latent1a)[2]="major_sampleID"

latent2a=spread(latent1,key=2, value=3)
latent2a[is.na(latent2a)]=0
md1=left_join(md1,latent2a,by="SEACell",copy=F)
md1=left_join(md1,latent1a[,c(1,2)],by="SEACell",copy=F)

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
# select 393 metacells
combined=subset(x = combined, subset = ATAC_cell >= 3 & RNA_cell >=3)
saveRDS(combined, paste0(bigDATA,"metacell.combined.feature.Rds"))

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

saveRDS(combined, paste0(bigDATA,"metacell.combined.feature.Rds"))

#=================================================================
#define which CREs are transcribed, read > 0 in any sample
combined=readRDS(paste0(bigDATA,"metacell.combined.feature.Rds"))
data1=data.frame(rowSums(data.frame(combined@assays[["tCRE"]]@counts)))
data1$peakID=rownames(data1)
colnames(data1)[1]="tCRE_count"
write.table(data1[,c(2,1)],gzfile(paste0(metacell,"output/tCRE_count_in_all_CRE.tsv.gz")), col.names=T, row.names=F, sep="\t",quote=F)

#===============================================================================
#take the single cell pseudotime and take the mean for metacell
pseudotime=read.delim(paste0(primary_folder,"Data_and_code/sc5nRNA/output/pseudotime_table_scRNA.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
pseudotime=left_join(pseudotime[,c(1,37:41)],latent0[,c(1,2,12)],by=c("cellID" = "X"),copy=F)
pseudotime1=pseudotime%>%group_by(SEACell)%>%dplyr::summarise(pseudotime=mean(pseudotime))
pseudotime1=pseudotime1[order(pseudotime1$pseudotime),]
pseudotime1$T1=1:nrow(pseudotime1)
write.table(pseudotime1, gzfile(paste0(metacell,"output/pseudotime.scRNA.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

#=====================================================================================================
# take the scRNA UMAP and combine to metacell
cluster1=c("iPS_0_undifferentiaed","iPS_1_neuroepithelial","NSC_0_iPS_like", "NSC_1_differentiating", "NSC_2_astrocyte_like", "neuron_progenitor_n_schwann", "neuron_immature", "neuron_mature","neuron_mature_2", "OPC_like")
All.RNA <- readRDS(paste0(bigDATA,"20220331_All.rna.flt.rds"))
umap1=data.frame(All.RNA@reductions[["umap_harmony"]]@cell.embeddings)
umap1$cellID=rownames(umap1)
umap1=left_join(umap1,latent0[,c(1,2,12)],by=c("cellID"="X"),copy=F)
umap2=umap1%>%group_by(SEACell,cluster)%>%dplyr::summarise(umapharmony_1=mean(umapharmony_1),umapharmony_2=mean(umapharmony_2))
umap2$group="meta cell"
umap1$group="single cell"
colnames(umap2)[1]="cellID"
umap3=rbind(umap1[,c(3,4,6,1,2)],umap2[,c(1,2,5,3,4)])
write.table(umap3,gzfile(paste0(metacell,"output/harmony.scRNA.umap.with.metacell.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

write.table(umap3,gzfile(paste0(path_fig1_data,"harmony.scRNA.umap.with.metacell.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
# ->for fig 1b

# take the snATAC UMAP and combine to metacell
All.atac <- readRDS("/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/All.atac.rds")
Aumap1=data.frame(All.atac@reductions[["umap_integrated"]]@cell.embeddings)
Aumap1$cellID=rownames(Aumap1)
Aumap1=left_join(Aumap1,latent0[,c(1,2,12)],by=c("cellID"="X"),copy=F)
Aumap2=Aumap1%>%group_by(SEACell,cluster)%>%dplyr::summarise(UMAP_1=mean(UMAP_1),UMAP_2=mean(UMAP_2))
Aumap2$group="meta cell"
Aumap1$group="single cell"
colnames(Aumap2)[1]="cellID"
Aumap3=rbind(Aumap1[,c(3,4,6,1,2)],Aumap2[,c(1,2,5,3,4)])
Aumap3$cluster=factor(Aumap3$cluster, levels=cluster1)
write.table(Aumap3,gzfile(paste0(metacell,"output/harmony.snATAC.umap.with.metacell.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

write.table(Aumap3,gzfile(paste0(path_fig1_data,"harmony.snATAC.umap.with.metacell.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
# ->for fig 1b


#===============================================================================
# expression level for each cell cluster
combined=readRDS(paste0(bigDATA,"metacell.combined.feature.Rds"))
avg_exp = AverageExpression(combined, group.by = "cluster",  assays = "RNA", layer = "data")
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

ATACm = as.data.frame(AverageExpression(object = combined, group.by = "cluster", assays = "aCRE", layer = "data")$aCRE)
ATACm$peakID=rownames(ATACm)
TSSm = as.data.frame(AverageExpression(object = combined, group.by = "cluster", assays = "tCRE", layer = "data")$tCRE)
TSSm$peakID=rownames(TSSm)
write.table(ATACm, gzfile(paste0(metacell,"output/CRE_ATAC_per_cluster.matrix.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(TSSm, gzfile(paste0(metacell,"output/CRE_TSS_per_cluster.matrix.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

tCRE_counts <- as.data.frame(AggregateExpression(object = combined, assays = "tCRE", group.by = "major_sampleID")$tCRE)
tCRE_counts$peakID=rownames(tCRE_counts)
write.table(tCRE_counts, gzfile(paste0(metacell,"output/CRE_TSS_per_sample.matrix.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================



saveRDS(combined, paste0(bigDATA,"metacell.combined.feature.Rds"))



