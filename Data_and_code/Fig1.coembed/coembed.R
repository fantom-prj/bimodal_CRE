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
RNA_folder=paste0(primary_folder,"Data_and_code/sc5nRNA/")
ATAC_folder=paste0(primary_folder,"Data_and_code/scATAC/")
scDART_SEACell=paste0(primary_folder,"Data_and_code/Fig1.coembed/scDART_n_SEACell/")
path_fig1_data=paste0(primary_folder,"Fig1/data/")
resource_folder=paste0(primary_folder,"Data_and_code/resources/")
peak_folder=paste0(primary_folder,"Data_and_code/scATAC/merged_peak/")

#===============================================================================
All.rna <- readRDS(paste0(RNA_folder,"20220331_All.rna.flt.rds"))
All.atac <- readRDS(paste0(ATAC_folder,"All.atac.rds"))
All.atac[["RNA"]] <- NULL

Idents(All.atac) <- "celltype_manual_LT"
Idents(All.rna) <- "celltype_manual"

#take gene region for bedtools
need2=rownames(All.rna)
gencode=read.delim(paste0(resource_folder,"cellranger_genes.gtf.gz"),skip=5,header=F, stringsAsFactors = F, check.names=F)
gencode=gencode[which(gencode$V3 == "gene"),]
gencode$gene_name=sapply(strsplit(gencode$V9,";"),"[",4)
gencode$gene_name=gsub(" gene_name ","",gencode$gene_name)
gencode1=gencode[which(gencode$gene_name %in% need2),]

setdiff(need2,gencode1$gene_name)
#"TBCE"      "LINC01238" "CYB561D2"  "MATR3"     "HSPA14"    "GGT1"      "TMSB15B"   gene name are redundant, introduce "TBCE.1"      "LINC01238.1" "CYB561D2.1"  "MATR3.1"     "HSPA14.1"    "GGT1.1"      "TMSB15B.1"  , removed.
gencode2=gencode1[,c(1,4,5,10,6,7)]
gencode2=gencode2[-which(gencode2$gene_name %in% c("TBCE","LINC01238","CYB561D2","MATR3","HSPA14","GGT1","TMSB15B")),]
gencode2$V4=gencode2$V4-1
gencode2$V4[which(gencode2$V7 == "+")]=gencode2$V4[which(gencode2$V7 == "+")]-2000
gencode2$V5[which(gencode2$V7 == "-")]=gencode2$V5[which(gencode2$V7 == "-")]+2000
options(scipen=999)
write.table(gencode2[order(gencode2$V1,gencode2$V4),],gzfile(paste0(resource_folder,"cellranger_genes19524.extended_2000up.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)

#======================================
#bedtools intersect
setwd(resource_folder)
system("bedtools intersect -wa -wb -a ",peak_folder,"filtered_peaks.merged.all.bed.gz -b cellranger_genes19524.extended_2000up.bed.gz | gzip > atac.peak.genes19524.extended_2000up.bed.gz")

#======================================
# perform co-embed by Seurat & generate input files for scDART (for each cell cluster)
peak=read.delim(paste0(resource_folder,"atac.peak.genes19524.extended_2000up.bed.gz"), header=F, stringsAsFactors=F, check.names=F)
peak$count=1

identity=unique(Idents(All.rna))
###
for (i in 1: length(identity)){
  path1=paste0(scDART_SEACell,identity[i])
  setwd(path1)
  
  iPS_1_neuroepithelial.All.atac = subset(x = All.atac, idents = identity[i])
  iPS_1_neuroepithelial.All.rna = subset(x = All.rna, idents = identity[i])
  iPS_1_neuroepithelial.All.rna <- FindVariableFeatures(iPS_1_neuroepithelial.All.rna)
  
  genes.use <- VariableFeatures(iPS_1_neuroepithelial.All.rna)
  
  transfer.anchors <- FindTransferAnchors(reference = iPS_1_neuroepithelial.All.rna, 
                                          reference.assay = "RNA",
                                          normalization.method = "LogNormalize",
                                          features = VariableFeatures(object = iPS_1_neuroepithelial.All.rna), 
                                          query = iPS_1_neuroepithelial.All.atac, 
                                          query.assay = "GeneActivity", 
                                          reduction = "cca")
  
  refdata <- GetAssayData(iPS_1_neuroepithelial.All.rna, assay = "RNA", slot = "data")[genes.use, ]
  imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = iPS_1_neuroepithelial.All.atac[["lsi"]], dims = 2:30)
  iPS_1_neuroepithelial.All.atac[["RNA"]] <- imputation
  
  tmp.x=iPS_1_neuroepithelial.All.rna
  DefaultAssay(tmp.x) <- "RNA"
  tmp.x[["SCT"]] <- NULL
  tmp.x$dataType <- "scRNA"
  tmp.y=iPS_1_neuroepithelial.All.atac
  DefaultAssay(tmp.y) <- "RNA"
  tmp.y[["GeneActivity"]] <- NULL
  tmp.y[["peaks"]] <- NULL
  tmp.y$dataType <- "scATAC"
  coembed <- merge(x = tmp.x, y = tmp.y)
  coembed <- FindVariableFeatures(coembed)
  coembed <- ScaleData(coembed, do.scale = FALSE)
  coembed <- RunPCA(coembed)
  coembed <- RunUMAP(coembed, dims = 1:30)
  a1=DimPlot(coembed, reduction = "umap", group.by = 'dataType') + ggtitle(identity[i])
  jpeg("coembed.jpg", width = 5, height = 4, units = "in", res=600)
  print(a1)
  dev.off()
  n2000=data.frame(coembed@meta.data)
  n2000.umap=data.frame(coembed@reductions[["umap"]]@cell.embeddings)
  n2000.pca=data.frame(coembed@reductions[["pca"]]@cell.embeddings)
  n2000$seq="ATAC"
  n2000$seq[which(is.na(n2000$TSS_fragments))]="RNA"
  n2000$cluster=identity[i]
  n2000=cbind(n2000[,c(69,70)],n2000.umap)
  n2000.pca=cbind(n2000[,c(1,2)],n2000.pca)
  write.table(n2000,"seurat.coembed.umap.tsv",col.names=T, row.names=T, sep="\t", quote=F)
  write.table(n2000.pca,"seurat.coembed.pca.tsv",col.names=T, row.names=T, sep="\t", quote=F)
  
  peak1=unique(peak[which(peak$V10 %in% genes.use),c(4,10,13)])
  peak1=spread(peak1, key=2, value=3)
  peak1[is.na(peak1)]=0
  row.names(peak1)=peak1$V4
  peak1=peak1[,-1]
  write.csv(peak1, "region2gene.csv", quote=F, row.names=T, col.names=T)
  
  need=which(rownames(iPS_1_neuroepithelial.All.atac[["peaks"]]) %in% unique(row.names(peak1)))
  iPS_1_neuroepithelial.All.atac[["peaks"]] <- subset(iPS_1_neuroepithelial.All.atac[["peaks"]], features = rownames(iPS_1_neuroepithelial.All.atac[["peaks"]])[need])
  ATAC_count=data.frame(t(iPS_1_neuroepithelial.All.atac@assays[["peaks"]]@counts), check.names=F)
  write.csv(ATAC_count, "ATAC.count.csv", quote=F, row.names=T, col.names=T)
  
  ATAC_label=iPS_1_neuroepithelial.All.atac@meta.data
  write.table(ATAC_label[,36], "ATAC.label.txt", quote=F, row.names=F, col.names=F)
  
  need=which(rownames(iPS_1_neuroepithelial.All.rna[["RNA"]]) %in% colnames(peak1))
  iPS_1_neuroepithelial.All.rna[["RNA"]] <- subset(iPS_1_neuroepithelial.All.rna[["RNA"]], features = rownames(iPS_1_neuroepithelial.All.rna[["RNA"]])[need])
  df= data.frame(t(iPS_1_neuroepithelial.All.rna[["RNA"]]@counts), check.names = FALSE)
  write.csv(df, "RNA.count.csv", quote=F, row.names=T, col.names=T)
  
  RNAlabel=iPS_1_neuroepithelial.All.rna@meta.data
  write.table(RNAlabel[,34], "RNA.label.txt", quote=F, row.names=F, col.names=F)
  
  rm(iPS_1_neuroepithelial.All.atac)
  rm(iPS_1_neuroepithelial.All.rna)
  gc()}

#===============================================================================
#run scDART per cell cluster
#run [cluster_name]_scdart.py in "Fig1.coembed"

#===============================================================================
# compare co-embed between Seurat and scDART

i=1
path1=paste0(scDART_SEACell,identity[i])
setwd(path1)
umap=read.delim("seurat.coembed.umap.tsv",header=T, stringsAsFactors = F, check.names = F)
umap$cellID=rownames(umap)
pca=read.delim("seurat.coembed.pca.tsv",header=T, stringsAsFactors = F, check.names = F)
umap=cbind(umap,pca[,c(3,4)])
atac.pca=read.csv("ATAC.PCA.csv",header=T, stringsAsFactors = F, check.names = F)
rna.pca=read.csv("RNA.PCA.csv",header=T, stringsAsFactors = F, check.names = F)
atac.label=fread("ATAC.count.csv",header=F, select=1, skip=1)
rna.label=fread("RNA.count.csv",header=F, select=1, skip=1)
atac.pca=cbind(atac.label,atac.pca[,c(2,3)])
rna.pca=cbind(rna.label,rna.pca[,c(2,3)])
scDART=rbind(atac.pca,rna.pca)
colnames(scDART)=c("cellID","scDART_PC1","scDART_PC2")
umap=left_join(umap, scDART, by="cellID",copy=F)

for (i in 2: length(identity)){
  path1=paste0(scDART_SEACell,identity[i])
  setwd(path1)
  umap1=read.delim("seurat.coembed.umap.tsv",header=T, stringsAsFactors = F, check.names = F)
  umap1$cellID=rownames(umap1)
  pca1=read.delim("seurat.coembed.pca.tsv",header=T, stringsAsFactors = F, check.names = F)
  umap1=cbind(umap1,pca1[,c(3,4)])
  atac.pca=read.csv("ATAC.PCA.csv",header=T, stringsAsFactors = F, check.names = F)
  rna.pca=read.csv("RNA.PCA.csv",header=T, stringsAsFactors = F, check.names = F)
  atac.label=fread("ATAC.count.csv",header=F, select=1, skip=1)
  rna.label=fread("RNA.count.csv",header=F, select=1, skip=1)
  atac.pca=cbind(atac.label,atac.pca[,c(2,3)])
  rna.pca=cbind(rna.label,rna.pca[,c(2,3)])
  scDART=rbind(atac.pca,rna.pca)
  colnames(scDART)=c("cellID","scDART_PC1","scDART_PC2")
  umap1=left_join(umap1, scDART, by="cellID",copy=F)
  umap=rbind(umap,umap1)
}

write.table(umap,gzfile(paste0(scDART_SEACell,"umap_for_all_clusters.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)
mumap=umap%>%group_by(cluster,seq)%>%dplyr::summarise(count=n())

a1=ggplot() + 
  labs(x="UMAP_1", y ="UMAP_2", title="Co-embed by Seurat", color="single cell")+
  geom_point(data=umap[sample(nrow(umap),nrow(umap)),], mapping=aes(x=UMAP_1, color=as.factor(seq), y=UMAP_2), size=0.25, alpha=0.3)+
  scale_color_manual(values=c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF","black", "gold"))+
  facet_wrap(facets=vars(cluster), ncol=4,nrow=3, scale="fixed")+
  theme(  panel.background = element_rect(fill = "white", color = "white", size = 0.5, linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = "grey90"), 
          axis.text.x = element_text(color ="black", angle=45, hjust=1, vjust=1),
          axis.text.y = element_text(color ="black"),
          plot.title = element_text(hjust = 0.5,  color ="black", size=11),
          plot.subtitle = element_text(color="black", hjust=0.5),
          text = element_text(size=10, family="Arial"),
          strip.text = element_text(size=6, family="Arial"),
          legend.background = element_rect(fill="white", color="black"),
          legend.title = element_text(hjust=0.5),
          legend.text = element_text(lineheight = 0.8),
          axis.line = element_line(size = 0.5, colour = "black"),
          axis.ticks = element_line(size = 0.5,colour = "black"),
          strip.background =element_rect(fill="grey90"),
          legend.key = element_blank(),
          legend.position = c(0.88, 0.12),
          plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))
jpeg(paste0(scDART_SEACell,"coembed.seurat.jpg"), width = 4, height = 4, units = "in", res=600)
print (a1)
dev.off() 

b1=ggplot() + 
  labs(x="scDART_PC1", y ="scDART_PC2", title="Co-embed by scDART", color="single cell")+
  geom_point(data=umap[sample(nrow(umap),nrow(umap)),], mapping=aes(x=scDART_PC1, color=as.factor(seq), y=scDART_PC2), size=0.25, alpha=0.3)+
  scale_color_manual(values=c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF","black", "gold"))+
  facet_wrap(facets=vars(cluster), ncol=4,nrow=3, scale="fixed")+
  theme(  panel.background = element_rect(fill = "white", color = "white", size = 0.5, linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = "grey90"), 
          axis.text.x = element_text(color ="black", angle=45, hjust=1, vjust=1),
          axis.text.y = element_text(color ="black"),
          plot.title = element_text(hjust = 0.5,  color ="black", size=11),
          plot.subtitle = element_text(color="black", hjust=0.5),
          text = element_text(size=10, family="Arial"),
          strip.text = element_text(size=6, family="Arial"),
          legend.background = element_rect(fill="white", color="black"),
          legend.title = element_text(hjust=0.5),
          legend.text = element_text(lineheight = 0.8),
          axis.line = element_line(size = 0.5, colour = "black"),
          axis.ticks = element_line(size = 0.5,colour = "black"),
          strip.background =element_rect(fill="grey90"),
          legend.key = element_blank(),
          legend.position = c(0.88, 0.12),
          plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))
jpeg(paste0(scDART_SEACell,"coembed.scDART.jpg"), width = 4, height = 4, units = "in", res=600)
print (b1)
dev.off() 

#===============================================================================
#Put back the PCA into each co-embeded cluster

library(Seurat)
library(SeuratData)
library(SeuratDisk)

All.rna <- readRDS(paste0(RNA_folder,"20220331_All.rna.flt.rds"))
All.atac <- readRDS(paste0(ATAC_folder,"All.atac.rds"))

DefaultAssay(All.rna) <- "RNA"
DefaultAssay(All.atac) <- "RNA"
All.rna$dataType <- "RNA"
All.atac$dataType <- "ATAC"
All.rna[["SCT"]] <- NULL
All.atac[["GeneActivity"]] <- NULL
All.atac[["peaks"]] <- NULL

Idents(All.atac) <- "celltype_manual_LT"
Idents(All.rna) <- "celltype_manual"

#umap=read.delim("/analysisdata/fantom6/Interactome/single_cell_wallace/scDART/umap_for_all_clusters.tsv", header=T, stringsAsFactor=F, check.names=F)

identity=unique(Idents(All.rna))
###
for (i in 1: length(identity)){
  path1=paste0(scDART_SEACell,identity[i])
  setwd(path1)
  
  iPS_1_neuroepithelial.All.atac = subset(x = All.atac, idents = identity[i])
  iPS_1_neuroepithelial.All.rna = subset(x = All.rna, idents = identity[i])
  
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
  
  coembed <- merge(x = iPS_1_neuroepithelial.All.atac, y = iPS_1_neuroepithelial.All.rna)
  umap=PCA[colnames(coembed),c(1,2)]
  pca=PCA[colnames(coembed),]
  
  umap_coordinates_mat <- as(umap, "matrix")
  pca_coordinates_mat <- as(pca, "matrix")
  coembed[['pca']] <- CreateDimReducObject(embeddings = pca_coordinates_mat, key = "pca_", global = T, assay = "RNA")
  coembed[['umap']] <- CreateDimReducObject(embeddings = umap_coordinates_mat, key = "umap_", global = T, assay = "RNA")
  
  DimPlot(coembed, reduction = "pca", group.by = 'dataType') + ggtitle(identity[i])
  DimPlot(coembed, reduction = "umap", group.by = 'dataType') + ggtitle(identity[i])
  
  saveRDS(coembed, file = "coembed.rds")
  SaveH5Seurat(coembed, filename = "coembed.h5Seurat")
  Convert("coembed.h5Seurat", dest = "h5ad")}


#===================================================================================================
# schwann cluster is too small, add a cluster that combine progenitor and schwann

peak=read.delim("/analysisdata/fantom6/Interactome/single_cell_wallace/cellranger_used_gene/atac.peak.genes19524.extended_2000up.bed", header=F, stringsAsFactors=F, check.names=F)
peak$count=1

All.rna <- readRDS(paste0(RNA_folder,"20220331_All.rna.flt.rds"))
All.atac <- readRDS(paste0(ATAC_folder,"All.atac.rds"))
All.atac[["RNA"]] <- NULL

Idents(All.atac) <- "celltype_manual_LT"
Idents(All.rna) <- "celltype_manual"

identity="neuron_progenitor_n_schwann"
###
path1=paste0(scDART_SEACell,identity[1])
setwd(path1)

iPS_1_neuroepithelial.All.atac = subset(x = All.atac, idents = c("neuron_progenitor","neuron_schwann_like"))
iPS_1_neuroepithelial.All.rna = subset(x = All.rna, idents = c("neuron_progenitor","neuron_schwann_like"))
iPS_1_neuroepithelial.All.rna <- FindVariableFeatures(iPS_1_neuroepithelial.All.rna)

genes.use <- VariableFeatures(iPS_1_neuroepithelial.All.rna)

transfer.anchors <- FindTransferAnchors(reference = iPS_1_neuroepithelial.All.rna, 
                                        reference.assay = "RNA",
                                        normalization.method = "LogNormalize",
                                        features = VariableFeatures(object = iPS_1_neuroepithelial.All.rna), 
                                        query = iPS_1_neuroepithelial.All.atac, 
                                        query.assay = "GeneActivity", 
                                        reduction = "cca")

refdata <- GetAssayData(iPS_1_neuroepithelial.All.rna, assay = "RNA", slot = "data")[genes.use, ]
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = iPS_1_neuroepithelial.All.atac[["lsi"]], dims = 2:30)
iPS_1_neuroepithelial.All.atac[["RNA"]] <- imputation

tmp.x=iPS_1_neuroepithelial.All.rna
DefaultAssay(tmp.x) <- "RNA"
tmp.x[["SCT"]] <- NULL
tmp.x$dataType <- "scRNA"
tmp.y=iPS_1_neuroepithelial.All.atac
DefaultAssay(tmp.y) <- "RNA"
tmp.y[["GeneActivity"]] <- NULL
tmp.y[["peaks"]] <- NULL
tmp.y$dataType <- "scATAC"
coembed <- merge(x = tmp.x, y = tmp.y)
coembed <- FindVariableFeatures(coembed)
coembed <- ScaleData(coembed, do.scale = FALSE)
coembed <- RunPCA(coembed)
coembed <- RunUMAP(coembed, dims = 1:30)
a1=DimPlot(coembed, reduction = "umap", group.by = 'dataType') + ggtitle("neuron_progenitor_n_schwann")
jpeg("coembed.jpg", width = 5, height = 4, units = "in", res=600)
print(a1)
dev.off()
n2000=data.frame(coembed@meta.data)
n2000.umap=data.frame(coembed@reductions[["umap"]]@cell.embeddings)
n2000.pca=data.frame(coembed@reductions[["pca"]]@cell.embeddings)
n2000$seq="ATAC"
n2000$seq[which(is.na(n2000$TSS_fragments))]="RNA"
n2000$cluster="neuron_progenitor_n_schwann"
n2000=cbind(n2000[,c(69,70)],n2000.umap)
n2000.pca=cbind(n2000[,c(1,2)],n2000.pca)
write.table(n2000,"seurat.coembed.umap.tsv",col.names=T, row.names=T, sep="\t", quote=F)
write.table(n2000.pca,"seurat.coembed.pca.tsv",col.names=T, row.names=T, sep="\t", quote=F)

peak1=unique(peak[which(peak$V10 %in% genes.use),c(4,10,13)])
peak1=spread(peak1, key=2, value=3)
peak1[is.na(peak1)]=0
row.names(peak1)=peak1$V4
peak1=peak1[,-1]
write.csv(peak1, "region2gene.csv", quote=F, row.names=T, col.names=T)

need=which(rownames(iPS_1_neuroepithelial.All.atac[["peaks"]]) %in% unique(row.names(peak1)))
iPS_1_neuroepithelial.All.atac[["peaks"]] <- subset(iPS_1_neuroepithelial.All.atac[["peaks"]], features = rownames(iPS_1_neuroepithelial.All.atac[["peaks"]])[need])
ATAC_count=data.frame(t(data.frame(iPS_1_neuroepithelial.All.atac@assays[["peaks"]]@counts, check.names = FALSE)), check.names=F)
write.csv(ATAC_count, "ATAC.count.csv", quote=F, row.names=T, col.names=T)

ATAC_label=iPS_1_neuroepithelial.All.atac@meta.data
write.table(ATAC_label[,36], "ATAC.label.txt", quote=F, row.names=F, col.names=F)

need=which(rownames(iPS_1_neuroepithelial.All.rna[["RNA"]]) %in% colnames(peak1))
iPS_1_neuroepithelial.All.rna[["RNA"]] <- subset(iPS_1_neuroepithelial.All.rna[["RNA"]], features = rownames(iPS_1_neuroepithelial.All.rna[["RNA"]])[need])
df= data.frame(t(data.frame(iPS_1_neuroepithelial.All.rna[["RNA"]]@counts, check.names = FALSE)), check.names = FALSE)
write.csv(df, "RNA.count.csv", quote=F, row.names=T, col.names=T)

RNAlabel=iPS_1_neuroepithelial.All.rna@meta.data
write.table(RNAlabel[,34], "RNA.label.txt", quote=F, row.names=F, col.names=F)

rm(iPS_1_neuroepithelial.All.atac)
rm(iPS_1_neuroepithelial.All.rna)
gc()

##########

library(Seurat)
library(SeuratData)
library(SeuratDisk)

All.rna <- readRDS(paste0(RNA_folder,"20220331_All.rna.flt.rds"))
All.atac <- readRDS(paste0(ATAC_folder,"All.atac.rds"))
DefaultAssay(All.rna) <- "RNA"
DefaultAssay(All.atac) <- "RNA"
All.rna$dataType <- "RNA"
All.atac$dataType <- "ATAC"
All.rna[["SCT"]] <- NULL
All.atac[["GeneActivity"]] <- NULL
All.atac[["peaks"]] <- NULL

Idents(All.atac) <- "celltype_manual_LT"
Idents(All.rna) <- "celltype_manual"

identity="neuron_progenitor_n_schwann"
###
path1=paste0(scDART_SEACell,identity[1])
setwd(path1)

iPS_1_neuroepithelial.All.atac = subset(x = All.atac, idents = c("neuron_progenitor","neuron_schwann_like"))
iPS_1_neuroepithelial.All.rna = subset(x = All.rna, idents = c("neuron_progenitor","neuron_schwann_like"))

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

coembed <- merge(x = iPS_1_neuroepithelial.All.atac, y = iPS_1_neuroepithelial.All.rna)
umap=PCA[colnames(coembed),c(1,2)]
pca=PCA[colnames(coembed),]

umap_coordinates_mat <- as(umap, "matrix")
pca_coordinates_mat <- as(pca, "matrix")
coembed[['pca']] <- CreateDimReducObject(embeddings = pca_coordinates_mat, key = "pca_", global = T, assay = "RNA")
coembed[['umap']] <- CreateDimReducObject(embeddings = umap_coordinates_mat, key = "umap_", global = T, assay = "RNA")

saveRDS(coembed, file = "coembed.rds")
SaveH5Seurat(coembed, filename = "coembed.h5Seurat")
Convert("coembed.h5Seurat", dest = "h5ad")

#===============================================================================
# plot heatmap for all clusters, including the combined one

#scDART PCA cluster
library(pheatmap)
library("RColorBrewer")

files=list.files(pattern="ATAC.latent.csv" , path=scDART_SEACell ,recursive =T)
identity = sapply(strsplit(files,"/"),"[",1)
for (i in 1: length(identity)){
  path1=paste0(scDART_SEACell,identity[i])
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
  latent.label=data.frame(c(rep("ATAC",nrow(ATAC.latent)),rep("RNA",nrow(RNA.latent))))
  colnames(latent.label)="Tech"
  rownames(latent.label)=rownames(latent)
  ann_colors = list(Tech = c(ATAC="#00A087FF","RNA"="#7E6148FF"))
  
  out1=pheatmap(latent, fontsize_col = 8, fontsize_row = 3.5, angle_col=45 , 
                color = colorRampPalette(c("navy","white", "firebrick3"))(25), 
                main = paste0("Latent dimension of ",identity[i]), 
                #clustering_method = "ward.D2",
                annotation_row = latent.label,
                annotation_colors = ann_colors,
                #scale = "row",
                #treeheight_row = 25,
                cluster_cols = F,
                show_rownames = F)
  out2=pheatmap(PCA, fontsize_col = 8, fontsize_row = 3.5, angle_col=45 , 
                color = colorRampPalette(c("navy","white", "firebrick3"))(25), 
                main = paste0("PCA of ",identity[i]), 
                #clustering_method = "ward.D2",
                annotation_row = latent.label,
                annotation_colors = ann_colors,
                #scale = "row",
                #treeheight_row = 25,
                cluster_cols = F,
                show_rownames = F)
  jpeg(file="scDART.latent.heatmap.jpeg", width = 3.5, height = 4, units = "in", res=600)
  print(out1)
  dev.off()
  jpeg(file="scDART.PCA.heatmap.jpeg", width = 3.5, height = 4, units = "in", res=600)
  print(out2)
  dev.off()}

#===============================================================================
#run SEACell
#run [cluster_name]_SEACell.py in "Fig1.coembed"

#Parse SEACell result
files=list.files(pattern="ATAC.latent.csv" , path=scDART_SEACell, recursive =T)
identity = sapply(strsplit(files,"/"),"[",1)
identity=identity[-8]

path1=paste0(scDART_SEACell,identity[1])
setwd(path1)

seacell=read.delim("metadata.tsv", header=F, skip=1, stringsAsFactors =F)
seacell$metacell=paste0(identity[i],"_meta",1:nrow(seacell))
seacell=data.frame(separate_rows(seacell, V1, sep="\t"))

ATAC.PCA=read.csv("ATAC.PCA.csv")
RNA.PCA=read.csv("RNA.PCA.csv")
ATAC.label=fread("ATAC.count.csv",header=F, select=1, skip=1)
RNA.label=fread("RNA.count.csv",header=F, select=1, skip=1)

ATAC.PCA=data.frame(cbind(ATAC.label,ATAC.PCA[,c(2:9)]))
RNA.PCA=data.frame(cbind(RNA.label,RNA.PCA[,c(2:9)]))
ATAC.PCA$dataType="ATAC"
RNA.PCA$dataType="RNA"
PCA=rbind(ATAC.PCA, RNA.PCA)
PCA=left_join(PCA, seacell,by="V1", copy=F)
rownames(PCA)=PCA$V1
colnames(PCA)[11]="SEACell"
PCA$cluster=identity[1]

for (i in 2: length(identity)){
  path1=paste0(scDART_SEACell,identity[i])
  setwd(path1)
  
  seacell=read.delim("metadata.tsv", header=F, skip=1, stringsAsFactors =F)
  seacell$metacell=paste0(identity[i],"_meta",1:nrow(seacell))
  seacell=data.frame(separate_rows(seacell, V1, sep="\t"))
  
  ATAC.PCA=read.csv("ATAC.PCA.csv")
  RNA.PCA=read.csv("RNA.PCA.csv")
  ATAC.label=fread("ATAC.count.csv",header=F, select=1, skip=1)
  RNA.label=fread("RNA.count.csv",header=F, select=1, skip=1)
  
  ATAC.PCA=data.frame(cbind(ATAC.label,ATAC.PCA[,c(2:9)]))
  RNA.PCA=data.frame(cbind(RNA.label,RNA.PCA[,c(2:9)]))
  ATAC.PCA$dataType="ATAC"
  RNA.PCA$dataType="RNA"
  PCA1=rbind(ATAC.PCA, RNA.PCA)
  PCA1=left_join(PCA1, seacell,by="V1", copy=F)
  rownames(PCA1)=PCA1$V1
  colnames(PCA1)[11]="SEACell"
  PCA1$cluster=identity[i]
  PCA=rbind(PCA, PCA1)}
write.table(PCA,gzfile(paste0(scDART_SEACell,"scDART_SEACell_for_all_clusters.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)
write.table(PCA,gzfile(paste0(path_fig1_data,"scDART_SEACell_for_all_clusters.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

# -> for fig ex1b-c

#===============================================================================





