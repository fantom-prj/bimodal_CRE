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
library(Seurat)
set.seed(2022)
library(Signac)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(GenomicRanges)
library(future)
library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
set.seed(2017)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TFBSTools)

###color#####
library("ggsci")
mypal = pal_npg("nrc", alpha = 1)(9)
mypal
library("scales")
show_col(mypal)
###############

#primary_folder=[primary_folder]
primary_folder="/analysisdata/fantom6/Interactome/single_cell_wallace/Figure/"
chromvar_folder=paste0(primary_folder,"Data_and_code/Fig3.ChromVar/")
metacell=paste0(primary_folder,"Data_and_code/Fig1.metacell/")
bigDATA=paste0(primary_folder,"Data_and_code/DATA/")
path_fig3_data=paste0(primary_folder,"Fig3/data/")

#===============================================================================
#modify format of jaspar table
jas=read.delim(paste0(chromvar_folder,"jaspar2024/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"), header=F)
jas1=jas[grep(">",jas$V1),]
jas1$V1=gsub(">","",jas1$V1)

write.table(jas[,1], paste0(chromvar_folder,"jaspar2024/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.adjust.txt"), col.names=F, row.names=F, sep="\t", quote=F)
write.table(jas1, paste0(chromvar_folder,"jaspar2024/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.ID.tsv"), col.names=F, row.names=F, sep="\t", quote=F)

#===============================================================================
#chromVar for different subset of CREs
scACREv3=read.delim(paste0(primary_folder,"Data_and_code/scATAC/merged_peak/all_aCRE.PE.final.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
scACREv3=scACREv3[which(scACREv3$tCRE == "yes"),]

pfm0 <- readJASPARMatrix(paste0(chromvar_folder,"jaspar2024/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.adjust.txt"), matrixClass = "PFM")
combined=readRDS(paste0(bigDATA,"metacell.combined.feature.Rds"))
motif.df=read.delim(paste0(chromvar_folder,"jaspar2024/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.ID.tsv"), header=F, stringsAsFactors = F)

#===
# TSS of all tCRE, non transcribed CRE removed to have better GC background
# filter tCRE to those with signal
need.rna=which(rownames(combined[["tCRE"]]) %in% scACREv3$peakID)
combined[["tCRE"]] <- subset(combined[["tCRE"]], features = rownames(combined[["tCRE"]])[need.rna])
#remove peak that is zero across all cells
non_zero_frag_peaks <- which(rowSums(combined[["tCRE"]]) != 0)
combined[["tCRE"]] <- subset(combined[["tCRE"]], features = rownames(combined[["tCRE"]])[non_zero_frag_peaks])

peak_data <- combined[["tCRE"]]@counts
peak_ranges <- combined[["tCRE"]]@ranges
rse <- SummarizedExperiment(assays = list(peak_data), rowRanges = peak_ranges)
rse <- addGCBias(rse, genome = BSgenome.Hsapiens.UCSC.hg38)
names(rse@assays) <- c("counts")

head(rse)
head(rowData(rse))
head(assay(rse))
rowSums(assay(rse))

motif_ix <- matchMotifs(pfm0, rse, genome = BSgenome.Hsapiens.UCSC.hg38, p.cutoff = 1e-04)
bg <- getBackgroundPeaks(object = rse)
dev <- computeDeviations(object = rse, annotations = motif_ix, background_peaks = bg)
variability <- computeVariability(dev)
variability=left_join(variability, motif.df, by=c("name"="V1"))
plotVariability(variability, labels = variability$V2, use_plotly = FALSE) 

path22=paste0(chromvar_folder,"output/")
dev_df=data.frame(t(dev@assays@data@listData[["deviations"]]))
z_df=data.frame(t(dev@assays@data@listData[["z"]]))
write.table(dev_df,gzfile(paste0(path22,"all_tCRE_dev.tsv,gz")), col.names=T, row.names=T, sep="\t", quote=F)
write.table(z_df,gzfile(paste0(path22,"all_tCRE_z.tsv.gz")), col.names=T, row.names=T, sep="\t", quote=F)
write.table(variability,gzfile(paste0(path22,"all_tCRE_var.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#=====
# ATAC of all aCRE
peak_data <- combined[["aCRE"]]@counts
peak_ranges <- combined[["aCRE"]]@ranges
rse1 <- SummarizedExperiment(assays = list(peak_data),rowRanges = peak_ranges)
rse1 <- addGCBias(rse1, genome = BSgenome.Hsapiens.UCSC.hg38)
names(rse1@assays) <- c("counts")
motif_ix1 <- matchMotifs(pfm0, rse1, genome = BSgenome.Hsapiens.UCSC.hg38, p.cutoff = 1e-04)
bg1 <- getBackgroundPeaks(object = rse1)
dev1 <- computeDeviations(object = rse1, annotations = motif_ix1, background_peaks = bg)
variability1 <- computeVariability(dev1)
plotVariability(variability1, use_plotly = FALSE) 
variability1=left_join(variability1, motif.df, by=c("name"="V1"))
dev_df1=data.frame(t(dev1@assays@data@listData[["deviations"]]))
z_df1=data.frame(t(dev1@assays@data@listData[["z"]]))
write.table(dev_df1,gzfile(paste0(path22,"all_aCRE_dev.tsv.gz")), col.names=T, row.names=T, sep="\t", quote=F)
write.table(z_df1,gzfile(paste0(path22,"all_aCRE_z.tsv.gz")), col.names=T, row.names=T, sep="\t", quote=F)
write.table(variability1,gzfile(paste0(path22,"all_aCRE_var.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
# limit to pen, removing unclass and CTCF-CREs (pooled enhancer promoter)
scACREv3=read.delim(paste0(primary_folder,"Data_and_code/scATAC/merged_peak/all_aCRE.PE.final.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
scACREv3=scACREv3[which(scACREv3$tCRE == "yes"),]
scACREv3a=scACREv3[-which(scACREv3$promoter_type_CT == "promoter-like" & is.na(scACREv3$gene)),]
scACREv3a$peakID=gsub("_","-",scACREv3a$peakID)

scACREv3a=scACREv3a[which(scACREv3a$promoter_type_CT != "unclassed"),]
scACREv3a=scACREv3a[which(scACREv3a$promoter_type_CT != "CTCF-alone"),]
scACREv3a%>%group_by(promoter_type_CT)%>%dplyr::summarise(count=n())

combined.pen=combined 
need.atac=which(rownames(combined.pen[["aCRE"]]) %in% scACREv3a$peakID)
need.rna=which(rownames(combined.pen[["tCRE"]]) %in% scACREv3a$peakID)
combined.pen[["aCRE"]] <- subset(combined.pen[["aCRE"]], features = rownames(combined.pen[["aCRE"]])[need.atac])
combined.pen[["tCRE"]] <- subset(combined.pen[["tCRE"]], features = rownames(combined.pen[["tCRE"]])[need.rna])

#===
# TSS
#remove peak that is zero across all cells
non_zero_frag_peaks <- which(rowSums(combined.pen[["tCRE"]]) != 0)
combined.pen[["tCRE"]] <- subset(combined.pen[["tCRE"]], features = rownames(combined.pen[["tCRE"]])[non_zero_frag_peaks])

peak_data <- combined.pen[["tCRE"]]@counts
peak_ranges <- combined.pen[["tCRE"]]@ranges
rse <- SummarizedExperiment(assays = list(peak_data),rowRanges = peak_ranges)
rse <- addGCBias(rse, genome = BSgenome.Hsapiens.UCSC.hg38)
names(rse@assays) <- c("counts")
motif_ix <- matchMotifs(pfm0, rse, genome = BSgenome.Hsapiens.UCSC.hg38, p.cutoff = 1e-04)
bg <- getBackgroundPeaks(object = rse)
dev <- computeDeviations(object = rse, annotations = motif_ix, background_peaks = bg)
variability <- computeVariability(dev)
plotVariability(variability, use_plotly = FALSE) 
variability=left_join(variability, motif.df, by=c("name"="V1"))
dev_df=data.frame(t(dev@assays@data@listData[["deviations"]]))
z_df=data.frame(t(dev@assays@data@listData[["z"]]))
write.table(dev_df,gzfile(paste0(path22,"pen_tCRE_dev.tsv.gz")), col.names=T, row.names=T, sep="\t", quote=F)
write.table(z_df,gzfile(paste0(path22,"pen_tCRE_z.tsv.gz")), col.names=T, row.names=T, sep="\t", quote=F)
write.table(variability,gzfile(paste0(path22,"pen_tCRE_var.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#====
# ATAC
peak_data <- combined.pen[["aCRE"]]@counts
peak_ranges <- combined.pen[["aCRE"]]@ranges
rse1 <- SummarizedExperiment(assays = list(peak_data),rowRanges = peak_ranges)
rse1 <- addGCBias(rse1, genome = BSgenome.Hsapiens.UCSC.hg38)
names(rse1@assays) <- c("counts")
motif_ix1 <- matchMotifs(pfm0, rse1, genome = BSgenome.Hsapiens.UCSC.hg38, p.cutoff = 1e-04)
bg1 <- getBackgroundPeaks(object = rse1)
dev1 <- computeDeviations(object = rse1, annotations = motif_ix1, background_peaks = bg1)
variability1 <- computeVariability(dev1)
plotVariability(variability1, use_plotly = FALSE) 
variability1=left_join(variability1, motif.df, by=c("name"="V1"))
dev_df1=data.frame(t(dev1@assays@data@listData[["deviations"]]))
z_df1=data.frame(t(dev1@assays@data@listData[["z"]]))
write.table(dev_df1,gzfile(paste0(path22,"pen_aCRE_dev.tsv.gz")), col.names=T, row.names=T, sep="\t", quote=F)
write.table(z_df1,gzfile(paste0(path22,"pen_aCRE_z.tsv.gz")), col.names=T, row.names=T, sep="\t", quote=F)
write.table(variability1,gzfile(paste0(path22,"pen_aCRE_var.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
# promoter alone
combined.p=combined.pen
need.atac.p=which(rownames(combined.p[["aCRE"]]) %in% scACREv3a$peakID[which(scACREv3a$promoter_type_CT == "promoter-like")])
need.rna.p=which(rownames(combined.p[["tCRE"]]) %in% scACREv3a$peakID[which(scACREv3a$promoter_type_CT == "promoter-like")])
combined.p[["aCRE"]] <- subset(combined.p[["aCRE"]], features = rownames(combined.p[["aCRE"]])[need.atac.p])
combined.p[["tCRE"]] <- subset(combined.p[["tCRE"]], features = rownames(combined.p[["tCRE"]])[need.rna.p])

#===
# TSS
#remove peak that is zero across all cells
non_zero_frag_peaks <- which(rowSums(combined.p[["tCRE"]]) != 0)
combined.p[["tCRE"]] <- subset(combined.p[["tCRE"]], features = rownames(combined.p[["tCRE"]])[non_zero_frag_peaks])

peak_data <- combined.p[["tCRE"]]@counts
peak_ranges <- combined.p[["tCRE"]]@ranges
rse <- SummarizedExperiment(assays = list(peak_data),rowRanges = peak_ranges)
rse <- addGCBias(rse, genome = BSgenome.Hsapiens.UCSC.hg38)
names(rse@assays) <- c("counts")
motif_ix <- matchMotifs(pfm0, rse, genome = BSgenome.Hsapiens.UCSC.hg38, p.cutoff = 1e-04)
bg <- getBackgroundPeaks(object = rse)
dev <- computeDeviations(object = rse, annotations = motif_ix, background_peaks = bg)
variability <- computeVariability(dev)
plotVariability(variability, use_plotly = FALSE) 
variability=left_join(variability, motif.df, by=c("name"="V1"))
dev_df=data.frame(t(dev@assays@data@listData[["deviations"]]))
z_df=data.frame(t(dev@assays@data@listData[["z"]]))
write.table(dev_df,gzfile(paste0(path22,"p_tCRE_dev.tsv.gz")), col.names=T, row.names=T, sep="\t", quote=F)
write.table(z_df,gzfile(paste0(path22,"p_tCRE_z.tsv.gz")), col.names=T, row.names=T, sep="\t", quote=F)
write.table(variability,gzfile(paste0(path22,"p_tCRE_var.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#====
# ATAC
peak_data <- combined.p[["aCRE"]]@counts
peak_ranges <- combined.p[["aCRE"]]@ranges
rse1 <- SummarizedExperiment(assays = list(peak_data),rowRanges = peak_ranges)
rse1 <- addGCBias(rse1, genome = BSgenome.Hsapiens.UCSC.hg38)
names(rse1@assays) <- c("counts")
motif_ix1 <- matchMotifs(pfm0, rse1, genome = BSgenome.Hsapiens.UCSC.hg38, p.cutoff = 1e-04)
bg1 <- getBackgroundPeaks(object = rse1)
dev1 <- computeDeviations(object = rse1, annotations = motif_ix1, background_peaks = bg1)
variability1 <- computeVariability(dev1)
plotVariability(variability1, use_plotly = FALSE) 
variability1=left_join(variability1, motif.df, by=c("name"="V1"))
dev_df1=data.frame(t(dev1@assays@data@listData[["deviations"]]))
z_df1=data.frame(t(dev1@assays@data@listData[["z"]]))
write.table(dev_df1,gzfile(paste0(path22,"p_aCRE_dev.tsv.gz")), col.names=T, row.names=T, sep="\t", quote=F)
write.table(z_df1,gzfile(paste0(path22,"p_aCRE_z.tsv.gz")), col.names=T, row.names=T, sep="\t", quote=F)
write.table(variability1,gzfile(paste0(path22,"p_aCRE_var.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
#enhancer alone
combined.en=combined.pen
need.atac.en=which(rownames(combined.en[["aCRE"]]) %in% scACREv3a$peakID[which(scACREv3a$promoter_type_CT == "enhancer-like")])
need.rna.en=which(rownames(combined.en[["tCRE"]]) %in% scACREv3a$peakID[which(scACREv3a$promoter_type_CT == "enhancer-like")])
combined.en[["aCRE"]] <- subset(combined.en[["aCRE"]], features = rownames(combined.en[["aCRE"]])[need.atac.en])
combined.en[["tCRE"]] <- subset(combined.en[["tCRE"]], features = rownames(combined.en[["tCRE"]])[need.rna.en])

#====
# TSS
#remove peak that is zero across all cells
non_zero_frag_peaks <- which(rowSums(combined.en[["tCRE"]]) != 0)
combined.en[["tCRE"]] <- subset(combined.en[["tCRE"]], features = rownames(combined.en[["tCRE"]])[non_zero_frag_peaks])

peak_data <- combined.en[["tCRE"]]@counts
peak_ranges <- combined.en[["tCRE"]]@ranges
rse <- SummarizedExperiment(assays = list(peak_data),rowRanges = peak_ranges)
rse <- addGCBias(rse, genome = BSgenome.Hsapiens.UCSC.hg38)
names(rse@assays) <- c("counts")
motif_ix <- matchMotifs(pfm0, rse, genome = BSgenome.Hsapiens.UCSC.hg38, p.cutoff = 1e-04)
bg <- getBackgroundPeaks(object = rse)
dev <- computeDeviations(object = rse, annotations = motif_ix, background_peaks = bg)
variability <- computeVariability(dev)
plotVariability(variability, use_plotly = FALSE) 
variability=left_join(variability, motif.df, by=c("name"="V1"))
dev_df=data.frame(t(dev@assays@data@listData[["deviations"]]))
z_df=data.frame(t(dev@assays@data@listData[["z"]]))
write.table(dev_df,gzfile(paste0(path22,"en_tCRE_dev.tsv.gz")), col.names=T, row.names=T, sep="\t", quote=F)
write.table(z_df,gzfile(paste0(path22,"en_tCRE_z.tsv.gz")), col.names=T, row.names=T, sep="\t", quote=F)
write.table(variability,gzfile(paste0(path22,"en_tCRE_var.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#====
# ATAC
peak_data <- combined.en[["aCRE"]]@counts
peak_ranges <- combined.en[["aCRE"]]@ranges
rse1 <- SummarizedExperiment(assays = list(peak_data),rowRanges = peak_ranges)
rse1 <- addGCBias(rse1, genome = BSgenome.Hsapiens.UCSC.hg38)
names(rse1@assays) <- c("counts")
motif_ix1 <- matchMotifs(pfm0, rse1, genome = BSgenome.Hsapiens.UCSC.hg38, p.cutoff = 1e-04)
bg1 <- getBackgroundPeaks(object = rse1)
dev1 <- computeDeviations(object = rse1, annotations = motif_ix1, background_peaks = bg1)
variability1 <- computeVariability(dev1)
plotVariability(variability1, use_plotly = FALSE) 
variability1=left_join(variability1, motif.df, by=c("name"="V1"))
dev_df1=data.frame(t(dev1@assays@data@listData[["deviations"]]))
z_df1=data.frame(t(dev1@assays@data@listData[["z"]]))
write.table(dev_df1,gzfile(paste0(path22,"en_aCRE_dev.tsv.gz")), col.names=T, row.names=T, sep="\t", quote=F)
write.table(z_df1,gzfile(paste0(path22,"en_aCRE_z.tsv.gz")), col.names=T, row.names=T, sep="\t", quote=F)
write.table(variability1,gzfile(paste0(path22,"en_aCRE_var.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===============
#all aCRE with promoter/enhancer. -> tCRE can use pen.tCRE directly
scACREv3=read.delim(paste0(primary_folder,"Data_and_code/scATAC/merged_peak/all_aCRE.PE.final.tsv"), header=T, stringsAsFactors = F, check.names = F)
scACREv3a=scACREv3[which(scACREv3$promoter_type_CT != "unclassed"),]
scACREv3a=scACREv3a[which(scACREv3a$promoter_type_CT != "CTCF-alone"),]

combined.pen=combined 
need.atac=which(rownames(combined.pen[["aCRE"]]) %in% scACREv3a$peakID)
combined.pen[["aCRE"]] <- subset(combined.pen[["aCRE"]], features = rownames(combined.pen[["aCRE"]])[need.atac])

peak_data <- combined.pen[["aCRE"]]@counts
peak_ranges <- combined.pen[["aCRE"]]@ranges
rse1 <- SummarizedExperiment(assays = list(peak_data),rowRanges = peak_ranges)
rse1 <- addGCBias(rse1, genome = BSgenome.Hsapiens.UCSC.hg38)
names(rse1@assays) <- c("counts")
motif_ix1 <- matchMotifs(pfm0, rse1, genome = BSgenome.Hsapiens.UCSC.hg38, p.cutoff = 1e-04)
bg1 <- getBackgroundPeaks(object = rse1)
dev1 <- computeDeviations(object = rse1, annotations = motif_ix1, background_peaks = bg1)
variability1 <- computeVariability(dev1)
plotVariability(variability1, use_plotly = FALSE) 
variability1=left_join(variability1, motif.df, by=c("name"="V1"))
dev_df1=data.frame(t(dev1@assays@data@listData[["deviations"]]))
z_df1=data.frame(t(dev1@assays@data@listData[["z"]]))
write.table(dev_df1,gzfile(paste0(path22,"pen2_aCRE_dev.tsv.gz")), col.names=T, row.names=T, sep="\t", quote=F)
write.table(z_df1,gzfile(paste0(path22,"pen2_aCRE_z.tsv.gz")), col.names=T, row.names=T, sep="\t", quote=F)
write.table(variability1,gzfile(paste0(path22,"pen2_aCRE_var.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#result stored in paste0(chromvar_folder,"output/")

#===============================================================================


# Main version: separate enhancer promoter, transcribed CREs alone
# Joint version is also included as "Both"
# Joint version with all CREs is also included as "All"
#===============================================================================
# parse all the result
path22=paste0(chromvar_folder,"output/")

files=list.files(path=path22, pattern="var")
files.names=sapply(strsplit(files,"_var"),"[", 1)

data=data.frame()
for (i in 1:length(files)){
  data1=read.delim(paste0(path22, files[i]), header=T, stringsAsFactors = F, check.names = F)
  data1$set=files.names[i]
  data=rbind(data, data1)}
data=data[which(data$set != "pen2_aCRE"),]

data$set1=sapply(strsplit(data$set,"_"),"[",1)
data$set2=sapply(strsplit(data$set,"_"),"[",2)
data$set1=gsub("en","Enhancer",data$set1)
data$set1=gsub("p","Promoter",data$set1)
data$set1=gsub("PromoterEnhancer","Both",data$set1)
data$set1=gsub("all","All",data$set1)
data$set2=gsub("aCRE","ATAC-signal",data$set2)
data$set2=gsub("tCRE","TSS-signal",data$set2)

data$label=paste0(data$name,"_",data$V2)
data2 = data %>% group_by(set) %>% arrange(desc(variability), .by_group = TRUE) %>% mutate(rank = row_number()) %>% ungroup()
write.table(data2, gzfile(paste0(path_fig3_data,"variability_ChromVar_noFP.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
# -> for fig ex3b

#===============================================================================
# add pseudotime
motif.df6=read.delim(paste0(metacell,"output/TF_motif/motif.base.TF.expression.per.max02.cell3.meta.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
length(unique(motif.df6$motifID))

pseudotime1=read.delim(paste0(metacell,"output/pseudotime.scRNA.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
pseudotime1$cluster=sapply(strsplit(pseudotime1$SEACell,"_meta"),"[",1)
pseudotime1=pseudotime1%>%group_by(cluster)%>%dplyr::mutate(Tc=mean(T1))%>%arrange(Tc,T1)
pseudotime1$T1=1:nrow(pseudotime1)
pseudotime1$T2=ceiling(pseudotime1$T1/4)

#for trajectory 1
pseudotime2=pseudotime1[-grep("astrocyte",pseudotime1$SEACell),]
pseudotime2=pseudotime2[-grep("OPC",pseudotime2$SEACell),]
pseudotime2$T1=1:nrow(pseudotime2)
pseudotime2$T2=ceiling(pseudotime2$T1/3)

#for trajectory 2
pseudotime3=pseudotime1[-grep("euron",pseudotime1$SEACell),]
pseudotime3$T1=1:nrow(pseudotime3)
pseudotime3$T2=ceiling(pseudotime3$T1/3)

#===============================================================================
# gather z score of the 4 runs
path22=paste0(chromvar_folder,"output/")

aCRE.en=read.delim(paste0(path22,"en_aCRE_z.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
tCRE.en=read.delim(paste0(path22,"en_tCRE_z.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
aCRE.p=read.delim(paste0(path22,"p_aCRE_z.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
tCRE.p=read.delim(paste0(path22,"p_tCRE_z.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)

aCRE.en$metacell=rownames(aCRE.en)
aCRE.en=reshape2::melt(aCRE.en, id=880)
aCRE.en=aCRE.en%>%group_by(variable)%>%dplyr::mutate(variance=var(value))
aCRE.en$group1="ATAC-signal"
aCRE.en$group2="Enhancer"
tCRE.en$metacell=rownames(tCRE.en)
tCRE.en=reshape2::melt(tCRE.en, id=880)
tCRE.en=tCRE.en%>%group_by(variable)%>%dplyr::mutate(variance=var(value))
tCRE.en$group1="TSS-signal"
tCRE.en$group2="Enhancer"

aCRE.p$metacell=rownames(aCRE.p)
aCRE.p=reshape2::melt(aCRE.p, id=880)
aCRE.p=aCRE.p%>%group_by(variable)%>%dplyr::mutate(variance=var(value))
aCRE.p$group1="ATAC-signal"
aCRE.p$group2="Promoter"
tCRE.p$metacell=rownames(tCRE.p)
tCRE.p=reshape2::melt(tCRE.p, id=880)
tCRE.p=tCRE.p%>%group_by(variable)%>%dplyr::mutate(variance=var(value))
tCRE.p$group1="TSS-signal"
tCRE.p$group2="Promoter"

#==============
#trajectory 1 version
overall=rbind(aCRE.en,aCRE.p,tCRE.en,tCRE.p)
overall=overall[which(overall$metacell %in% pseudotime2$SEACell),]
overall=left_join(overall, pseudotime2[,c(1,4,6)], by=c("metacell"="SEACell"),copy=F)
overall=left_join(overall, motif.df6[,c(1,2,7,8,9)], by=c("metacell"="metacell", "variable"="motifID"), copy=F)
colnames(overall)[c(10,11)]=c("expression_TF1","expression_TF2")
overall=overall[which(overall$variable %in% motif.df6$motifID),] #limit to 400 motif

write.table(overall, gzfile(paste0(path_fig3_data,"plot.table.all.noFP.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#==============
#trajectory 2 version
overall=rbind(aCRE.en,aCRE.p,tCRE.en,tCRE.p)
overall=overall[which(overall$metacell %in% pseudotime3$SEACell),]
overall=left_join(overall, pseudotime3[,c(1,4,6)], by=c("metacell"="SEACell"),copy=F)
overall=left_join(overall, motif.df6[,c(1,2,7,8,9)], by=c("metacell"="metacell", "variable"="motifID"), copy=F)
colnames(overall)[c(10,11)]=c("expression_TF1","expression_TF2")
overall=overall[which(overall$variable %in% motif.df6$motifID),] #limit to 400 motif

write.table(overall, gzfile(paste0(path_fig3_data,"plot.table.all.trajectory2.noFP.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)


#===============================================================================
# perform correlation
path22=paste0(chromvar_folder,"output/")
assay=c("en_aCRE","en_tCRE","p_aCRE","p_tCRE","pen_aCRE","pen_tCRE","all_aCRE","all_tCRE")
g1=c("Enhancer","Enhancer","Promoter","Promoter","Both","Both","All","All")
g2=c("ATAC","TSS","ATAC","TSS","ATAC","TSS","ATAC","TSS")
overall.all=data.frame()

for (i in 1:length(assay)){
variability=read.delim(paste0(path22,assay[i],"_var.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data1=read.delim(paste0(path22,assay[i],"_z.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data1$metacell=rownames(data1)
data1=reshape2::melt(data1, id=880)
data1=left_join(data1, variability[,c(1,2,6,7)], by=c("variable"="name"),copy=F)
data1$group1=g1[i]
data1$group2=g2[i]
overall.all=rbind(overall.all,data1)}

overall.all=left_join(overall.all, pseudotime2[,c(1,4,6)], by=c("metacell"="SEACell"),copy=F)
overall.all=left_join(overall.all, motif.df6[,c(1,7:9)], by=c("metacell"="metacell", "variable"="motifID"), copy=F)
overall.all=overall.all[which(overall.all$variable %in% motif.df6$motifID),] # limit to 400 motifs
colnames(overall.all)[c(2,3,6,11,12)]=c("motifID","chromvar_z","motifName","expression_TF1","expression_TF2")

overall.all2=overall.all%>%group_by(group1, group2, motifID,motifName,variability,p_value_adj,T2)%>%dplyr::summarise(chromvar_z=mean(chromvar_z), TF1=mean(expression_TF1), TF2=mean(expression_TF2)) %>%dplyr::mutate (rel.TF1=TF1*max(chromvar_z)/max(TF1), rel.TF2=TF2*max(chromvar_z)/max(TF2))

#============================
#split TF1 and TF2 for correlation
overall.all3c=reshape2::melt(overall.all2[,c(1:10)], id=c(1:8))
colnames(overall.all3c)[c(9,10)]=c("TF","TF_expression")

#summary for chromVar variability cut off
overall.all3d=overall.all3c%>%group_by(group1,group2,motifID,motifName,variability,p_value_adj,TF)%>%dplyr::summarise(max_TF_exp=max(TF_expression))
overall.all3d$label=paste0(overall.all3d$motifID,"_",overall.all3d$motifName,"_",overall.all3d$TF)
overall.all3d$motif=paste0(overall.all3d$motifID,"_",overall.all3d$motifName)
overall.all3d=overall.all3d[which(!is.na(overall.all3d$max_TF_exp)),]
overall.all3d$label[-grep("\\:",overall.all3d$motif)]=overall.all3d$motif[-grep("\\:",overall.all3d$motif)]

overall.all3d$group3=paste0(overall.all3d$group1,"_",overall.all3d$group2)

overall.all3e=spread(overall.all3d[,c(11,9,10,8,5)],key=1, value=5)
overall.all3e1=spread(overall.all3d[,c(11,9,10,8,6)],key=1, value=5)
overall.all3e=left_join(overall.all3e,overall.all3e1[,c(1,4:11)], by="label", copy=F, suffix=c("_var","_padj"))
overall.all3e_ff=unique(overall.all3e[,c(2,4:11)])
overall.all3e$var_ATAC_all="no"
overall.all3e$var_ATAC_all[which(overall.all3e$All_ATAC_padj <0.01 & overall.all3e$All_ATAC_var>median(overall.all3e_ff$All_ATAC_var))]="yes"
overall.all3e$var_TSS_all="no"
overall.all3e$var_TSS_all[which(overall.all3e$All_TSS_padj <0.01 & overall.all3e$All_TSS_var>median(overall.all3e_ff$All_TSS_var))]="yes"

overall.all3e$var_ATAC_both="no"
overall.all3e$var_ATAC_both[which(overall.all3e$Both_ATAC_padj <0.01 & overall.all3e$Both_ATAC_var>median(overall.all3e_ff$Both_ATAC_var))]="yes"
overall.all3e$var_TSS_both="no"
overall.all3e$var_TSS_both[which(overall.all3e$Both_TSS_padj <0.01 & overall.all3e$Both_TSS_var>median(overall.all3e_ff$Both_TSS_var))]="yes"

overall.all3e$var_ATAC_enhancer="no"
overall.all3e$var_ATAC_enhancer[which(overall.all3e$Enhancer_ATAC_padj <0.01 & overall.all3e$Enhancer_ATAC_var>median(overall.all3e_ff$Enhancer_ATAC_var))]="yes"
overall.all3e$var_TSS_enhancer="no"
overall.all3e$var_TSS_enhancer[which(overall.all3e$Enhancer_TSS_padj <0.01 & overall.all3e$Enhancer_TSS_var>median(overall.all3e_ff$Enhancer_TSS_var))]="yes"

overall.all3e$var_ATAC_promoter="no"
overall.all3e$var_ATAC_promoter[which(overall.all3e$Promoter_ATAC_padj <0.01 & overall.all3e$Promoter_ATAC_var>median(overall.all3e_ff$Promoter_ATAC_var))]="yes"
overall.all3e$var_TSS_promoter="no"
overall.all3e$var_TSS_promoter[which(overall.all3e$Promoter_TSS_padj <0.01 & overall.all3e$Promoter_TSS_var>median(overall.all3e_ff$Promoter_TSS_var))]="yes"

colSums(overall.all3e[,c(20:27)] == "yes")
#-> all yes according to the cut off max>=0.2 (per meta cell), >=3 meta cells, require both TF

#====================
#get the correlation for all the expressed TF first, filter later
length(unique(overall.all3c$motifID)) #400
overall.all3c$motif=paste0(overall.all3c$motifID,"_",overall.all3c$motifName)
overall.all3c$label=paste0(overall.all3c$motif,"_",overall.all3c$TF)
overall.all3c$label[-grep("\\:",overall.all3c$motif)]=overall.all3c$motif[-grep("\\:",overall.all3c$motif)]
length(unique(overall.all3c$label)) #423
write.table(overall.all3c,gzfile(paste0(path_fig3_data,"chromvar_z_expression_400motif_noFP.tsv.gz")), col.names = T, row.names = F, sep="\t", quote=F)
overall.all3c=read.delim(paste0(path_fig3_data,"chromvar_z_expression_400motif_noFP.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

#===============================================================================
#for the 400 motifs, expression level of all the TFs are present, metacell without expression -> 0

# only include 400 motif with expression support
together2=overall.all3c[which(!is.na(overall.all3c$TF_expression)),]%>%group_by(group1,group2,label,motifID,motifName,TF)%>%dplyr::summarise(acf_exVall_max=max(ccf(ts(chromvar_z),ts(TF_expression), 4, type="correlation")$acf[5:9]),
                                                                                                                                             acf_lag_exVall_max=which.max(ccf(ts(chromvar_z),ts(TF_expression), 4, type="correlation")$acf[5:9]) +4 - 5,               
                                                                                                                                             acf_exVall_min=min(ccf(ts(chromvar_z),ts(TF_expression), 4, type="correlation")$acf[5:9]),
                                                                                                                                             acf_lag_exVall_min=which.min(ccf(ts(chromvar_z),ts(TF_expression), 4, type="correlation")$acf[5:9]) +4 - 5)
together2$acf_exVall_opt=together2$acf_exVall_max
together2$acf_lag_exVall_opt=together2$acf_lag_exVall_max
together2$acf_exVall_opt[which(abs(together2$acf_exVall_max)<abs(together2$acf_exVall_min))]=together2$acf_exVall_min[which(abs(together2$acf_exVall_max)<abs(together2$acf_exVall_min))]
together2$acf_lag_exVall_opt[which(abs(together2$acf_exVall_max)<abs(together2$acf_exVall_min))]=together2$acf_lag_exVall_min[which(abs(together2$acf_exVall_max)<abs(together2$acf_exVall_min))]

together2=left_join(together2, unique(overall.all3c[,c("group1","group2","label","TF","variability","p_value_adj")]), by=c("group1","group2","label","TF"),copy=F)
colnames(together2)[c(13,14)]=c("chromVar_variability","chromVar_p_value_adj")

together2$group2[which(together2$group2=="ATAC")]="ATAC-signal"
together2$group2[which(together2$group2=="TSS")]="TSS-signal"

data2=read.delim(paste0(path_fig3_data,"variability_ChromVar_noFP.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data2a=data2[which(data2$name %in% unique(together2$motifID)),]
data2a=data2a%>%group_by(set)%>%dplyr::mutate(scaled_var=variability/max(variability))
data2a=data2a%>%group_by(set)%>%dplyr::mutate(zscaled_var=(variability-mean(variability))/sd(variability))
data2a=data2a%>%group_by(set)%>%dplyr::mutate(median_var=median(variability), median_scaled_var=median(scaled_var), median_zscaled_var=median(zscaled_var))

together2=left_join(together2, data2a[,c("name","set1","set2","median_var","zscaled_var","median_zscaled_var")], by=c("group1"="set1","group2"="set2","motifID"="name"),copy=F)
together2$chromVar_Var="no"
together2$chromVar_Var[which(together2$chromVar_p_value_adj <0.01 & together2$chromVar_variability>together2$median_var)]="yes" #use median of 400 motif as cutoff
together2$pass_both = "No"
together2$pass_both[which(abs(together2$acf_exVall_opt) >= 0.5 & together2$chromVar_Var == "yes")] = "Yes"

write.table(together2,gzfile(paste0(path_fig3_data,"auto_correlation_to_p_en.z.summarise_noFP.tsv.gz")), col.names = T, row.names = F, sep="\t", quote=F)
# -> for fig ex3c

#===============================================================================
together2=read.delim(paste0(path_fig3_data,"auto_correlation_to_p_en.z.summarise_noFP.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
#together2=together2%>%group_by(group1,group2,motifID)%>%slice_max(abs(acf_exVall_opt))
together2b <- together2 %>%pivot_wider(id_cols = label, names_from = c(group1, group2), values_from = c(acf_exVall_opt), names_glue = "{group1}_{group2}")
together2c <- together2 %>%pivot_wider(id_cols = label, names_from = c(group1, group2), values_from = c(pass_both), names_glue = "{group1}_{group2}")

colnames(together2b)[c(2:9)]=paste0(colnames(together2b)[c(2:9)],"_corr")
colnames(together2b)[c(2:9)]=gsub("-signal","",colnames(together2b)[c(2:9)])
colnames(together2c)[c(2:9)]=paste0("varNcorr05_",colnames(together2c)[c(2:9)])
colnames(together2c)[c(2:9)]=gsub("-signal","",colnames(together2c)[c(2:9)])

outputfile=left_join(overall.all3e, together2b, by="label", copy=F)
outputfile=left_join(outputfile, together2c, by="label", copy=F)
colSums(outputfile[,c(36:43)]=="Yes")
length(unique(outputfile$motif[which(rowSums(outputfile[,c(40:43)] == "Yes")>0)])) #247

outputfile$motifID=sapply(strsplit(outputfile$motif,"_"),"[",1)
outputfilea=outputfile[-grep("TF2",outputfile$label),]
outputfilea=left_join(outputfilea,unique(motif.df6[,c(1,5)]), by="motifID", copy=F)
outputfileb=outputfile[grep("TF2",outputfile$label),]
outputfileb=left_join(outputfileb,unique(motif.df6[,c(1,6)]), by="motifID", copy=F)
colnames(outputfilea)[45]="TF_symbol"
colnames(outputfileb)[45]="TF_symbol"
outputfile=rbind(outputfilea,outputfileb)
write.table(outputfile,gzfile(paste0(path_fig3_data,"variance_corr_p_en_400motif.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
#prepare literature review on known TF function, for FE later

motif.df6=read.delim(paste0(metacell,"output/TF_motif/motif.base.TF.expression.per.max02.cell3.meta.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
TF_involve=union(unique(motif.df6$Approved.symbol.TF1), unique(motif.df6$Approved.symbol.TF2)) #343 TF, 400 motif
motif.df6d=unique(motif.df6[,c(1:6)])
motif.df6e1=unique(motif.df6[,c(5,7,8)])
motif.df6e2=unique(motif.df6[which(!is.na(motif.df6$motifName2)),c(6,7,9)])
colnames(motif.df6e2)=colnames(motif.df6e1)
motif.df6e=unique(rbind(motif.df6e1,motif.df6e2))
motif.df6f=spread(motif.df6e, key=1, value=3)
write.table(motif.df6d,gzfile(paste0(path_fig3_data,"expressed_400motif.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(motif.df6f,gzfile(paste0(path_fig3_data,"expressed_342TFx391meta.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

motif.df7=read.delim(paste0(metacell,"output/TF_motif/motif.info.with_TF.expression.repress.max02.tsv"),header=T, stringsAsFactors = F, check.names = F)
motif.df7$label=paste0(motif.df7$V1,"_",motif.df7$V2)
motif.df7$label2=gsub("motifName","TF",motif.df7$variable)
motif.df7$label[grep("::",motif.df7$label)]=paste0(motif.df7$label[grep("::",motif.df7$label)],"_",motif.df7$label2[grep("::",motif.df7$label)])

close=c("REST","EZH2","ZNF217","KLF3","BACH1","FOXG1","MAFG","TGIF2","INSM1")

motif.df8=read.delim(paste0(metacell,"output/TF_motif/motif.info.with_TF.expression.tsv"),header=T, stringsAsFactors = F, check.names = F)
motif.df8$label=paste0(motif.df8$V1,"_",motif.df8$V2)
motif.df8$label2=gsub("motifName","TF",motif.df8$variable)
motif.df8$label[grep("::",motif.df8$label)]=paste0(motif.df8$label[grep("::",motif.df8$label)],"_",motif.df8$label2[grep("::",motif.df8$label)])

motif.df8$TF_activator="No"
motif.df8$TF_activator[grep("activator",motif.df8$TF_activity)]="Yes"
motif.df8$TF_activator[grep("activator",motif.df8$TF_activity2)]="Yes"
motif.df8$TF_repressor="No"
motif.df8$TF_repressor[grep("repressor",motif.df8$TF_activity)]="Yes"
motif.df8$TF_repressor[grep("repressor",motif.df8$TF_activity2)]="Yes"
motif.df8$TF_pioneer="No"
motif.df8$TF_pioneer[grep("yes",motif.df8$`pioneer factor`)]="Yes"
motif.df8$TF_pioneer[grep("maybe",motif.df8$`pioneer factor`)]="Yes"
motif.df8$TF_closer="No"
motif.df8$TF_closer[which(motif.df8$Approved.symbol %in% close)]="Yes"

motif.df8%>%group_by(TF_activator,TF_repressor,TF_pioneer,TF_closer)%>%dplyr::summarise(count=n())

motif.df8a=unique(motif.df8[,c(1,2,29:32)])%>%group_by(V1,V2)%>%dplyr::summarise(
  TF_activator=paste(unique(TF_activator),collapse=";"),TF_repressor=paste(unique(TF_repressor),collapse=";"),TF_pioneer=paste(unique(TF_pioneer),collapse=";"),TF_closer=paste(unique(TF_closer),collapse=";"))
motif.df8a$TF_activator[grep("Yes",motif.df8a$TF_activator)]="Yes"
motif.df8a$TF_activator[grep(";",motif.df8a$TF_activator)]="No"
motif.df8a$TF_repressor[grep("Yes",motif.df8a$TF_repressor)]="Yes"
motif.df8a$TF_repressor[grep(";",motif.df8a$TF_repressor)]="No"
motif.df8a$TF_pioneer[grep("Yes",motif.df8a$TF_pioneer)]="Yes"
motif.df8a$TF_pioneer[grep(";",motif.df8a$TF_pioneer)]="No"
motif.df8a$TF_closer[grep("Yes",motif.df8a$TF_closer)]="Yes"
motif.df8a$TF_closer[grep(";",motif.df8a$TF_closer)]="No"
motif.df8a%>%group_by(TF_activator,TF_repressor,TF_pioneer)%>%dplyr::summarise(count=n())

motif.df8b=motif.df8a[which(motif.df8a$V1 %in% motif.df6$motifID),]
write.table(motif.df8b,gzfile(paste0(path_fig3_data,"motif.info.with_TF.expression_n400.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
# fixed cutoff file
motif.df8b=read.delim(paste0(path_fig3_data,"motif.info.with_TF.expression_n400.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
outputfile=read.delim(paste0(path_fig3_data,"variance_corr_p_en_400motif.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
outputfile=left_join(outputfile,motif.df8b, by=c("motifID"="V1"),copy=F)

#revise colnames
colnames(outputfile)[46]="motifName"

outputfile$'ATAC+'="No"
outputfile$'ATAC+'[which(outputfile$Promoter_ATAC_corr >= 0.5 & outputfile$varNcorr05_Promoter_ATAC == "Yes")]="Yes"
outputfile$'ATAC+'[which(outputfile$Enhancer_ATAC_corr >= 0.5 & outputfile$varNcorr05_Enhancer_ATAC == "Yes")]="Yes"

outputfile$'ATAC-'="No"
outputfile$'ATAC-'[which(outputfile$Promoter_ATAC_corr <= (-0.5) & outputfile$varNcorr05_Promoter_ATAC == "Yes")]="Yes"
outputfile$'ATAC-'[which(outputfile$Enhancer_ATAC_corr <= (-0.5) & outputfile$varNcorr05_Enhancer_ATAC == "Yes")]="Yes"

outputfile$'TSS+'="No"
outputfile$'TSS+'[which(outputfile$Promoter_TSS_corr >= 0.5 & outputfile$varNcorr05_Promoter_TSS == "Yes")]="Yes"
outputfile$'TSS+'[which(outputfile$Enhancer_TSS_corr >= 0.5 & outputfile$varNcorr05_Enhancer_TSS == "Yes")]="Yes"

outputfile$'TSS-'="No"
outputfile$'TSS-'[which(outputfile$Promoter_TSS_corr <= (-0.5) & outputfile$varNcorr05_Promoter_TSS == "Yes")]="Yes"
outputfile$'TSS-'[which(outputfile$Enhancer_TSS_corr <= (-0.5) & outputfile$varNcorr05_Enhancer_TSS == "Yes")]="Yes"

#select representative TF1 or TF2, choose TF1 or 2 for overall result
outputfile <- outputfile %>% group_by(motifID) %>% dplyr::slice_max(order_by = rowSums(across(`ATAC+`:`TSS-`) == "Yes"), with_ties = FALSE)

outputfile1=outputfile[which(rowSums(outputfile[,c(51:54)] == "Yes")>0),]
write.table(outputfile,gzfile(paste0(path_fig3_data,"TF_grouping_by_chromvar_correlation_and_literature.n400.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(outputfile1,gzfile(paste0(path_fig3_data,"TF_grouping_by_chromvar_correlation_and_literature.n247.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
# -> for fig 3f

#===============================================================================
#modal specific motif, with different cutoffs
together2=read.delim(paste0(path_fig3_data,"auto_correlation_to_p_en.z.summarise_noFP.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
together2=together2[which(together2$group1 %in% c("Promoter","Enhancer")),]
together2$effect="Null"
together2$effect[which(together2$acf_exVall_opt>=0.5 & together2$chromVar_variability>together2$median_var & together2$chromVar_p_value_adj < 0.01)]="Pos"
together2$effect[which(together2$acf_exVall_opt<=(-0.5) & together2$chromVar_variability>together2$median_var & together2$chromVar_p_value_adj < 0.01)]="Neg"
together2$group3=paste0(together2$group2,"_",together2$group1)
together2k=spread(together2[,c(4:6,21,20)],key=4, value=5)
#select representative TF1 or TF2, choose TF1 or 2 for overall result
together2k <- together2k %>% group_by(motifID) %>% dplyr::slice_max(order_by = rowSums(across(`ATAC-signal_Enhancer`:`TSS-signal_Promoter`) != "Null"), with_ties = FALSE)

together2k$en_match="No"
together2k$en_match[which(together2k$`ATAC-signal_Enhancer` == together2k$`TSS-signal_Enhancer` & together2k$`ATAC-signal_Enhancer`!="Null")]="Yes"
together2k$p_match="No"
together2k$p_match[which(together2k$`ATAC-signal_Promoter` == together2k$`TSS-signal_Promoter` & together2k$`ATAC-signal_Promoter`!="Null")]="Yes"
together2k$cross_match="No"
together2k$cross_match[which(together2k$`ATAC-signal_Enhancer` == together2k$`TSS-signal_Promoter` & together2k$`ATAC-signal_Enhancer`!="Null")]="Yes"
together2k$cross_match[which(together2k$`ATAC-signal_Promoter` == together2k$`TSS-signal_Enhancer` & together2k$`ATAC-signal_Promoter`!="Null")]="Yes"
length(which(together2k$en_match == "Yes" | together2k$p_match == "Yes" | together2k$cross_match == "Yes")) #77
length(which(rowSums(together2k[,c(4:7)] == "Null")<4)) #247

corr_cutoff=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)
var_cutoff=c(0.25,0.5,0.75)
data6=data.frame(matrix(ncol=4,nrow=24))
colnames(data6)=c("var_cutoff","corr_cutoff","union","coupled")
for (j in 1:3){
  together2=together2%>%group_by(group3)%>%mutate(cutoff=quantile(chromVar_variability[!duplicated(motifID)], var_cutoff[j])) #keep only 400 unique variability in calculating the median
  for (i in 1:8){
    together2$effect="Null"
    together2$effect[which(together2$acf_exVall_opt>=corr_cutoff[i] & together2$chromVar_variability>together2$cutoff & together2$chromVar_p_value_adj < 0.01)]="Pos"
    together2$effect[which(together2$acf_exVall_opt<=(-corr_cutoff[i]) & together2$chromVar_variability>together2$cutoff & together2$chromVar_p_value_adj < 0.01)]="Neg"
    together2k=spread(together2[,c(4:6,21,20)],key=4, value=5)
    together2k <- together2k %>% group_by(motifID) %>% dplyr::slice_max(order_by = rowSums(across(`ATAC-signal_Enhancer`:`TSS-signal_Promoter`) != "Null"), n=1, with_ties = FALSE)
    together2k$en_match="No"
    together2k$en_match[which(together2k$`ATAC-signal_Enhancer` == together2k$`TSS-signal_Enhancer` & together2k$`ATAC-signal_Enhancer`!="Null")]="Yes"
    together2k$p_match="No"
    together2k$p_match[which(together2k$`ATAC-signal_Promoter` == together2k$`TSS-signal_Promoter` & together2k$`ATAC-signal_Promoter`!="Null")]="Yes"
    together2k$cross_match="No"
    together2k$cross_match[which(together2k$`ATAC-signal_Enhancer` == together2k$`TSS-signal_Promoter` & together2k$`ATAC-signal_Enhancer`!="Null")]="Yes"
    together2k$cross_match[which(together2k$`ATAC-signal_Promoter` == together2k$`TSS-signal_Enhancer` & together2k$`ATAC-signal_Promoter`!="Null")]="Yes"
    data6$var_cutoff[(j-1)*8+i] = var_cutoff[j]
    data6$corr_cutoff[(j-1)*8+i] = corr_cutoff[i]
    data6$union[(j-1)*8+i]=length(which(rowSums(together2k[,c(4:7)] == "Null")<4))
    data6$coupled[(j-1)*8+i]=length(which(together2k$en_match == "Yes" | together2k$p_match == "Yes" | together2k$cross_match == "Yes"))
  }}
data6$modal_specific=(data6$union-data6$coupled)/data6$union
data6$var_cutoff=paste0("Var > P",data6$var_cutoff*100)
data6$var_cutoff=gsub("P50","Median",data6$var_cutoff)
data6$var_cutoff=factor(data6$var_cutoff, levels=c("Var > P25","Var > Median","Var > P75"))

write.table(data6,gzfile(paste0(path_fig3_data,"ATAC_RNA.venn.with.cuttoff.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
# -> for fig ex3f

#===============================================================================
# for overall FE with literature with multiple cutoff
together2=read.delim(paste0(path_fig3_data,"auto_correlation_to_p_en.z.summarise_noFP.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
together2=together2[which(together2$group1 %in% c("Promoter","Enhancer")),]
together2$effect="Null"
together2$effect[which(together2$acf_exVall_opt>=0.5 & together2$chromVar_variability>together2$median_var & together2$chromVar_p_value_adj < 0.01)]="Pos"
together2$effect[which(together2$acf_exVall_opt<=(-0.5) & together2$chromVar_variability>together2$median_var & together2$chromVar_p_value_adj < 0.01)]="Neg"
together2$group3=paste0(together2$group2,"_",together2$group1)

motif.df8b=read.delim(paste0(path_fig3_data,"motif.info.with_TF.expression_n400.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

corr_cutoff=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)
var_cutoff=c(0.25,0.5,0.75)
data7=data.frame()
for (j in 1:3){
  together2=together2%>%group_by(group3)%>%mutate(cutoff=quantile(chromVar_variability[!duplicated(motifID)], var_cutoff[j])) #keep only 400 unique variability in calculating the median
  for (i in 1:8){
    together2$effect="Null"
    together2$effect[which(together2$acf_exVall_opt>=corr_cutoff[i] & together2$chromVar_variability>together2$cutoff & together2$chromVar_p_value_adj < 0.01)]="Pos"
    together2$effect[which(together2$acf_exVall_opt<=(-corr_cutoff[i]) & together2$chromVar_variability>together2$cutoff & together2$chromVar_p_value_adj < 0.01)]="Neg"
    together2k=spread(together2[,c(4:6,21,20)],key=4, value=5)
    together2k <- together2k %>% group_by(motifID) %>% dplyr::slice_max(order_by = rowSums(across(`ATAC-signal_Enhancer`:`TSS-signal_Promoter`) != "Null"), n=1, with_ties = FALSE)
    together2k=left_join(together2k,motif.df8b, by=c("motifID"="V1"),copy=F)
    together2k$`ATAC+`="No"
    together2k$`ATAC+`[which(together2k$`ATAC-signal_Enhancer`=="Pos" | together2k$`ATAC-signal_Promoter` == "Pos")]="Yes"
    together2k$`ATAC-`="No"
    together2k$`ATAC-`[which(together2k$`ATAC-signal_Enhancer`=="Neg" | together2k$`ATAC-signal_Promoter` == "Neg")]="Yes"
    together2k$`TSS+`="No"
    together2k$`TSS+`[which(together2k$`TSS-signal_Enhancer`=="Pos" | together2k$`TSS-signal_Promoter` == "Pos")]="Yes"
    together2k$`TSS-`="No"
    together2k$`TSS-`[which(together2k$`TSS-signal_Enhancer`=="Neg" | together2k$`TSS-signal_Promoter` == "Neg")]="Yes"
    outputfile2=together2k[which(rowSums(together2k[,c(13:16)] == "No")<4),] 
    outputfile2=outputfile2[which(rowSums(outputfile2[,c(9:12)] == "No")<4),] 
    FE1=outputfile2[,c("motifID","ATAC+","ATAC-","TSS+","TSS-")]
    FE1[is.na(FE1)]="No"
    FEsum=reshape2::melt(FE1,id=1)
    FEsum2=reshape2::melt(outputfile2[,c("motifID","TF_activator","TF_repressor","TF_pioneer","TF_closer")], id=1)
    FEsum=left_join(FEsum,FEsum2, by="motifID", copy=F, relationship = "many-to-many")
    FEsum2=FEsum%>%group_by(variable.x,variable.y,value.x, value.y)%>%dplyr::summarise(count=n())
    FEsum2$key=paste0(FEsum2$value.x,"_",FEsum2$value.y)
    FEsum3=spread(FEsum2[,c(1,2,5,6)], key=4, value=3)
    FEsum3[is.na(FEsum3)]=0
    for(k in 1:nrow(FEsum3)){
      GSEATasting <- matrix(c(FEsum3$No_No[k], FEsum3$No_Yes[k], FEsum3$Yes_No[k], FEsum3$Yes_Yes[k]), nrow = 2, dimnames = list(oligo1 = c("no", "yes"), oligo2 = c("no", "yes")))
      FEsum3$p.val[k] = fisher.test(GSEATasting, alternative = "two.sided")$p.value
      FEsum3$OR[k] = fisher.test(GSEATasting, alternative = "two.sided")$estimate}
    FEsum3$var_cutoff=var_cutoff[j]
    FEsum3$corr_cutoff=corr_cutoff[i]
    data7=rbind(data7, FEsum3)
  }}

data7$logOR=log(data7$OR)
data7$label="***"
data7$label[which(data7$p.val>=0.001)]="**"
data7$label[which(data7$p.val>=0.01)]="*"
data7$label[which(data7$p.val>=0.05)]="n.s."
data7$label1=paste0("log(OR)=",signif(data7$logOR,2),data7$label)

data7$var_cutoff=paste0("Var > P",data7$var_cutoff*100)
data7$var_cutoff=gsub("P50","Median",data7$var_cutoff)
data7$var_cutoff=factor(data7$var_cutoff, levels=c("Var > P25","Var > Median","Var > P75"))

write.table(data7,gzfile(paste0(path_fig3_data,"TF_grouping_by_chromvar_correlation_and_literature.FE.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
# -> for fig 3h

#===============================================================================
#supplementary heatmap showing all cell types (n=10)
#take the mean scaled zscore of each cluster
#z is a matrix -> z score here is across the meta cell
#-> need to apply another zscore across signal

path22=paste0(chromvar_folder,"output/")
cluster_info=read.delim(paste0(primary_folder,"Fig1/cluster_info.tsv"),header=T, stringsAsFactors = F, check.names = F)
files=list.files(path=path22, pattern="_z")
files.names=sapply(strsplit(files,"_z"),"[", 1)

motif.df=read.delim(paste0(chromvar_folder,"jaspar2024/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.ID.tsv"), header=F, stringsAsFactors = F)

ddata=data.frame()
for (i in 1:length(files)){
  ddata1=read.delim(paste0(path22, files[i]), header=T, stringsAsFactors = F, check.names = F)
  ddata1$cluster=sapply(strsplit(rownames(ddata1),"_meta"),"[",1)
  ddata2=reshape2::melt(ddata1, id=880)
  ddata3=ddata2%>%group_by(cluster, variable)%>%dplyr::summarise(value=mean(value))
  ddata3=left_join(ddata3, motif.df, by=c("variable"="V1"), copy=F)
  ddata3$V2=paste0(ddata3$variable,"_",ddata3$V2)
  ddata3=left_join(ddata3, cluster_info[,c(1,3)], by="cluster", copy=F)
  ddata3$set=files.names[i]
  ddata=rbind(ddata, ddata3)}
ddata=ddata[which(ddata$set %in% c("en_aCRE","en_tCRE","p_aCRE","p_tCRE")),]
ddata$label2=paste0(ddata$label,"_",ddata$set)

outputfile1=read.delim(paste0(path_fig3_data,"TF_grouping_by_chromvar_correlation_and_literature.n247.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
ddata=ddata[which(ddata$V2 %in% outputfile1$motif),]
ddata$group=sapply(strsplit(ddata$set,"_"),"[",2)
ddata=ddata%>%group_by(group)%>%dplyr::mutate(scaledz=(value - mean(value, na.rm=TRUE)) / sd(value, na.rm=TRUE))
#normalized value between aCRE and tCRE

write.table(ddata,gzfile(paste0(path_fig3_data,"heatmap_z_247motif.plot.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
# -> for fig ex3d

# select the best TF from each cluster
ddatb=ddata%>%group_by(label)%>%slice_max(order_by = scaledz, n = 2, with_ties = FALSE)
need=c("MA0142.1_Pou5f1::Sox2","MA1122.2_TFDP1","MA0839.2_CREB3L1","MA0489.3_Jun","MA1109.2_NEUROD1","MA0669.1_NEUROG2","MA0756.3_ONECUT2","MA0491.3_JUND","MA0809.3_TEAD4","MA1120.2_SOX13","MA0103.4_ZEB1","MA1644.2_NFYC")
ddata1=ddata[which(ddata$V2 %in% need),]
write.table(ddata1,gzfile(paste0(path_fig3_data,"heatmap_z_12motif.plot.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
# -> for fig 3d

#===============================================================================
# joint promoter-enhancer ChromVar result
#===============================================================================
# repeat data6 using the promoter enhancer combined version

together2=read.delim(paste0(path_fig3_data,"auto_correlation_to_p_en.z.summarise_noFP.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
together2=together2[which(together2$group1 %in% c("Both")),]
together2$effect="Null"
together2$effect[which(together2$acf_exVall_opt>=0.5 & together2$chromVar_variability>together2$median_var & together2$chromVar_p_value_adj < 0.01)]="Pos"
together2$effect[which(together2$acf_exVall_opt<=(-0.5) & together2$chromVar_variability>together2$median_var & together2$chromVar_p_value_adj < 0.01)]="Neg"
together2$group3=paste0(together2$group2,"_",together2$group1)
together2k=spread(together2[,c(4:6,21,20)],key=4, value=5)
together2k <- together2k %>% group_by(motifID) %>% dplyr::slice_max(order_by = rowSums(across(`ATAC-signal_Both`:`TSS-signal_Both`) != "Null"), with_ties = FALSE)

together2k$match="No"
together2k$match[which(together2k$`ATAC-signal_Both` == together2k$`TSS-signal_Both` & together2k$`ATAC-signal_Both`!="Null")]="Yes"
length(which(together2k$match == "Yes" )) #33
length(which(rowSums(together2k[,c(4:5)] == "Null")<2)) #181

corr_cutoff=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)
var_cutoff=c(0.25,0.5,0.75)
data6=data.frame(matrix(ncol=4,nrow=24))
colnames(data6)=c("var_cutoff","corr_cutoff","union","coupled")
for (j in 1:3){
  together2=together2%>%group_by(group3)%>%mutate(cutoff=quantile(chromVar_variability[!duplicated(motifID)], var_cutoff[j])) #keep only 400 unique variability in calculating the median
  for (i in 1:8){
    together2$effect="Null"
    together2$effect[which(together2$acf_exVall_opt>=corr_cutoff[i] & together2$chromVar_variability>together2$cutoff & together2$chromVar_p_value_adj < 0.01)]="Pos"
    together2$effect[which(together2$acf_exVall_opt<=(-corr_cutoff[i]) & together2$chromVar_variability>together2$cutoff & together2$chromVar_p_value_adj < 0.01)]="Neg"
    together2k=spread(together2[,c(4:6,21,20)],key=4, value=5)
    together2k <- together2k %>% group_by(motifID) %>% dplyr::slice_max(order_by = rowSums(across(`ATAC-signal_Both`:`TSS-signal_Both`) != "Null"), with_ties = FALSE)
    together2k$match="No"
    together2k$match[which(together2k$`ATAC-signal_Both` == together2k$`TSS-signal_Both` & together2k$`ATAC-signal_Both`!="Null")]="Yes"
    data6$var_cutoff[(j-1)*8+i] = var_cutoff[j]
    data6$corr_cutoff[(j-1)*8+i] = corr_cutoff[i]
    data6$union[(j-1)*8+i]=length(which(rowSums(together2k[,c(4:5)] == "Null")<2))
    data6$coupled[(j-1)*8+i]=length(which(together2k$match == "Yes" ))
  }}
data6$modal_specific=(data6$union-data6$coupled)/data6$union
data6$var_cutoff=paste0("Var > P",data6$var_cutoff*100)
data6$var_cutoff=gsub("P50","Median",data6$var_cutoff)
data6$var_cutoff=factor(data6$var_cutoff, levels=c("Var > P25","Var > Median","Var > P75"))
write.table(data6,gzfile(paste0(path_fig3_data,"ATAC_RNA.venn.with.cuttoff_joint_pen.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
# -> for fig ex3g

#==============================================
motif.df8b=read.delim(paste0(path_fig3_data,"motif.info.with_TF.expression_n400.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
outputfile=read.delim(paste0(path_fig3_data,"variance_corr_p_en_400motif.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
outputfile=left_join(outputfile,motif.df8b, by=c("motifID"="V1"),copy=F)

#revise colnames
colnames(outputfile)[46]="motifName"

outputfile$'ATAC+'="No"
outputfile$'ATAC+'[which(outputfile$Both_ATAC_corr >= 0.5 & outputfile$varNcorr05_Both_ATAC == "Yes")]="Yes"

outputfile$'ATAC-'="No"
outputfile$'ATAC-'[which(outputfile$Both_ATAC_corr <= (-0.5) & outputfile$varNcorr05_Both_ATAC == "Yes")]="Yes"

outputfile$'TSS+'="No"
outputfile$'TSS+'[which(outputfile$Both_TSS_corr >= 0.5 & outputfile$varNcorr05_Both_TSS == "Yes")]="Yes"

outputfile$'TSS-'="No"
outputfile$'TSS-'[which(outputfile$Both_TSS_corr <= (-0.5) & outputfile$varNcorr05_Both_TSS == "Yes")]="Yes"

#select representative TF1 or TF2, choose TF1 or 2 for overall result
outputfile <- outputfile %>% group_by(motifID) %>% dplyr::slice_max(order_by = rowSums(across(`ATAC+`:`TSS-`) == "Yes"), with_ties = FALSE)

outputfile1=outputfile[which(rowSums(outputfile[,c(51:54)] == "Yes")>0),]
write.table(outputfile1,gzfile(paste0(path_fig3_data,"TF_grouping_by_joint_pen_chromvar_correlation_and_literature.n181.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
# -> for fig ex3f

#===============================================================================
together2=read.delim(paste0(path_fig3_data,"auto_correlation_to_p_en.z.summarise_noFP.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
together2=together2[which(together2$group1 %in% c("Both")),]
together2$effect="Null"
together2$effect[which(together2$acf_exVall_opt>=0.5 & together2$chromVar_variability>together2$median_var & together2$chromVar_p_value_adj < 0.01)]="Pos"
together2$effect[which(together2$acf_exVall_opt<=(-0.5) & together2$chromVar_variability>together2$median_var & together2$chromVar_p_value_adj < 0.01)]="Neg"
together2$group3=paste0(together2$group2,"_",together2$group1)

motif.df8b=read.delim(paste0(path_fig3_data,"motif.info.with_TF.expression_n400.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

corr_cutoff=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)
var_cutoff=c(0.25,0.5,0.75)
data7=data.frame()
for (j in 1:3){
  together2=together2%>%group_by(group3)%>%mutate(cutoff=quantile(chromVar_variability[!duplicated(motifID)], var_cutoff[j])) #keep only 400 unique variability in calculating the median
  for (i in 1:8){
    together2$effect="Null"
    together2$effect[which(together2$acf_exVall_opt>=corr_cutoff[i] & together2$chromVar_variability>together2$cutoff & together2$chromVar_p_value_adj < 0.01)]="Pos"
    together2$effect[which(together2$acf_exVall_opt<=(-corr_cutoff[i]) & together2$chromVar_variability>together2$cutoff & together2$chromVar_p_value_adj < 0.01)]="Neg"
    together2k=spread(together2[,c(4:6,21,20)],key=4, value=5)
    together2k <- together2k %>% group_by(motifID) %>% dplyr::slice_max(order_by = rowSums(across(`ATAC-signal_Both`:`TSS-signal_Both`) != "Null"), n=1, with_ties = FALSE)
    together2k=left_join(together2k,motif.df8b, by=c("motifID"="V1"),copy=F)
    together2k$`ATAC+`="No"
    together2k$`ATAC+`[which(together2k$`ATAC-signal_Both`=="Pos" )]="Yes"
    together2k$`ATAC-`="No"
    together2k$`ATAC-`[which(together2k$`ATAC-signal_Both`=="Neg" )]="Yes"
    together2k$`TSS+`="No"
    together2k$`TSS+`[which(together2k$`TSS-signal_Both`=="Pos" )]="Yes"
    together2k$`TSS-`="No"
    together2k$`TSS-`[which(together2k$`TSS-signal_Both`=="Neg" )]="Yes"
    outputfile2=together2k[which(rowSums(together2k[,c(11:14)] == "No")<4),] 
    outputfile2=outputfile2[which(rowSums(outputfile2[,c(7:10)] == "No")<4),] 
    FE1=outputfile2[,c("motifID","ATAC+","ATAC-","TSS+","TSS-")]
    FE1[is.na(FE1)]="No"
    FEsum=reshape2::melt(FE1,id=1)
    FEsum2=reshape2::melt(outputfile2[,c("motifID","TF_activator","TF_repressor","TF_pioneer","TF_closer")], id=1)
    FEsum=left_join(FEsum,FEsum2, by="motifID", copy=F, relationship = "many-to-many")
    FEsum2=FEsum%>%group_by(variable.x,variable.y,value.x, value.y)%>%dplyr::summarise(count=n())
    FEsum2$key=paste0(FEsum2$value.x,"_",FEsum2$value.y)
    FEsum3=spread(FEsum2[,c(1,2,5,6)], key=4, value=3)
    FEsum3[is.na(FEsum3)]=0
    for(k in 1:nrow(FEsum3)){
      GSEATasting <- matrix(c(FEsum3$No_No[k], FEsum3$No_Yes[k], FEsum3$Yes_No[k], FEsum3$Yes_Yes[k]), nrow = 2, dimnames = list(oligo1 = c("no", "yes"), oligo2 = c("no", "yes")))
      FEsum3$p.val[k] = fisher.test(GSEATasting, alternative = "two.sided")$p.value
      FEsum3$OR[k] = fisher.test(GSEATasting, alternative = "two.sided")$estimate}
    FEsum3$var_cutoff=var_cutoff[j]
    FEsum3$corr_cutoff=corr_cutoff[i]
    data7=rbind(data7, FEsum3)
  }}

data7$logOR=log(data7$OR)
data7$label="***"
data7$label[which(data7$p.val>=0.001)]="**"
data7$label[which(data7$p.val>=0.01)]="*"
data7$label[which(data7$p.val>=0.05)]="n.s."
data7$label1=paste0("log(OR)=",signif(data7$logOR,2),data7$label)

data7$var_cutoff=paste0("Var > P",data7$var_cutoff*100)
data7$var_cutoff=gsub("P50","Median",data7$var_cutoff)
data7$var_cutoff=factor(data7$var_cutoff, levels=c("Var > P25","Var > Median","Var > P75"))

write.table(data7,gzfile(paste0(path_fig3_data,"TF_grouping_by_joint_pen_chromvar_correlation_and_literature.FE.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
# -> for fig ex3h

#===============================================================================


#===============================================================================
# enhancer or promoter specific, compared to pool chromvar?
together2=read.delim(paste0(path_fig3_data,"auto_correlation_to_p_en.z.summarise_noFP.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
together2$effect="Null"
together2$effect[which(together2$acf_exVall_opt>=0.5 & together2$chromVar_variability>together2$median_var & together2$chromVar_p_value_adj < 0.01)]="Pos"
together2$effect[which(together2$acf_exVall_opt<=(-0.5) & together2$chromVar_variability>together2$median_var & together2$chromVar_p_value_adj < 0.01)]="Neg"
together2$group3=paste0(together2$group2,"_",together2$group1)
together2k=spread(together2[,c(4:6,21,20)],key=4, value=5)
together2k$rank=rowSums(together2k[,c(6,7,10,11)] != "Null")
together2k <- together2k %>% group_by(motifID) %>% dplyr::slice_max(rank, with_ties = FALSE)

en_ATAC_spec=together2k$motifID[which(together2k$`ATAC-signal_Enhancer`!="Null" & together2k$`ATAC-signal_Both` == "Null")]
en_TSS_spec=together2k$motifID[which(together2k$`TSS-signal_Enhancer`!="Null" & together2k$`TSS-signal_Both` == "Null")]
en_specific=union(en_ATAC_spec,en_TSS_spec)

p_ATAC_spec=together2k$motifID[which(together2k$`ATAC-signal_Promoter`!="Null" & together2k$`ATAC-signal_Both` == "Null")]
p_TSS_spec=together2k$motifID[which(together2k$`TSS-signal_Promoter`!="Null" & together2k$`TSS-signal_Both` == "Null")]
p_specific=union(p_ATAC_spec,p_TSS_spec)

en_ATAC_same=together2k$motifID[which(together2k$`ATAC-signal_Enhancer`!="Null" & together2k$`ATAC-signal_Both` != "Null")]
en_TSS_same=together2k$motifID[which(together2k$`TSS-signal_Enhancer`!="Null" & together2k$`TSS-signal_Both` != "Null")]
en_same=union(en_ATAC_same,en_TSS_same)

p_ATAC_same=together2k$motifID[which(together2k$`ATAC-signal_Promoter`!="Null" & together2k$`ATAC-signal_Both` != "Null")]
p_TSS_same=together2k$motifID[which(together2k$`TSS-signal_Promoter`!="Null" & together2k$`TSS-signal_Both` != "Null")]
p_same=union(p_ATAC_same,p_TSS_same)

length(en_ATAC_spec) #43
length(en_TSS_spec) #60
length(en_specific) #91
length(p_ATAC_spec) #42
length(p_TSS_spec) #29
length(p_specific) #68

#-> enhancer TSS discovered more motif when enhancer is considered separately




