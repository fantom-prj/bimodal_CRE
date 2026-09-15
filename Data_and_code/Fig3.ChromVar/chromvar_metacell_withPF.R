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
chromvar_folder=paste0(primary_folder,"Data_and_code/ChromVar/")
metacell=paste0(primary_folder,"Data_and_code/metacell/")
bigDATA=paste0(primary_folder,"Data_and_code/DATA/")
path_fig3_data=paste0(primary_folder,"Fig3/data/")


#===============================================================================
#===============================================================================
# run ChromVar considering footprint
path22=paste0(chromvar_folder,"output_withFP/")

pfm0 <- readJASPARMatrix(paste0(chromvar_folder,"jaspar2024/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.adjust.txt"), matrixClass = "PFM")
combined=readRDS("/analysisdata/fantom6/Interactome/single_cell_wallace/post_meta/metacell.combined.feature.Rds")
motif.df=read.delim(paste0(chromvar_folder,"jaspar2024/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.ID.tsv"), header=F, stringsAsFactors = F)

scACREv3=read.delim(paste0(primary_folder,"Data_and_code/scATAC/merged_peak/all_aCRE.PE.final.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
scACREv3a=scACREv3[which(scACREv3$tCRE == "Yes"),]
scACREv3a=scACREv3a[-which(scACREv3a$promoter_type_CT == "promoter-like" & is.na(scACREv3a$gene)),]
scACREv3a$peakID=gsub("_","-",scACREv3a$peakID)

footprint_mat <- read.delim(paste0(primary_folder,"Data_and_code/Foot_Printing/output/Footprint_0.2_cutoff_400_CRE463966.matrix.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
rownames(footprint_mat) <- footprint_mat$peakID
footprint_mat=footprint_mat[,-c(1:5)]
colnames(footprint_mat)=sapply(strsplit(colnames(footprint_mat),"_"),"[",2)

combined.p=combined
need.atac.p=which(rownames(combined.p[["aCRE"]]) %in% scACREv3a$peakID[which(scACREv3a$promoter_type_CT == "promoter-like")])
need.rna.p=which(rownames(combined.p[["tCRE"]]) %in% scACREv3a$peakID[which(scACREv3a$promoter_type_CT == "promoter-like")])
combined.p[["aCRE"]] <- subset(combined.p[["aCRE"]], features = rownames(combined.p[["aCRE"]])[need.atac.p])
combined.p[["tCRE"]] <- subset(combined.p[["tCRE"]], features = rownames(combined.p[["tCRE"]])[need.rna.p])
non_zero_frag_peaks <- which(rowSums(combined.p[["tCRE"]]) != 0) #remove peak that is zero across all cells
combined.p[["tCRE"]] <- subset(combined.p[["tCRE"]], features = rownames(combined.p[["tCRE"]])[non_zero_frag_peaks])

# tCRE promoter
peak_data <- combined.p[["tCRE"]]@counts
peak_ranges <- combined.p[["tCRE"]]@ranges
rse <- SummarizedExperiment(assays = list(peak_data), rowRanges = peak_ranges)
rse <- addGCBias(rse, genome = BSgenome.Hsapiens.UCSC.hg38)
names(rse@assays) <- c("counts")
motif_ix <- matchMotifs(pfm0, rse, genome = BSgenome.Hsapiens.UCSC.hg38, p.cutoff = 1e-04)

#=== adjust to footprint (equal to footprint + motif)
peak_order <- rownames(rse) #filter to only promoter tCRE
footprint_mat_aligned <- footprint_mat[peak_order, ]
footprint_logical <- (footprint_mat_aligned == "Yes") # Convert to logical if it's currently "Yes"/"No" or 1/0
common_motifs <- intersect(colnames(motif_ix), colnames(footprint_logical)) #filter to 400 motif
motif_ix_filtered <- motif_ix[, common_motifs]
footprint_subset <- footprint_logical[rownames(rse), common_motifs]
assay(motif_ix_filtered) <- as.matrix(footprint_subset) # Overwrite the assay

bg <- getBackgroundPeaks(object = rse)
dev <- computeDeviations(object = rse, annotations = motif_ix_filtered, background_peaks = bg)
variability <- computeVariability(dev)
plotVariability(variability, use_plotly = FALSE) 
variability=left_join(variability, motif.df, by=c("name"="V1"))
dev_df=data.frame(t(dev@assays@data@listData[["deviations"]]))
z_df=data.frame(t(dev@assays@data@listData[["z"]]))
write.table(dev_df,gzfile(paste0(path22,"p_tCRE_dev.tsv.gz")), col.names=T, row.names=T, sep="\t", quote=F)
write.table(z_df,gzfile(paste0(path22,"p_tCRE_z.tsv.gz")), col.names=T, row.names=T, sep="\t", quote=F)
write.table(variability,gzfile(paste0(path22,"p_tCRE_var.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

# aCRE promoter
peak_data <- combined.p[["aCRE"]]@counts
peak_ranges <- combined.p[["aCRE"]]@ranges
rse1 <- SummarizedExperiment(assays = list(peak_data), rowRanges = peak_ranges)
rse1 <- addGCBias(rse1, genome = BSgenome.Hsapiens.UCSC.hg38)
names(rse1@assays) <- c("counts")
motif_ix1 <- matchMotifs(pfm0, rse1, genome = BSgenome.Hsapiens.UCSC.hg38, p.cutoff = 1e-04)

#=== adjust to footprint (equal to footprint + motif)
peak_order1 <- rownames(rse1) #filter to only promoter aCRE
footprint_mat_aligned1 <- footprint_mat[peak_order1, ]
footprint_logical1 <- (footprint_mat_aligned1 == "Yes") # Convert to logical if it's currently "Yes"/"No" or 1/0
common_motifs1 <- intersect(colnames(motif_ix1), colnames(footprint_logical1)) #filter to 400 motif
motif_ix_filtered1 <- motif_ix1[, common_motifs1]
footprint_subset1 <- footprint_logical1[rownames(rse1), common_motifs1]
assay(motif_ix_filtered1) <- as.matrix(footprint_subset1) # Overwrite the assay

bg1 <- getBackgroundPeaks(object = rse1)
dev1 <- computeDeviations(object = rse1, annotations = motif_ix_filtered1, background_peaks = bg1)
variability1 <- computeVariability(dev1)
plotVariability(variability1, use_plotly = FALSE) 
variability1=left_join(variability1, motif.df, by=c("name"="V1"))
dev_df1=data.frame(t(dev1@assays@data@listData[["deviations"]]))
z_df1=data.frame(t(dev1@assays@data@listData[["z"]]))
write.table(dev_df1,gzfile(paste0(path22,"p_aCRE_dev.tsv.gz")), col.names=T, row.names=T, sep="\t", quote=F)
write.table(z_df1,gzfile(paste0(path22,"p_aCRE_z.tsv.gz")), col.names=T, row.names=T, sep="\t", quote=F)
write.table(variability1,gzfile(paste0(path22,"p_aCRE_var.tsv.gz")), col.names=T, row.names=T, sep="\t", quote=F)

#===============================================================================
# enhancer
combined.en=combined
need.atac.en=which(rownames(combined.en[["aCRE"]]) %in% scACREv3a$peakID[which(scACREv3a$promoter_type_CT == "enhancer-like")])
need.rna.en=which(rownames(combined.en[["tCRE"]]) %in% scACREv3a$peakID[which(scACREv3a$promoter_type_CT == "enhancer-like")])
combined.en[["aCRE"]] <- subset(combined.en[["aCRE"]], features = rownames(combined.en[["aCRE"]])[need.atac.en])
combined.en[["tCRE"]] <- subset(combined.en[["tCRE"]], features = rownames(combined.en[["tCRE"]])[need.rna.en])
#remove peak that is zero across all cells
non_zero_frag_peaks <- which(rowSums(combined.en[["tCRE"]]) != 0)
combined.en[["tCRE"]] <- subset(combined.en[["tCRE"]], features = rownames(combined.en[["tCRE"]])[non_zero_frag_peaks])

peak_data <- combined.en[["tCRE"]]@counts
peak_ranges <- combined.en[["tCRE"]]@ranges
rse <- SummarizedExperiment(assays = list(peak_data),rowRanges = peak_ranges)
rse <- addGCBias(rse, genome = BSgenome.Hsapiens.UCSC.hg38)
names(rse@assays) <- c("counts")
motif_ix <- matchMotifs(pfm0, rse, genome = BSgenome.Hsapiens.UCSC.hg38, p.cutoff = 1e-04)

#=== adjust to footprint (equal to footprint + motif)
peak_order <- rownames(rse) #filter to only enhancer tCRE
footprint_mat_aligned <- footprint_mat[peak_order, ]
footprint_logical <- (footprint_mat_aligned == "Yes") # Convert to logical if it's currently "Yes"/"No" or 1/0
common_motifs <- intersect(colnames(motif_ix), colnames(footprint_logical)) #filter to 400 motif
motif_ix_filtered <- motif_ix[, common_motifs]
footprint_subset <- footprint_logical[rownames(rse), common_motifs]
assay(motif_ix_filtered) <- as.matrix(footprint_subset) # Overwrite the assay

bg <- getBackgroundPeaks(object = rse)
dev <- computeDeviations(object = rse, annotations = motif_ix_filtered, background_peaks = bg)
variability <- computeVariability(dev)
plotVariability(variability, use_plotly = FALSE) 
variability=left_join(variability, motif.df, by=c("name"="V1"))
dev_df=data.frame(t(dev@assays@data@listData[["deviations"]]))
z_df=data.frame(t(dev@assays@data@listData[["z"]]))
write.table(dev_df,gzfile(paste0(path22,"en_tCRE_dev.tsv.gz")), col.names=T, row.names=T, sep="\t", quote=F)
write.table(z_df,gzfile(paste0(path22,"en_tCRE_z.tsv.gz")), col.names=T, row.names=T, sep="\t", quote=F)
write.table(variability,gzfile(paste0(path22,"en_tCRE_var.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#aCRE
peak_data <- combined.en[["aCRE"]]@counts
peak_ranges <- combined.en[["aCRE"]]@ranges
rse1 <- SummarizedExperiment(assays = list(peak_data),rowRanges = peak_ranges)
rse1 <- addGCBias(rse1, genome = BSgenome.Hsapiens.UCSC.hg38)
names(rse1@assays) <- c("counts")
motif_ix1 <- matchMotifs(pfm0, rse1, genome = BSgenome.Hsapiens.UCSC.hg38, p.cutoff = 1e-04)

#=== adjust to footprint (equal to footprint + motif)
peak_order1 <- rownames(rse1) #filter to only enhancer aCRE
footprint_mat_aligned1 <- footprint_mat[peak_order1, ]
footprint_logical1 <- (footprint_mat_aligned1 == "Yes") # Convert to logical if it's currently "Yes"/"No" or 1/0
common_motifs1 <- intersect(colnames(motif_ix1), colnames(footprint_logical1)) #filter to 400 motif
motif_ix_filtered1 <- motif_ix1[, common_motifs1]
footprint_subset1 <- footprint_logical1[rownames(rse1), common_motifs1]
assay(motif_ix_filtered1) <- as.matrix(footprint_subset1) # Overwrite the assay

bg1 <- getBackgroundPeaks(object = rse1)
dev1 <- computeDeviations(object = rse1, annotations = motif_ix_filtered1, background_peaks = bg1)
variability1 <- computeVariability(dev1)
plotVariability(variability1, use_plotly = FALSE) 
variability1=left_join(variability1, motif.df, by=c("name"="V1"))
dev_df1=data.frame(t(dev1@assays@data@listData[["deviations"]]))
z_df1=data.frame(t(dev1@assays@data@listData[["z"]]))
write.table(dev_df1,gzfile(paste0(path22,"en_aCRE_dev.tsv.gz")), col.names=T, row.names=T, sep="\t", quote=F)
write.table(z_df1,gzfile(paste0(path22,"en_aCRE_z.tsv.gz")), col.names=T, row.names=T, sep="\t", quote=F)
write.table(variability1,gzfile(paste0(path22,"en_aCRE_var.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
#result stored in paste0(chromvar_folder,"output_withFP/")

#===============================================================================



#===============================================================================
# parsing chromVar results (with footprinting)
# gather variability of the 4 runs

path22=paste0(chromvar_folder,"output_withFP/")

files=list.files(path=path22, pattern="var")
files.names=sapply(strsplit(files,"_var"),"[", 1)

data=data.frame()
for (i in 1:length(files)){
  data1=read.delim(paste0(path22, files[i]), header=T, stringsAsFactors = F, check.names = F)
  data1$set=files.names[i]
  data=rbind(data, data1)}
data$set1=sapply(strsplit(data$set,"_"),"[",1)
data$set2=sapply(strsplit(data$set,"_"),"[",2)

data$set1=gsub("en","Enhancer",data$set1)
data$set1=gsub("p","Promoter",data$set1)
data$set2=gsub("aCRE","ATAC-signal",data$set2)
data$set2=gsub("tCRE","TSS-signal",data$set2)

data$label=paste0(data$name,"_",data$V2)
data2 <- data %>% group_by(set) %>% arrange(desc(variability), .by_group = TRUE) %>% mutate(rank = row_number()) %>% ungroup()

data2=data2%>%group_by(set)%>%dplyr::mutate(scaled_var=variability/max(variability))
data2=data2%>%group_by(set2)%>%dplyr::mutate(zscaled_var=(variability-mean(variability))/sd(variability))
data2=data2%>%group_by(set)%>%dplyr::mutate(median_var=median(variability), median_scaled_var=median(scaled_var), median_zscaled_var=median(zscaled_var))

write.table(data2, gzfile(paste0(path_fig3_data,"variability_ChromVar_withFP.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
# for figure ex3b
data2=read.delim(paste0(path_fig3_data,"variability_ChromVar_withFP.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data2$set1=factor(data2$set1, levels=c("Enhancer","Promoter"))
data2$V2[which(data2$rank >5)]=NA
data3=data2%>%group_by(set1, set2)%>%dplyr::summarise(median=median(variability),count=n())

ex3b=ggplot(data2[which(data2$set1 %in% c("Enhancer","Promoter")),], aes(x = rank, y = variability)) +
  facet_grid(cols=vars(set1), rows=vars(set2), scales = "free") +
  coord_cartesian(xlim=c(1,879)) +
  geom_errorbar(aes(ymin = bootstrap_lower_bound, ymax = bootstrap_upper_bound), linewidth = 0.2, color = "grey25") +
  geom_point(shape=20, size = 0.05, color = "black") +
  geom_hline(data=data3[which(data3$set1 %in% c("Enhancer","Promoter")),], aes(yintercept = median), linewidth=0.2, color="navy", linetype="dashed") +
  geom_text_repel(aes(label=V2), direction = "y", segment.size = 0.2 , nudge_x=250, hjust = 0, size=1.6, force=0.1) +
  labs(x = "Sorted 879 TF motifs", y = "Variability", title="Variability of TF motif activity") +
  theme1+ theme(axis.text.x = element_blank(), axis.ticks = element_blank())
pdf(paste0(path_fig3,"ex3b.variability_ChromVar_withFP.pdf"), width = 2.6, height = 2.3)
print (ex3b)
dev.off() 

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
path22=paste0(chromvar_folder,"output_withFP/")

aCRE.en=read.delim(paste0(path22,"en_aCRE_z.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
tCRE.en=read.delim(paste0(path22,"en_tCRE_z.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
aCRE.p=read.delim(paste0(path22,"p_aCRE_z.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
tCRE.p=read.delim(paste0(path22,"p_tCRE_z.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)

aCRE.en$metacell=rownames(aCRE.en)
aCRE.en=reshape2::melt(aCRE.en, id=401)
aCRE.en=aCRE.en%>%group_by(variable)%>%dplyr::mutate(variance=var(value))
aCRE.en$group1="ATAC-signal"
aCRE.en$group2="Enhancer"
tCRE.en$metacell=rownames(tCRE.en)
tCRE.en=reshape2::melt(tCRE.en, id=401)
tCRE.en=tCRE.en%>%group_by(variable)%>%dplyr::mutate(variance=var(value))
tCRE.en$group1="TSS-signal"
tCRE.en$group2="Enhancer"

aCRE.p$metacell=rownames(aCRE.p)
aCRE.p=reshape2::melt(aCRE.p, id=401)
aCRE.p=aCRE.p%>%group_by(variable)%>%dplyr::mutate(variance=var(value))
aCRE.p$group1="ATAC-signal"
aCRE.p$group2="Promoter"
tCRE.p$metacell=rownames(tCRE.p)
tCRE.p=reshape2::melt(tCRE.p, id=401)
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

write.table(overall, gzfile(paste0(path_fig3_data,"plot.table.all_withFP.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#==============
#trajectory 2 version
overall=rbind(aCRE.en,aCRE.p,tCRE.en,tCRE.p)
overall=overall[which(overall$metacell %in% pseudotime3$SEACell),]
overall=left_join(overall, pseudotime3[,c(1,4,6)], by=c("metacell"="SEACell"),copy=F)
overall=left_join(overall, motif.df6[,c(1,2,7,8,9)], by=c("metacell"="metacell", "variable"="motifID"), copy=F)
colnames(overall)[c(10,11)]=c("expression_TF1","expression_TF2")

write.table(overall, gzfile(paste0(path_fig3_data,"plot.table.all.trajectory2_withFP.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)


# dot plot #z-score version, match better with the variability result
#overall=read.delim("/analysisdata/fantom6/Interactome/single_cell_wallace/post_meta/plot.table.all.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
overall=read.delim(paste0(path_fig3_data,"plot.table.all.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
overall2=overall%>%group_by(group1, group2, variable,motifName,T2)%>%dplyr::summarise(value2=mean(value), TF1=mean(expression_TF1), TF2=mean(expression_TF2)) 
overall2=overall2 %>%group_by(variable)  %>%dplyr::mutate (rel.TF1=TF1*0.5/max(TF1), rel.TF2=TF2*0.5/max(TF2))

selection=c("Pou5f1::Sox2","TEAD4","FOS::JUND","NEUROD1","ONECUT2","ZEB1","MEIS3","Sox3")
overall1=overall2[which(overall2$motifName %in% selection),]
overall1$motifName=factor(overall1$motifName, levels=c("Pou5f1::Sox2","TEAD4","FOS::JUND","NEUROD1","ONECUT2","INSM1","ZEB1","MEIS3","Sox3"))
f3c= ggplot(overall1) +
  geom_point(mapping=aes(x =T2, y = value2, color=group2), size=0.2, shape=19) + 
  stat_smooth(mapping=aes(x =T2, y = rel.TF1*20), method = "lm", formula = y ~ poly(x, 21), se = FALSE, color="black", alpha=1, size=0.2) +
  stat_smooth(mapping=aes(x =T2, y = rel.TF2*20), method = "lm", formula = y ~ poly(x, 21), se = FALSE, color="grey", alpha=0.6, size=0.2) +
  scale_color_manual(values=c("#7E6148FF","#B09C85FF"))+
  facet_grid(row=vars(group1), col=vars(motifName), scales = "free_x")+
  labs(color=NULL, y="Relative TF motif activity & relative TF expression", x = "Ranked by pseudotime for trajectory 1", title = "TF motif activity by chromatin accessibility and transcription") +
  theme1+theme(axis.text.x = element_blank())
pdf(paste0(path_fig3,"f3c.Chromvar.Pseudotime.z.neuronal.trajectory.aCRE.vs.tCRE.expression.pdf"), width = 6, height = 2.5)
f3c
dev.off()

#trajectory2
boverall=read.delim(paste0(path_fig3_data,"plot.table.all.trajectory2.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
boverall2=boverall%>%group_by(group1, group2, variable,motifName,T2)%>%dplyr::summarise(value2=mean(value), TF1=mean(expression_TF1), TF2=mean(expression_TF2)) 
boverall2=boverall2 %>%group_by(variable)  %>%dplyr::mutate (rel.TF1=TF1*0.5/max(TF1), rel.TF2=TF2*0.5/max(TF2))

selection=c("TCF7L2","Hand1","FOS::JUND*","NEUROD1*","BCL11A","CREB3L1*","REST*","GATA6","OSR1*")
boverall1=boverall2[which(boverall2$motifName %in% selection),]
boverall1$motifName=factor(boverall1$motifName, levels=c("TCF7L2","Hand1","BCL11A","GATA6"))

xi= ggplot(boverall1) +
  #scale_y_continuous(limits=c(-0.4, 1.05), breaks=c(0,0.5,1))+
  geom_point(mapping=aes(x =T2, y = value2, color=group2), size=0.2, shape=19) + 
  stat_smooth(mapping=aes(x =T2, y = rel.TF1*20), method = "lm", formula = y ~ poly(x, 11), se = FALSE, color="black", alpha=1, size=0.2) +
  #stat_smooth(mapping=aes(x =T2, y = rel.TF2*20), method = "lm", formula = y ~ poly(x, 11), se = FALSE, color="grey", alpha=0.6, size=0.2) +
  scale_color_manual(values=c("#7E6148FF","#B09C85FF"))+
  scale_x_continuous(limits=c(43,84))+
  facet_grid(row=vars(group1), col=vars(motifName), scales = "free_x")+
  labs(color=NULL, y="Mean TF motif activity & TF relative expression", x = "Ranked by pseudotime for trajectory 2 (excluding iPSC)", title = "TF motif activity by chromatin accessibility and transcription") +
  theme1+theme(axis.text.x = element_blank())
pdf(paste0(path_fig3,"f3c2.Chromvar.Pseudotime.z.OPC.trajectory.aCRE.vs.tCRE.expression.pdf"), width = 3.5, height = 2.5)
xi
dev.off()



#===============================================================================
# perform correlation
path22=paste0(chromvar_folder,"output_withFP/")
tCRE.en_v=read.delim(paste0(path22,"en_tCRE_var.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
aCRE.en_v=read.delim(paste0(path22,"en_aCRE_var.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
tCRE.en_z=read.delim(paste0(path22,"en_tCRE_z.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
aCRE.en_z=read.delim(paste0(path22,"en_aCRE_z.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

tCRE.p_v=read.delim(paste0(path22,"p_tCRE_var.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
aCRE.p_v=read.delim(paste0(path22,"p_aCRE_var.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
tCRE.p_z=read.delim(paste0(path22,"p_tCRE_z.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
aCRE.p_z=read.delim(paste0(path22,"p_aCRE_z.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

tCRE.en_z$metacell=rownames(tCRE.en_z)
tCRE.en_z=reshape2::melt(tCRE.en_z, id=401)
tCRE.en_z=left_join(tCRE.en_z, tCRE.en_v[,c(1,2,6,7)], by=c("variable"="name"),copy=F)
tCRE.en_z$group1="TSS"
tCRE.en_z$group2="Enhancer"
aCRE.en_z$metacell=rownames(aCRE.en_z)
aCRE.en_z=reshape2::melt(aCRE.en_z, id=401)
aCRE.en_z=left_join(aCRE.en_z, aCRE.en_v[,c(1,2,6,7)], by=c("variable"="name"),copy=F)
aCRE.en_z$group1="ATAC"
aCRE.en_z$group2="Enhancer"
tCRE.p_z$metacell=rownames(tCRE.p_z)
tCRE.p_z=reshape2::melt(tCRE.p_z, id=401)
tCRE.p_z=left_join(tCRE.p_z, tCRE.p_v[,c(1,2,6,7)], by=c("variable"="name"),copy=F)
tCRE.p_z$group1="TSS"
tCRE.p_z$group2="Promoter"
aCRE.p_z$metacell=rownames(aCRE.p_z)
aCRE.p_z=reshape2::melt(aCRE.p_z, id=401)
aCRE.p_z=left_join(aCRE.p_z, aCRE.p_v[,c(1,2,6,7)], by=c("variable"="name"),copy=F)
aCRE.p_z$group1="ATAC"
aCRE.p_z$group2="Promoter"
overall.all=rbind(tCRE.en_z,aCRE.en_z,tCRE.p_z,aCRE.p_z)

overall.all=left_join(overall.all, pseudotime2[,c(1,4,6)], by=c("metacell"="SEACell"),copy=F)
overall.all=left_join(overall.all, motif.df6[,c(1,7:9)], by=c("metacell"="metacell", "variable"="motifID"), copy=F)
colnames(overall.all)[c(2,3,6,11,12)]=c("motifID","chromvar_z","motifName","expression_TF1","expression_TF2")

overall.all2=overall.all%>%group_by(group1, group2, motifID,motifName,variability,p_value_adj,T2)%>%dplyr::summarise(chromvar_z=mean(chromvar_z), TF1=mean(expression_TF1), TF2=mean(expression_TF2)) %>%dplyr::mutate (rel.TF1=TF1*max(chromvar_z)/max(TF1), rel.TF2=TF2*max(chromvar_z)/max(TF2))

length(unique(overall.all2$motifID[which(overall.all2$p_value_adj < 0.01 & overall.all2$variability > 2)])) #285

#============================
#split TF1 and TF2 for correlation
overall.all3c=reshape2::melt(overall.all2[,c(1:10)], id=c(1:8))
colnames(overall.all3c)[c(9,10)]=c("TF","TF_expression")

overall.all3d=overall.all3c%>%group_by(group1,group2,motifID,motifName,variability,p_value_adj,TF)%>%dplyr::summarise(max_TF_exp=max(TF_expression))
overall.all3d$label=paste0(overall.all3d$motifID,"_",overall.all3d$motifName,"_",overall.all3d$TF)
overall.all3d$motif=paste0(overall.all3d$motifID,"_",overall.all3d$motifName)
overall.all3d=overall.all3d[which(!is.na(overall.all3d$max_TF_exp)),]
overall.all3d$label[-grep("\\:",overall.all3d$motif)]=overall.all3d$motif[-grep("\\:",overall.all3d$motif)]

overall.all3d$group3=paste0(overall.all3d$group1,"_",overall.all3d$group2)

overall.all3e=spread(overall.all3d[,c(11,9,10,8,5)],key=1, value=5)
overall.all3e1=spread(overall.all3d[,c(11,9,10,8,6)],key=1, value=5)
overall.all3e=left_join(overall.all3e,overall.all3e1[,c(1,4:7)], by="label", copy=F, suffix=c("_var","_padj"))

overall.all3e$var_ATAC_enhancer="no"
overall.all3e$var_ATAC_enhancer[which(overall.all3e$ATAC_Enhancer_padj <0.01 & overall.all3e$ATAC_Enhancer_var>median(overall.all3e$ATAC_Enhancer_var))]="yes"
overall.all3e$var_TSS_enhancer="no"
overall.all3e$var_TSS_enhancer[which(overall.all3e$TSS_Enhancer_padj <0.01 & overall.all3e$TSS_Enhancer_var>median(overall.all3e$TSS_Enhancer_var))]="yes"
overall.all3e$var_ATAC_promoter="no"
overall.all3e$var_ATAC_promoter[which(overall.all3e$ATAC_Promoter_padj <0.01 & overall.all3e$ATAC_Promoter_var>median(overall.all3e$ATAC_Promoter_var))]="yes"
overall.all3e$var_TSS_promoter="no"
overall.all3e$var_TSS_promoter[which(overall.all3e$TSS_Promoter_padj <0.01 & overall.all3e$TSS_Promoter_var>median(overall.all3e$TSS_Promoter_var))]="yes"

#-> all yes according to the cut off max>=0.2 (per meta cell), >=3 meta cells, require both TF

overall.all3e=data.frame(overall.all3e)
row.names(overall.all3e)=overall.all3e$label
overall.all3f=overall.all3e[which(overall.all3e$var_ATAC_enhancer=="yes" | overall.all3e$var_TSS_enhancer == "yes"|overall.all3e$var_ATAC_promoter=="yes" | overall.all3e$var_TSS_promoter == "yes"),] 
length(unique(overall.all3f$motif)) #342 motif pass ChromVar variation

#====================
#get the correlation for all the expressed TF first, filter later
length(unique(overall.all3c$motifID)) #400
overall.all3c$motif=paste0(overall.all3c$motifID,"_",overall.all3c$motifName)
overall.all3c$label=paste0(overall.all3c$motif,"_",overall.all3c$TF)
overall.all3c$label[-grep("\\:",overall.all3c$motif)]=overall.all3c$motif[-grep("\\:",overall.all3c$motif)]
length(unique(overall.all3c$motifID)) #400
write.table(overall.all3c,gzfile(paste0(path_fig3_data,"chromvar_z_expression_400motif_withFP.tsv.gz")), col.names = T, row.names = F, sep="\t", quote=F)

#exclude NA from expression level when calculating correlation
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
together2=together2%>%group_by(group1,group2)%>%dplyr::mutate(chromVar_var_median=median(chromVar_variability))
together2$chromVar_Var="no"
together2$chromVar_Var[which(together2$chromVar_p_value_adj <0.01 & together2$chromVar_variability>together2$chromVar_var_median)]="yes"


together2$group1[which(together2$group1=="ATAC")]="ATAC-signal"
together2$group1[which(together2$group1=="TSS")]="TSS-signal"

var1=read.delim(paste0(path_fig3_data,"variability_ChromVar_withFP.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
together2=left_join(together2, var1[,c("name","set1","set2", "p_value_adj","zscaled_var","median_zscaled_var")], by=c("group1"="set2","group2"="set1","motifID"="name"),copy=F)

together2$pass_both = "No"
together2$pass_both[which(abs(together2$acf_exVall_opt) > 0.5 & together2$chromVar_Var == "yes")] = "Yes"
write.table(together2,gzfile(paste0(path_fig3_data,"auto_correlation_to_p_en.z.summarise_withFP.tsv.gz")), col.names = T, row.names = F, sep="\t", quote=F)
together2=read.delim(paste0(path_fig3_data,"auto_correlation_to_p_en.z.summarise_withFP.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

together2a=unique(together2[,c("group1", "group2","median_zscaled_var")])
ex3c=ggplot() + 
  labs(x="ChromVar Variability (Scaled)", y ="Abs(Correlation)", alpha="Pass chromVar\n& Correlation", title="Evaluable TF mmotifs")+
  facet_grid(rows=vars(group1), cols=vars(group2))+
  geom_point(data=together2, mapping=aes(x=zscaled_var, color=interaction(group1, group2), y=abs(acf_exVall_opt),  alpha = pass_both), size=0.2, shape=19)+
  scale_color_npg(guide=NULL)+
  scale_alpha_manual(values=c("Yes"=1, "No"=0.3))+
  geom_hline(yintercept=0.5, linetype="dashed", linewidth=0.2)+
  geom_vline(data=together2a,aes(xintercept=median_zscaled_var), linetype="dashed", linewidth=0.2)+
  theme1+theme(legend.position = "right")
pdf(paste0(path_fig3,"ex3c.scale_variability_n_correlation_xy_withFP.pdf"), width = 2.6, height = 1.8)
print (ex3c)
dev.off() 



