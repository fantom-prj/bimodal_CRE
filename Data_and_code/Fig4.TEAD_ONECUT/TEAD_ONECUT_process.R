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
library(shadowtext)
library(ggrastr)

###color#####
library("ggsci")
mypal = pal_npg("nrc", alpha = 1)(9)
mypal
library("scales")
show_col(mypal)

#####################

#primary_folder=[primary_folder]
primary_folder="/analysisdata/fantom6/Interactome/single_cell_wallace/Figure/"
FP_output=paste0(primary_folder,"Data_and_code/Fig4.Foot_Printing/output/")
primary_FP_path="/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/subset_BAM/BINDetect_compare_output/"
ABC_folder=paste0(primary_folder,"Data_and_code/Fig4.ABC/")
ChIP_path=paste0(primary_folder,"Data_and_code/Fig4.ChIP_seq/")
path_fig4_data=paste0(primary_folder,"Fig4/data/")

setwd(paste0(primary_folder,"Data_and_code/Fig4.TEAD_ONECUT"))

#===============================================================================
# 1. perform tobias footprint using ATAC-seq data summarized into 10 cell clusters
# according to Data_and_code/Fig4.Foot_Printing/footprint_process.R
# running parameters and all the summary tables are stored in Data_and_code/Fig4.Foot_Printing

#===============================================================================
# 2. Run analyses after ChIP-seq of TEAD4 and ONECUT2
# according to Data_and_code/Fig4.ChIP_seq/chipseq.process.R
# running parameters and results are stored in Data_and_code/Fig4.ChIP_seq

#===============================================================================
# 3. Run analyses after KD of TEAD2/4 and ONECUT2 followed by bulk ATAC-seq
# according to Data_and_code/Fig4.KD_ATAC/KD_ATAC.R
# running parameters and results are stored in Data_and_code/Fig4.KD_ATAC

#===============================================================================
# 4. Run ABC-model
# according to Data_and_code/Fig4.ABC/ABC_input_and_output.R
# running parameters and results are stored in Data_and_code/Fig4.ABC

#===============================================================================

#chip seq correlate with footprinting?
library(pROC)

TEAD4final=read.delim(paste0(ChIP_path,"final_data/TEAD4_ChIP_aCRE_final.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
TEAD4final=TEAD4final[which(TEAD4final$tobias_motif == "yes"),]

tead_FPb1=read.delim(paste0(primary_FP_path,"TEAD4_MA0809.3/beds/TEAD4_MA0809.3_out_iPS_0_undifferentiaed_CB_footprints_bound.bed"), header=F, stringsAsFactors = F)
tead_FPub1=read.delim(paste0(primary_FP_path,"TEAD4_MA0809.3/beds/TEAD4_MA0809.3_out_iPS_0_undifferentiaed_CB_footprints_unbound.bed"), header=F, stringsAsFactors = F)

tead_FP1=rbind(tead_FPb1,tead_FPub1)
tead_FP1=tead_FP1%>%group_by(V10)%>%dplyr::slice_max(V13)
tead_FP1=unique(tead_FP1[,c(10,13)])
tead_FP1$ChIP="No"
tead_FP1$ChIP[which(tead_FP1$V10 %in% TEAD4final$peakID)]="Yes"

roc_data1 <- roc(tead_FP1$ChIP, tead_FP1$V13, levels=c("Yes","No"))
auc1=auc(roc_data1) #0.937
write.table(tead_FP1,gzfile(paste0(path_fig4_data,"tobias_chip_tead4_cluster1.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)
# -> for fig ex4c

#===
tead_FPb2=read.delim(paste0(primary_FP_path,"TEAD4_MA0809.3/beds/TEAD4_MA0809.3_out_iPS_1_neuroepithelial_CB_footprints_bound.bed"), header=F, stringsAsFactors = F)
tead_FPub2=read.delim(paste0(primary_FP_path,"TEAD4_MA0809.3/beds/TEAD4_MA0809.3_out_iPS_1_neuroepithelial_CB_footprints_unbound.bed"), header=F, stringsAsFactors = F)
tead_FP2=rbind(tead_FPb2,tead_FPub2)
tead_FP2=tead_FP2%>%group_by(V10)%>%dplyr::slice_max(V13)
tead_FP2=unique(tead_FP2[,c(10,13)])
tead_FP2$ChIP="No"
tead_FP2$ChIP[which(tead_FP2$V10 %in% TEAD4final$peakID)]="Yes"
roc_data2 <- roc(tead_FP2$ChIP, tead_FP2$V13, levels=c("Yes","No"))
auc2=auc(roc_data2) #0.938
write.table(tead_FP2,gzfile(paste0(path_fig4_data,"tobias_chip_tead4_cluster2.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

#=======
#ONECUT2
ONECUTfinal=read.delim(paste0(ChIP_path,"final_data/ONECUT_ChIP_aCRE_final.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
ONECUTfinal=ONECUTfinal[which(ONECUTfinal$tobias_motif == "yes"),]

onecut_FPb1=read.delim(paste0(primary_FP_path,"ONECUT2_MA0756.3/beds/ONECUT2_MA0756.3_out_neuron_immature_CB_footprints_bound.bed"), header=F, stringsAsFactors = F)
onecut_FPub1=read.delim(paste0(primary_FP_path,"ONECUT2_MA0756.3/beds/ONECUT2_MA0756.3_out_neuron_immature_CB_footprints_unbound.bed"), header=F, stringsAsFactors = F)
onecut_FP1=rbind(onecut_FPb1,onecut_FPub1)
onecut_FP1=onecut_FP1%>%group_by(V10)%>%dplyr::slice_max(V13)
onecut_FP1=unique(onecut_FP1[,c(10,13)])
onecut_FP1$ChIP="No"
onecut_FP1$ChIP[which(onecut_FP1$V10 %in% ONECUTfinal$peakID)]="Yes"
roc_data3 <- roc(onecut_FP1$ChIP, onecut_FP1$V13, levels=c("Yes","No"))
auc3=auc(roc_data3) #0.9302
write.table(onecut_FP1,gzfile(paste0(path_fig4_data,"tobias_chip_onecut2_cluster7.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

#===
onecut_FPb2=read.delim(paste0(primary_FP_path,"ONECUT2_MA0756.3/beds/ONECUT2_MA0756.3_out_neuron_mature_CB_footprints_bound.bed"), header=F, stringsAsFactors = F)
onecut_FPub2=read.delim(paste0(primary_FP_path,"ONECUT2_MA0756.3/beds/ONECUT2_MA0756.3_out_neuron_mature_CB_footprints_unbound.bed"), header=F, stringsAsFactors = F)
onecut_FP2=rbind(onecut_FPb2,onecut_FPub2)
onecut_FP2=onecut_FP2%>%group_by(V10)%>%dplyr::slice_max(V13)
onecut_FP2=unique(onecut_FP2[,c(10,13)])
onecut_FP2$ChIP="No"
onecut_FP2$ChIP[which(onecut_FP2$V10 %in% ONECUTfinal$peakID)]="Yes"
roc_data4 <- roc(onecut_FP2$ChIP, onecut_FP2$V13, levels=c("Yes","No"))
auc4=auc(roc_data4) #0.926
write.table(onecut_FP2,gzfile(paste0(path_fig4_data,"tobias_chip_onecut2_cluster8.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

#===
onecut_FPb3=read.delim(paste0(primary_FP_path,"ONECUT2_MA0756.3/beds/ONECUT2_MA0756.3_out_neuron_mature_2_CB_footprints_bound.bed"), header=F, stringsAsFactors = F)
onecut_FPub3=read.delim(paste0(primary_FP_path,"ONECUT2_MA0756.3/beds/ONECUT2_MA0756.3_out_neuron_mature_2_CB_footprints_unbound.bed"), header=F, stringsAsFactors = F)
onecut_FP3=rbind(onecut_FPb3,onecut_FPub3)
onecut_FP3=onecut_FP3%>%group_by(V10)%>%dplyr::slice_max(V13)
onecut_FP3=unique(onecut_FP3[,c(10,13)])
onecut_FP3$ChIP="No"
onecut_FP3$ChIP[which(onecut_FP3$V10 %in% ONECUTfinal$peakID)]="Yes"
roc_data5 <- roc(onecut_FP3$ChIP, onecut_FP3$V13, levels=c("Yes","No"))
auc5=auc(roc_data5) #0.912

write.table(onecut_FP3,gzfile(paste0(path_fig4_data,"tobias_chip_onecut2_cluster9.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)
# -> for fig ex4c

#===============================================================================
# prepare heatmap at metacell level
combined=readRDS(paste0(primary_folder,"Data_and_code/DATA/metacell.combined.feature.Rds"))
#subset the data only included aCRE that i needed
acre=union(TEAD4final0$peakID,ONECUTfinal0$peakID)
need.atac=which(rownames(combined[["aCRE"]]) %in% acre)
need.rna=which(rownames(combined[["tCRE"]]) %in% acre)
combined[["aCRE"]] <- subset(combined[["aCRE"]], features = rownames(combined[["aCRE"]])[need.atac])
combined[["tCRE"]] <- subset(combined[["tCRE"]], features = rownames(combined[["tCRE"]])[need.rna])

zt=as.data.frame(combined@assays[["tCRE"]]@data)
za=as.data.frame(combined@assays[["aCRE"]]@data)
pseudotime1=read.delim(paste0(primary_folder,"Data_and_code/metacell/output/pseudotime.scRNA.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
pseudotime1$cluster=sapply(strsplit(pseudotime1$SEACell,"_meta"),"[",1)
cluster_info=read.delim(paste0(primary_folder,"Fig1/cluster_info.tsv"),header=T, stringsAsFactors = F, check.names = F)
pseudotime1=left_join(pseudotime1,cluster_info[,c(1,3)],by="cluster",copy=F)
dd=data.frame(unique(pseudotime1$label))
dd$rank=c(1,2,3,4,6,5,0,7,8,9)
colnames(dd)[1]="label"
pseudotime1=left_join(pseudotime1,dd,by="label",copy=F)
pseudotime1=pseudotime1[which(pseudotime1$rank!=0),]
pseudotime1=pseudotime1[order(pseudotime1$rank,pseudotime1$pseudotime),]
pseudotime1=pseudotime1[which(pseudotime1$SEACell %in% colnames(za)),]
#
za=za[,pseudotime1$SEACell] #filter and get the ranking of pseudotime
zt=zt[,pseudotime1$SEACell] #filter and get the ranking of pseudotime
za$peakID=rownames(za)
zt$peakID=rownames(zt)

#TEAD_both with and without transcription
TEAD4final0=left_join(TEAD4final0[,-5], za, by="peakID", copy=F)
TEAD4final0=left_join(TEAD4final0, zt, by="peakID", copy=F, suffix=c("_kkATAC","_kkRNA"))
write.table(TEAD4final0,gzfile(paste0(path_fig4_data,"heatmap.aCRE.TEAD4.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
# -> for fig ex4a

#ONECUT_both with and without transcription
ONECUTfinal0=left_join(ONECUTfinal0[,-5], za, by="peakID", copy=F)
ONECUTfinal0=left_join(ONECUTfinal0, zt, by="peakID", copy=F, suffix=c("_kkATAC","_kkRNA"))
write.table(ONECUTfinal0,gzfile(paste0(path_fig4_data,"heatmap.aCRE.ONECUT.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
# -> for fig ex4a

#===============================================================================
# identify ChIP-seq peak with matched motif
# see if those ChIP-seq peak overlap with ATAC peak defined from each cell type
iPSC.ATAC=read.delim(paste0(primary_folder,"Data_and_code/scATAC/sample_ATAC_peak/peaks.iPS.bed.gz"), header=F, skip=51, stringsAsFactors = F, check.names = F)
NSC.ATAC=read.delim(paste0(primary_folder,"Data_and_code/scATAC/sample_ATAC_peak/peaks.NSC.bed.gz"), header=F, skip=51, stringsAsFactors = F, check.names = F)
NRN.ATAC=read.delim(paste0(primary_folder,"Data_and_code/scATAC/sample_ATAC_peak/peaks.NRN.bed.gz"), header=F, skip=51, stringsAsFactors = F, check.names = F)
library(GenomicRanges)
griPSC <- GRanges(seqnames = iPSC.ATAC$V1, ranges = IRanges(start = iPSC.ATAC$V2, end = iPSC.ATAC$V3))
grNSC <- GRanges(seqnames = NSC.ATAC$V1, ranges = IRanges(start = NSC.ATAC$V2, end = NSC.ATAC$V3))
grNRN <- GRanges(seqnames = NRN.ATAC$V1, ranges = IRanges(start = NRN.ATAC$V2, end = NRN.ATAC$V3))

TEAD4m=read.delim(paste0(primary_folder,"Data_and_code/Fig4.ChIP_seq/final_data/TEAD4/TEAD4_2c_output_merge.bed"), header=F, stringsAsFactors = F, check.names = F)
TEAD4m$V4=paste0(TEAD4m$V1,"_",TEAD4m$V2,"_",TEAD4m$V3)
write.table(TEAD4m, paste0(primary_folder,"Data_and_code/Fig4.ChIP_seq/final_data/TEAD4/TEAD4_2c_output_merge.bed"), col.names=F, row.names=F, sep="\t", quote=F)

pfm0 = readJASPARMatrix(paste0(primary_folder,"Data_and_code/Fig3.ChromVar/jaspar2024/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.adjust.txt"), matrixClass = "PFM")
bed0 = import(paste0(primary_folder,"Data_and_code/Fig4.ChIP_seq/final_data/TEAD4/TEAD4_2c_output_merge.bed"))
seqlevelsStyle(bed0) = "UCSC"
rse = SummarizedExperiment(rowRanges = bed0)
genome(rse) = "hg38"
motif_ix = matchMotifs(pfm0, rse, genome = BSgenome.Hsapiens.UCSC.hg38, p.cutoff = 1e-04)
gr = rowRanges(rse)
mm_df = cbind(data.frame(ID = mcols(gr)$name), as.matrix(motifMatches(motif_ix)))
mm_df1=mm_df[,c("ID","MA0090.4","MA1121.2","MA0808.1","MA0809.3")]

overlaps1 <- findOverlaps(bed0, griPSC, minoverlap = 1)
mm_df1$overlap_ATACiPSC <- ifelse(seq_along(bed0) %in% queryHits(overlaps1), "Yes", "No")
overlaps2 <- findOverlaps(bed0, grNSC, minoverlap = 1)
mm_df1$overlap_ATACNSC <- ifelse(seq_along(bed0) %in% queryHits(overlaps2), "Yes", "No")
overlaps3 <- findOverlaps(bed0, grNRN, minoverlap = 1)
mm_df1$overlap_ATACNeuron <- ifelse(seq_along(bed0) %in% queryHits(overlaps3), "Yes", "No")

ONECUTm=read.delim(paste0(primary_folder,"Data_and_code/Fig4.ChIP_seq/final_data/ONECUT2/ONECUT_2c_output_merge.bed"), header=F, stringsAsFactors = F, check.names = F)
ONECUTm$V4=paste0(ONECUTm$V1,"_",ONECUTm$V2,"_",ONECUTm$V3)
write.table(ONECUTm, paste0(primary_folder,"Data_and_code/Fig4.ChIP_seq/final_data/ONECUT2/ONECUT_2c_output_merge.bed"), col.names=F, row.names=F, sep="\t", quote=F)

bed1 = import(paste0(primary_folder,"Data_and_code/Fig4.ChIP_seq/final_data/ONECUT2/ONECUT_2c_output_merge.bed"))
seqlevelsStyle(bed1) = "UCSC"
rse = SummarizedExperiment(rowRanges = bed1)
genome(rse) = "hg38"
motif_ix = matchMotifs(pfm0, rse, genome = BSgenome.Hsapiens.UCSC.hg38, p.cutoff = 1e-04)
gr = rowRanges(rse)
mm_df = cbind(data.frame(ID = mcols(gr)$name), as.matrix(motifMatches(motif_ix)))
mm_df2=mm_df[,c("ID","MA0679.3","MA0756.3","MA0757.2")]

overlaps1 <- findOverlaps(bed1, griPSC, minoverlap = 1)
mm_df2$overlap_ATACiPSC <- ifelse(seq_along(bed1) %in% queryHits(overlaps1), "Yes", "No")
overlaps2 <- findOverlaps(bed1, grNSC, minoverlap = 1)
mm_df2$overlap_ATACNSC <- ifelse(seq_along(bed1) %in% queryHits(overlaps2), "Yes", "No")
overlaps3 <- findOverlaps(bed1, grNRN, minoverlap = 1)
mm_df2$overlap_ATACNeuron <- ifelse(seq_along(bed1) %in% queryHits(overlaps3), "Yes", "No")

mm_df1$TF="TEAD4"
mm_df2$TF="ONECUT2"

write.table(mm_df1,gzfile(paste0(path_fig4_data,"TEAD_ChIP_peak_merged_withmotif_overlapATAC.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(mm_df2,gzfile(paste0(path_fig4_data,"ONECUT_ChIP_peak_merged_withmotif_overlapATAC.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

mm_df1=read.delim(paste0(path_fig4_data,"TEAD_ChIP_peak_merged_withmotif_overlapATAC.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
mm_df2=read.delim(paste0(path_fig4_data,"ONECUT_ChIP_peak_merged_withmotif_overlapATAC.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
mm_df1a=mm_df1[which(rowSums(mm_df1[,c(2:5)] == "TRUE")>0),-c(2:5)]
mm_df2a=mm_df2[which(rowSums(mm_df2[,c(2:4)] == "TRUE")>0),-c(2:4)]

both=rbind(mm_df1a,mm_df2a)
write.table(both,gzfile(paste0(path_fig4_data,"TEAD_ONECUT_ChIP_peak_merged_withmotif_overlapATAC.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

# -> for fig 4e

#===============================================================================
# integrate with ABC model to define TEAD-enhancer linked genes and ONECUT-enhancer linked genes
# GO by topgo, KEGG and Reactome by clusterProfiler
setwd(paste0(primary_folder,"Data_and_code/Fig4.TEAD_ONECUT"))
data0=read.delim(paste0(ABC_folder,"EnhancerPredictionsFull_threshold0.02_3cell_eRNAalone.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

maCREe=read.delim(paste0(path_fig4_data,"proportion_chromatin_state.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
maCREe$V4=gsub("-","_",maCREe$V4)
maCREe=maCREe[-grep("romoter",maCREe$promoter_type2),]

avg_expr=read.delim(paste0(primary_folder,"Data_and_code/sc5nRNA/output/gene_normalizedExp_perSample.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
avg_expr$geneName=rownames(avg_expr)

library(topGO)
options(connectionObserver = NULL)
library(org.Hs.eg.db)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)

uni=avg_expr$geneName[which(avg_expr$iPS>0.1)]
background_map <- bitr(uni, fromType = "SYMBOL", toType = "ENTREZID",OrgDb = org.Hs.eg.db)
background_entrez <- unique(background_map$ENTREZID)

ABC=c(0.02, 0.04, 0.08)
GO.table=data.frame()
ReK.table=data.frame()
#indirect by ABC
data0_TEAD4=data0[which(data0$peakID %in% maCREe$V4[which(maCREe$group == "TEAD in iPSC" & maCREe$variable == "FootPrint_ChIP")]),]

for (i in 1:length(ABC)){
quest=intersect(uni,unique(data0_TEAD4$TargetGene[which(data0_TEAD4$CellType=="iPSC" & data0_TEAD4$ABC.Score >= ABC[i])]))
geneList <- factor(as.integer(uni %in% quest))
names(geneList) <- uni
GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, nodeSize=5, annot=annFUN.org, mapping="org.Hs.eg.db", ID="symbol")
allGO = usedGO(object = GOdata)
resultTopGO.elim <- runTest(GOdata, algorithm = "elim", statistic = "Fisher")
GO.table1 = data.table(GenTable( GOdata, Fisher.elim = resultTopGO.elim, topNodes = length(allGO)))
GO.table1$ABCcutoff=ABC[i]
GO.table1$n_gene=length(quest)
GO.table=rbind(GO.table,GO.table1)

gene_map <- bitr(quest, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
genes_entrez <- unique(gene_map$ENTREZID)
kk <- enrichKEGG(gene = genes_entrez, organism = "hsa", universe = background_entrez, pvalueCutoff  = 0.5, qvalueCutoff  = 0.5)
kk_df <- as.data.frame(kk)
kk_df=kk_df[,-c(1:2)]
kk_df$set="KEGG"
react <- enrichPathway(gene = genes_entrez, organism = "human", universe = background_entrez, pvalueCutoff  = 0.5, qvalueCutoff  = 0.5, readable = TRUE)
reactome_df <- as.data.frame(react)
reactome_df$set="Reactome"
ReK.table1=rbind(kk_df,reactome_df)
ReK.table1$ABCcutoff=ABC[i]
ReK.table1$n_gene=length(quest)
ReK.table=rbind(ReK.table,ReK.table1)}

GO.table=GO.table%>%group_by(ABCcutoff)%>%mutate(FDR=p.adjust(as.numeric(Fisher.elim), method="BH"))
write.table(GO.table, gzfile("GO/TEAD4_ChIP_indirect.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)
write.table(ReK.table, gzfile("GO/KEGG_REACTOME_TEAD4_ChIP_indirect.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)

#
ABC=c(0.02, 0.04, 0.08)
GO.table=data.frame()
ReK.table=data.frame()
data0_TEAD4=data0[which(data0$peakID %in% maCREe$V4[which(maCREe$group == "TEAD in iPSC" & maCREe$variable == "FootPrint_ChIP_KD")]),]
for (i in 1:length(ABC)){
  quest=intersect(uni,unique(data0_TEAD4$TargetGene[which(data0_TEAD4$CellType=="iPSC" & data0_TEAD4$ABC.Score >= ABC[i])]))
  geneList <- factor(as.integer(uni %in% quest))
  names(geneList) <- uni
  GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, nodeSize=5, annot=annFUN.org, mapping="org.Hs.eg.db", ID="symbol")
  allGO = usedGO(object = GOdata)
  resultTopGO.elim <- runTest(GOdata, algorithm = "elim", statistic = "Fisher")
  GO.table1 = data.table(GenTable( GOdata, Fisher.elim = resultTopGO.elim, topNodes = length(allGO)))
  GO.table1$ABCcutoff=ABC[i]
  GO.table1$n_gene=length(quest)
  GO.table=rbind(GO.table,GO.table1)
  
  gene_map <- bitr(quest, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  genes_entrez <- unique(gene_map$ENTREZID)
  kk <- enrichKEGG(gene = genes_entrez, organism = "hsa", universe = background_entrez, pvalueCutoff  = 0.5, qvalueCutoff  = 0.5)
  kk_df <- as.data.frame(kk)
  kk_df=kk_df[,-c(1:2)]
  kk_df$set="KEGG"
  react <- enrichPathway(gene = genes_entrez, organism = "human", universe = background_entrez, pvalueCutoff  = 0.5, qvalueCutoff  = 0.5, readable = TRUE)
  reactome_df <- as.data.frame(react)
  reactome_df$set="Reactome"
  ReK.table1=rbind(kk_df,reactome_df)
  ReK.table1$ABCcutoff=ABC[i]
  ReK.table1$n_gene=length(quest)
  ReK.table=rbind(ReK.table,ReK.table1)}
GO.table=GO.table%>%group_by(ABCcutoff)%>%mutate(FDR=p.adjust(as.numeric(Fisher.elim), method="BH"))
write.table(GO.table, gzfile("GO/TEAD4_ChIP_KD_indirect.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)
write.table(ReK.table, gzfile("GO/KEGG_REACTOME_TEAD4_ChIP_KD_indirect.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)

#=========
#ONECUT
uni=avg_expr$geneName[which(avg_expr$Neuron>0.1)]
background_map <- bitr(uni, fromType = "SYMBOL", toType = "ENTREZID",OrgDb = org.Hs.eg.db)
background_entrez <- unique(background_map$ENTREZID)

ABC=c(0.02, 0.04, 0.08)
GO.table=data.frame()
ReK.table=data.frame()
data0_ONECUT=data0[which(data0$peakID %in% maCREe$V4[which(maCREe$group == "ONECUT in Neuron" & maCREe$variable == "FootPrint_ChIP")]),]
for (i in 1:length(ABC)){
  quest=intersect(uni,unique(data0_ONECUT$TargetGene[which(data0_ONECUT$CellType=="Neuron" & data0_ONECUT$ABC.Score >= ABC[i])]))
  geneList <- factor(as.integer(uni %in% quest))
  names(geneList) <- uni
  GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, nodeSize=5, annot=annFUN.org, mapping="org.Hs.eg.db", ID="symbol")
  allGO = usedGO(object = GOdata)
  resultTopGO.elim <- runTest(GOdata, algorithm = "elim", statistic = "Fisher")
  GO.table1 = data.table(GenTable( GOdata, Fisher.elim = resultTopGO.elim, topNodes = length(allGO)))
  GO.table1$ABCcutoff=ABC[i]
  GO.table1$n_gene=length(quest)
  GO.table=rbind(GO.table,GO.table1)
  
  gene_map <- bitr(quest, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  genes_entrez <- unique(gene_map$ENTREZID)
  kk <- enrichKEGG(gene = genes_entrez, organism = "hsa", universe = background_entrez, pvalueCutoff  = 1, qvalueCutoff  = 1)
  kk_df <- as.data.frame(kk)
  kk_df=kk_df[,-c(1:2)]
  kk_df$set="KEGG"
  react <- enrichPathway(gene = genes_entrez, organism = "human", universe = background_entrez, pvalueCutoff  = 1, qvalueCutoff  = 1, readable = TRUE)
  reactome_df <- as.data.frame(react)
  reactome_df$set="Reactome"
  ReK.table1=rbind(kk_df,reactome_df)
  ReK.table1$ABCcutoff=ABC[i]
  ReK.table1$n_gene=length(quest)
  ReK.table=rbind(ReK.table,ReK.table1)}
GO.table=GO.table%>%group_by(ABCcutoff)%>%mutate(FDR=p.adjust(as.numeric(Fisher.elim), method="BH"))
write.table(GO.table, gzfile("GO/ONECUT2_ChIP_indirect.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)
write.table(ReK.table, gzfile("GO/KEGG_REACTOME_ONECUT_ChIP_indirect.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)

#
ABC=c(0.02, 0.04, 0.08)
GO.table=data.frame()
ReK.table=data.frame()
data0_ONECUT=data0[which(data0$peakID %in% maCREe$V4[which(maCREe$group == "ONECUT in Neuron" & maCREe$variable == "FootPrint_ChIP_KD")]),]
for (i in 1:length(ABC)){
  quest=intersect(uni,unique(data0_ONECUT$TargetGene[which(data0_ONECUT$CellType=="Neuron" & data0_ONECUT$ABC.Score >= ABC[i])]))
  geneList <- factor(as.integer(uni %in% quest))
  names(geneList) <- uni
  GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, nodeSize=5, annot=annFUN.org, mapping="org.Hs.eg.db", ID="symbol")
  allGO = usedGO(object = GOdata)
  resultTopGO.elim <- runTest(GOdata, algorithm = "elim", statistic = "Fisher")
  GO.table1 = data.table(GenTable( GOdata, Fisher.elim = resultTopGO.elim, topNodes = length(allGO)))
  GO.table1$ABCcutoff=ABC[i]
  GO.table1$n_gene=length(quest)
  GO.table=rbind(GO.table,GO.table1)
  
  gene_map <- bitr(quest, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  genes_entrez <- unique(gene_map$ENTREZID)
  kk <- enrichKEGG(gene = genes_entrez, organism = "hsa", universe = background_entrez, pvalueCutoff  = 1, qvalueCutoff  = 1)
  kk_df <- as.data.frame(kk)
  kk_df=kk_df[,-c(1:2)]
  kk_df$set="KEGG"
  react <- enrichPathway(gene = genes_entrez, organism = "human", universe = background_entrez, pvalueCutoff  = 1, qvalueCutoff  = 1, readable = TRUE)
  reactome_df <- as.data.frame(react)
  reactome_df$set="Reactome"
  ReK.table1=rbind(kk_df,reactome_df)
  ReK.table1$ABCcutoff=ABC[i]
  ReK.table1$n_gene=length(quest)
  ReK.table=rbind(ReK.table,ReK.table1)}

GO.table=GO.table%>%group_by(ABCcutoff)%>%mutate(FDR=p.adjust(as.numeric(Fisher.elim), method="BH"))
write.table(GO.table, gzfile("GO/ONECUT2_ChIP_KD_indirect.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)
write.table(ReK.table, gzfile("GO/KEGG_REACTOME_ONECUT_ChIP_KD_indirect.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)

#========================
setwd(paste0(primary_folder,"Data_and_code/Fig4.TEAD_ONECUT/GO"))
files=list.files(pattern="_indirect.tsv.gz")
files.names=gsub("_indirect.tsv.gz","",files)
files.names=gsub("KEGG_REACTOME_","",files.names)
ReaKEG=data.frame()
for (i in 1:4){
data=read.delim(files[i], header=T, stringsAsFactors = F)
data=data[which(data$p.adjust<0.05 & data$set != "GOBP" & data$ABCcutoff==0.02),]
if (nrow(data)>0){
data$test=files.names[i]
ReaKEG=rbind(ReaKEG, data)}}

GO=data.frame()
for (i in 5:8){
data1=read.delim(files[i], header=T, stringsAsFactors = F)
data1=data1[which(data1$FDR <0.05 & data1$ABCcutoff==0.02),]
if (nrow(data1)>0){
  data1$test=files.names[i]
  GO=rbind(GO, data1)}}
colnames(GO)[c(1,2,9)]=c("ID","Description","qvalue")
GO$set="GO_BP"
final=bind_rows(ReaKEG,GO)
write.table(final, gzfile("indirect_all_summary.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)
write.table(final, gzfile(paste0(path_fig4_data,"indirect_all_summary.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
# -> for fig ex4f

#===============================================================================
# direct promoter
maCREe=read.delim(paste0(path_fig4_data,"proportion_chromatin_state.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
maCREe$V4=gsub("-","_",maCREe$V4)

aCREanno2=read.delim(paste0(primary_folder,"Data_and_code/Fig2.CRE_feature/CRE_to_gene/filtered_peaks.merged.all.table5gene.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
aCREanno2$V4=gsub("-","_",aCREanno2$V4) #-> link gene to promoter directly
aCREanno2=aCREanno2[which(!is.na(aCREanno2$HGNC)),]
aCREanno2=left_join(aCREanno2, unique(maCREe[,c(1,3)]), by="V4", copy=F)

#direct
GO.table=data.frame()
ReK.table=data.frame()
cell=c("iPS","Neuron")
label=c("TEAD in iPSC", "ONECUT in Neuron")

for (x in 1:2){
  for (z in c("FootPrint_ChIP","FootPrint_ChIP_KD")){
    uni=avg_expr$geneName[which(avg_expr[[cell[x]]]>0.1)]
    background_map <- bitr(uni, fromType = "SYMBOL", toType = "ENTREZID",OrgDb = org.Hs.eg.db)
    background_entrez <- unique(background_map$ENTREZID)
    
    need=maCREe$V4[which(maCREe$group == label[x] & maCREe$variable == z)]
    aCREanno3=aCREanno2[which(aCREanno2$V4 %in% need),]
    direct=intersect(uni,unique(aCREanno3$geneName))
    
    geneList <- factor(as.integer(uni %in% direct))
    names(geneList) <- uni
    GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, nodeSize=5, annot=annFUN.org, mapping="org.Hs.eg.db", ID="symbol")
    allGO = usedGO(object = GOdata)
    resultTopGO.elim <- runTest(GOdata, algorithm = "elim", statistic = "Fisher")
    GO.table1 = data.table(GenTable( GOdata, Fisher.elim = resultTopGO.elim, topNodes = length(allGO)))
    GO.table1$label=label[x]
    GO.table1$variable=z
    GO.table1$n_gene=length(direct)
    GO.table=rbind(GO.table,GO.table1)
    
    gene_map <- bitr(direct, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    genes_entrez <- unique(gene_map$ENTREZID)
    kk <- enrichKEGG(gene = genes_entrez, organism = "hsa", universe = background_entrez, pvalueCutoff  = 1, qvalueCutoff  = 1)
    kk_df <- as.data.frame(kk)
    kk_df=kk_df[,-c(1:2)]
    kk_df$set="KEGG"
    react <- enrichPathway(gene = genes_entrez, organism = "human", universe = background_entrez, pvalueCutoff  = 1, qvalueCutoff  = 1, readable = TRUE)
    reactome_df <- as.data.frame(react)
    reactome_df$set="Reactome"
    ReK.table1=rbind(kk_df,reactome_df)
    ReK.table1$label=label[x]
    ReK.table1$variable=z
    ReK.table1$n_gene=length(direct)
    ReK.table=rbind(ReK.table,ReK.table1)}}
  
GO.table=GO.table%>%group_by(label,variable)%>%mutate(FDR=p.adjust(as.numeric(Fisher.elim), method="BH"))
GO.table$label=gsub("TEAD","TEAD4",GO.table$label)
GO.table$label=gsub("ONECUT","ONECUT2",GO.table$label)
ReK.table$label=gsub("TEAD","TEAD4",ReK.table$label)
ReK.table$label=gsub("ONECUT","ONECUT2",ReK.table$label)
write.table(GO.table, gzfile("GO/Both_ChIP_direct.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)
write.table(ReK.table, gzfile("GO/KEGG_REACTOME_Both_ChIP_direct.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)
  
#===============================================================================
setwd(paste0(primary_folder,"Data_and_code/Fig4.TEAD_ONECUT/GO"))

ReaKEG=read.delim("KEGG_REACTOME_Both_ChIP_direct.tsv.gz", header=T, stringsAsFactors = F)
#ReaKEG=ReaKEG[which(ReaKEG$p.adjust<0.05),]

GO=read.delim("Both_ChIP_direct.tsv.gz", header=T, stringsAsFactors = F)
#GO=GO[which(GO$FDR <0.05),] #-> nothing

colnames(GO)[c(1,2,10)]=c("ID","Description","qvalue")
GO$set="GO_BP"
final=bind_rows(ReaKEG,GO)
final$variable=gsub("FootPrint_","",final$variable)
final$label=sapply(strsplit(final$label, " in"),"[",1)
final$test=paste0(final$label,"_",final$variable)

indirect=read.delim("indirect_all_summary.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
final=final[which(final$ID %in% indirect$ID),]
final=left_join(final, indirect[,c(1,15,16)],by=c("ID","test"),copy=F, suffix=c("","_2"))
final=final[which(!is.na(final$n_gene_2)),]
write.table(final[,c(1:13,16,21,17:20)], gzfile("direct_all_summary.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)
write.table(final[,c(1:13,16,21,17:20)], gzfile(paste0(path_fig4_data,"direct_all_summary.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
# -> for fig ex4h

#===============================================================================
# summary on ABC results
data0=read.delim(paste0(ABC_folder,"EnhancerPredictionsFull_threshold0.02_3cell_eRNAalone.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data0$peakID=gsub("_","-",data0$peakID)
#limit to enhancer alone
aCRE_cell=read.delim(paste0(peak_folder,"all_aCRE.TSS.access.markalone_cuttag.tsv.gz"), header=T, stringsAsFactors = F)

data0a=data0[which(data0$peakID %in% aCRE_cell$peakID[which(aCRE_cell$promoter_type_CT=="enhancer-like")]),]
data0b=data0a%>%group_by(CellType,peakID)%>%dplyr::summarise(ABCscore=sum(ABC.Score))
data0c=spread(data0b, key=1, value=3)
aCRE_cell=left_join(aCRE_cell,data0c, by="peakID",copy=F)

aCRE_celle=aCRE_cell[which(aCRE_cell$promoter_type_CT == "enhancer-like"),]
aCRE_celle.ips=aCRE_celle[,c(1,grep("iPS",colnames(aCRE_celle)))]
aCRE_celle.nsc=aCRE_celle[,c(1,grep("NSC",colnames(aCRE_celle)))]
aCRE_celle.nrn=aCRE_celle[,c(1,grep("Neuron",colnames(aCRE_celle)))]
colnames(aCRE_celle.ips)[c(3,4,11)]=c("TSS","state","ABC_score")
colnames(aCRE_celle.nsc)=colnames(aCRE_celle.ips)
colnames(aCRE_celle.nrn)=colnames(aCRE_celle.ips)
aCRE_celle.ips$sample="iPSC"
aCRE_celle.nsc$sample="NSC"
aCRE_celle.nrn$sample="Neuron"
final=rbind(aCRE_celle.ips[,c(1,3,4,11,12)],aCRE_celle.nsc[,c(1,3,4,11,12)],aCRE_celle.nrn[,c(1,3,4,11,12)])
final$state[which(is.na(final$state))]="Others"
final$state[which(final$state=="CTCF")]="Others"
final$TSS[which(final$TSS==1)]="Transcribed"
final$TSS[which(final$TSS==0)]="Non-transcribed"
final1=final%>%group_by(sample, state)%>%dplyr::summarise(count=n())%>%mutate(percent=count/sum(count), group="All_enhancer")
final2=final[which(!is.na(final$ABC_score)),]%>%group_by(sample, state)%>%dplyr::summarise(count=n())%>%mutate(percent=count/sum(count),group="Predicted_enhancer")
final3=rbind(final1, final2)
final3$state=gsub("_enhancer","",final3$state)

final1b=final%>%group_by(sample, TSS)%>%dplyr::summarise(count=n())%>%mutate(percent=count/sum(count), group="All_enhancer")
final2b=final[which(!is.na(final$ABC_score)),]%>%group_by(sample, TSS)%>%dplyr::summarise(count=n())%>%mutate(percent=count/sum(count),group="Predicted_enhancer")
final3b=rbind(final1b, final2b)

write.table(final3, gzfile(paste0(path_fig4_data,"ABC_summary_state.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(final3b, gzfile(paste0(path_fig4_data,"ABC_summary_transcription.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
# -> for fig ex4f&g











