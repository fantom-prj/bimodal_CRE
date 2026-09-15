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
peak_folder=paste0(primary_folder,"Data_and_code/scATAC/merged_peak/")
motif_folder=paste0(primary_folder,"Data_and_code/Fig2.motif_enrichment/")
path_fig2_data=paste0(primary_folder,"Fig2/data/")

#===============================================================================
setwd(motif_folder)

# TF motif enrichment
# prepare bed files
aCRE_cell=read.delim(paste0(peak_folder,"all_aCRE.TSS.access.markalone_cuttag.tsv.gz"), header=T, stringsAsFactors = F)
aCRE_cell$chr=sapply(strsplit(aCRE_cell$peakID,"-"),"[",1)
aCRE_cell$start=as.numeric(sapply(strsplit(aCRE_cell$peakID,"-"),"[",2))
aCRE_cell$end=as.numeric(sapply(strsplit(aCRE_cell$peakID,"-"),"[",3))
aCRE_cell$V6="."
aCRE_cell1=aCRE_cell[which(aCRE_cell$promoter_type_CT == "enhancer-like"),]

path14=paste0(motif_folder,"F6_meme/input_bed/")

need=c("Primed_enhancer","Repressed_enhancer","Active_enhancer","Promoter","Bivalent_promoter")
need1=c("en_primed","en_repressed","en_active","p_active","p_bivalent")
for (j in 1:length(need)){
  aCRE_celli=aCRE_cell[which(aCRE_cell$iPSC_CT == need[j]),]
  aCRE_cells=aCRE_cell[which(aCRE_cell$NSC_CT == need[j]),]
  aCRE_celln=aCRE_cell[which(aCRE_cell$Neuron_CT == need[j]),]
  
  aCRE_cellia=aCRE_celli[which(aCRE_celli$iPSC_t ==0),]
  aCRE_cellit=aCRE_celli[which(aCRE_celli$iPSC_t ==1),]
  aCRE_cellsa=aCRE_cells[which(aCRE_cells$NSC_t ==0),]
  aCRE_cellst=aCRE_cells[which(aCRE_cells$NSC_t ==1),]
  aCRE_cellna=aCRE_celln[which(aCRE_celln$Neuron_t ==0),]
  aCRE_cellnt=aCRE_celln[which(aCRE_celln$Neuron_t ==1),]
  
  write.table(aCRE_cellia[order(aCRE_cellia$chr,aCRE_cellia$start),c(12,14,15,1,5,16)],gzfile(paste0(path14,"iPSC/iPSC_",need1[j],"xT.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
  write.table(aCRE_cellsa[order(aCRE_cellsa$chr,aCRE_cellsa$start),c(12,14,15,1,6,16)],gzfile(paste0(path14,"NSC/NSC_",need1[j],"xT.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
  write.table(aCRE_cellna[order(aCRE_cellna$chr,aCRE_cellna$start),c(12,14,15,1,7,16)],gzfile(paste0(path14,"NRN/NRN_",need1[j],"xT.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
  write.table(aCRE_cellit[order(aCRE_cellit$chr,aCRE_cellit$start),c(12,14,15,1,5,16)],gzfile(paste0(path14,"iPSC/iPSC_",need1[j],"T.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
  write.table(aCRE_cellst[order(aCRE_cellst$chr,aCRE_cellst$start),c(12,14,15,1,6,16)],gzfile(paste0(path14,"NSC/NSC_",need1[j],"T.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
  write.table(aCRE_cellnt[order(aCRE_cellnt$chr,aCRE_cellnt$start),c(12,14,15,1,7,16)],gzfile(paste0(path14,"NRN/NRN_",need1[j],"T.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)}

# run motif_folder/F6_meme/1.bed2fasta.sh & motif_folder/F6_meme/2.run_sea.sh
# results stored in motif_folder/F6_meme/out

#===============================================================================
#parse MEME results
path15=paste0(motif_folder,"F6_meme/out/")
files = list.files(path=path15, pattern="sea.tsv",recursive = T)
files.names=gsub("SEA_","",files)
files.names=gsub("/sea.tsv","",files.names)
df = data.frame()
for (i in 1:length(files)){
  sea_set = read.table(paste0(path15, files[i]), header=T) %>% select(ID, ALT_ID, ENR_RATIO, LOG_QVALUE) %>% arrange(ID) %>% mutate(Set = files.names[i])
  df = rbind(df, sea_set)}

df = df %>% mutate(Set = factor(Set)) %>% mutate(logQval = -LOG_QVALUE)
write.table(df, paste0(motif_folder, "F6_meme/F6_meme.SEAall.tsv.gz"), sep = "\t", row.names = F, quote = F)

#===================
# extract non-transcribed promoter and enhancer
expressed=read.delim(paste0(primary_folder,"Fig3/data/expressed_400motif.tsv.gz"))

enrich=read.delim(paste0(motif_folder,"F6_meme/F6_meme.SEAall.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
unique(enrich$Set)
enrich=enrich[grep("xT",enrich$Set),]
enrich$Set=gsub("xT","",enrich$Set)
enrich$Set1=sapply(strsplit(enrich$Set,"_"),"[",1)
enrich$Set2=sapply(strsplit(enrich$Set,"_"),"[",2)
enrich$Set3=sapply(strsplit(enrich$Set,"_"),"[",3)
enrich$Set4=paste0(enrich$Set3,"_",enrich$Set2,"_",enrich$Set1)

enrich1=enrich[-grep("active",enrich$Set),] #select from non-active CREs
#enrich=enrich[c(grep("repressed",enrich$Set),grep("primed",enrich$Set)),]
need=unique(enrich1$ID[which(enrich1$ENR_RATIO >3 & enrich1$logQval >2)])
enricha=enrich[which(enrich$ID %in% need),]
enricha$ID2=paste0(enricha$ALT_ID,"_",enricha$ID)
enricha=enricha%>%group_by(ID2)%>%dplyr::mutate(scaled_value=ENR_RATIO/max(ENR_RATIO))
enricha=enricha[which(enricha$ID %in% unique(expressed$motifID)),]

#select generally stronger ones if same motif names
enrichc=enricha%>%group_by(ALT_ID,ID)%>%dplyr::summarise(ENR_RATIO=max(ENR_RATIO))%>%slice_max(ENR_RATIO)
enricha=enricha[which(enricha$ID %in% enrichc$ID),]

write.table(enricha,gzfile(paste0(path_fig2_data,"F_enrichemnt_for_xT.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)
#-> for fig 2h

#===============================================================================
# identify TF motif from the list of repressed enhancer
library(motifmatchr)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SummarizedExperiment)
set.seed(2017)
library(rtracklayer)
library(GenomicRanges)

both_table=read.delim(paste0(path_fig2_data,"CRE_463866_sample_expression_ATAC_TSS.plot.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
mboth_table=reshape2::melt(both_table, id=c(1:3,6:13))
mboth_table$sample=factor(mboth_table$sample, levels=c("iPSC","NSC","Neuron"))
mboth_table5=spread(unique(mboth_table[which(mboth_table$promotertypeCT=="enhancer-like"),c(1,3,9,10)]),key=3, value=4)

mboth_table5[mboth_table5=="Genomic"]=0
mboth_table5[mboth_table5=="CTCF"]=0
mboth_table5[mboth_table5=="Primed"]=2 #blue
mboth_table5[mboth_table5=="Repressed"]=3 #green
mboth_table5[mboth_table5=="Active"]=4 #red
mboth_table5[,c(3:5)] <- mboth_table5[,c(3:5)] %>% mutate_if(is.character, as.numeric)
mboth_table5$active_count=rowSums(mboth_table5[,c(3:5)]==4)

need=unique(both_table$peakID[which(both_table$state == "Repressed" & mboth_table$promotertypeCT=="enhancer-like")]) # repressed enhancer state in any cell type
mboth_table2=mboth_table5[which(mboth_table5$peakID %in% need),]
write.table(mboth_table2,gzfile(paste0(motif_folder,"repressedEnhancer_motif/all_repressed_enhancer_state.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

need2=unique(both_table$peakID[which(both_table$state == "Primed" & mboth_table$promotertypeCT=="enhancer-like")]) # primed enhancer state in any cell type
mboth_table3=mboth_table5[which(mboth_table5$peakID %in% need2),]

pfm0 <- readJASPARMatrix(paste0(primary_folder,"Data_and_code/Fig3.ChromVar/jaspar2024/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.adjust.txt"), matrixClass = "PFM")
bed0 <- import(paste0(peak_folder,"peaks.merged.all.bed.gz"))
bed1 <- bed0[bed0@elementMetadata@listData[["name"]] %in% need, ]
seqlevelsStyle(bed1) <- "UCSC"
rse <- SummarizedExperiment(rowRanges = bed1)
genome(rse) <- "hg38"
motif_ix <- matchMotifs(pfm0, rse, genome = BSgenome.Hsapiens.UCSC.hg38, p.cutoff = 1e-04)
gr <- rowRanges(rse)
mm_df <- cbind(
  data.frame(ID = mcols(gr)$name),
  as.matrix(motifMatches(motif_ix)))
motif_id=read.delim(paste0(primary_folder,"Data_and_code/Fig3.ChromVar/jaspar2024/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.ID.tsv"), header=F, stringsAsFactors = F)
mmm_df=reshape2::melt(mm_df,id=1)
mmm_df=left_join(mmm_df,motif_id,by=c("variable"="V1"),copy=F)
mmm_df=mmm_df[which(mmm_df$value == "TRUE"),]
write.table(mmm_df,gzfile(paste0(motif_folder,"repressedEnhancer_motif/all_repressed_enhancer_motifmatch.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

mm_df1=mm_df[,c("ID","MA1581.2","MA0657.2")]
colnames(mm_df1)[c(2:3)]=c("ZBTB6","KLF13")
mboth_table2=left_join(mboth_table2, mm_df1, by=c("peakID"="ID"),copy=F)


mboth_table_enhancer=right_join(unique(mboth_table[,c(1:11)]),mboth_table5[,c(1,6)],by="peakID",copy=F)
mboth_table_enhancer=left_join(mboth_table_enhancer, mm_df1, by=c("peakID"="ID"),copy=F)

#=
kk0=mboth_table_enhancer%>%group_by(sample,group2)%>%dplyr::summarise(state="All", count=n(),active_count = length(which(active_count > 0)))
kk0$percent=kk0$active_count/kk0$count*100
kk1=mboth_table_enhancer%>%group_by(sample,group2,state)%>%dplyr::summarise(count=n(),active_count = length(which(active_count > 0)))
kk1$percent=kk1$active_count/kk1$count*100
#need3=unique(both_table$peakID[which(both_table$state == "Repressed" & mboth_table$promotertypeCT=="enhancer-like" & mboth_table$sample=="Neuron")]) # repressed enhancer state in neuron
kk2=mboth_table_enhancer%>%group_by(sample,group2,state,ZBTB6)%>%dplyr::summarise(count=n(),active_count = length(which(active_count > 0)))
kk2$percent=kk2$active_count/kk2$count*100
#need4=unique(both_table$peakID[which(both_table$state == "Repressed" & mboth_table$promotertypeCT=="enhancer-like" & mboth_table$sample=="iPSC")]) # repressed enhancer state in iPSC
kk3=mboth_table_enhancer%>%group_by(sample,group2,state,KLF13)%>%dplyr::summarise(count=n(),active_count = length(which(active_count > 0)))
kk3$percent=kk3$active_count/kk3$count*100

kk0$value="TRUE"
kk1$value="TRUE"
kk2=kk2[which(kk2$sample == "Neuron" & kk2$state == "Repressed"),c(1:3,5:7,4)]
kk3=kk3[which(kk3$sample == "iPSC" & kk3$state == "Repressed"),c(1:3,5:7,4)]
colnames(kk2)[7]="value"
colnames(kk3)[7]="value"
kk2$state="Repressed_ZBTB6"
kk3$state="Repressed_KLF13"
kk2$value=as.character(kk2$value)
kk3$value=as.character(kk3$value)
kk=rbind(kk0,kk1,kk2,kk3)
kk$group2[which(kk$group2 == "With transcription")]="Transcribed enhancer"
kk$group2[which(kk$group2 == "Without transcription")]="Non-transcribed enhancer"

kk$other_count=kk$count-kk$active_count

write.table(kk,gzfile(paste0(path_fig2_data,"heatmap.markalone_ZBTB6_KLF13.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)
#-> for fig 2i

#========




















# code not used in the study
#===============================================================================
#===============================================================================
# extract transcribed promoter and enhancer
expressed=read.delim(paste0(primary_folder,"Fig3/data/expressed_400motif.tsv.gz"))

enrich=read.delim(paste0(motif_folder,"F6_meme/F6_meme.SEAall.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
unique(enrich$Set)
enrich=enrich[-grep("xT",enrich$Set),]
enrich$Set=gsub("T","",enrich$Set)
enrich$Set1=sapply(strsplit(enrich$Set,"_"),"[",1)
enrich$Set2=sapply(strsplit(enrich$Set,"_"),"[",2)
enrich$Set3=sapply(strsplit(enrich$Set,"_"),"[",3)
enrich$Set4=paste0(enrich$Set3,"_",enrich$Set2,"_",enrich$Set1)

enrich1=enrich[-grep("active",enrich$Set),] #select from non-active CREs
#enrich=enrich[c(grep("repressed",enrich$Set),grep("primed",enrich$Set)),]
need=unique(enrich1$ID[which(enrich1$ENR_RATIO >3 & enrich1$logQval >2 & enrich1$Set3=="repressed")])
enricha=enrich[which(enrich$ID %in% need),]
enricha$ID2=paste0(enricha$ALT_ID,"_",enricha$ID)
enricha=enricha%>%group_by(ID2)%>%dplyr::mutate(scaled_value=ENR_RATIO/max(ENR_RATIO))
enricha=enricha[which(enricha$ID %in% unique(expressed$motifID)),]

#select generally stronger ones if same motif names
enrichc=enricha%>%group_by(ALT_ID,ID)%>%dplyr::summarise(ENR_RATIO=max(ENR_RATIO))%>%slice_max(ENR_RATIO)
enricha=enricha[which(enricha$ID %in% enrichc$ID),]

#
enrichbscale=data.frame(spread(enricha[,c(2,10,12)], key=2, value=3))
rownames(enrichbscale)=enrichbscale$ALT_ID
enrichbscale=enrichbscale[,-1]
row_order <- rownames(enrichbscale)[hclust(dist(enrichbscale))$order]

enricha$ALT_ID=factor(enricha$ALT_ID,levels=row_order)
enricha$Set5=paste0(enricha$Set3,"_",enricha$Set2)
enricha$Set5=gsub("_en","\nenhancer",enricha$Set5)
enricha$Set5=gsub("_p","\npromoter",enricha$Set5)
enricha$Set5=factor(enricha$Set5, levels=c("active\npromoter","bivalent\npromoter","active\nenhancer","repressed\nenhancer","primed\nenhancer"))
enricha$Set1=factor(enricha$Set1, levels=c("iPSC","NSC","NRN"))
ggplot() + 
  facet_wrap(~Set5, ncol=5)+
  labs(x=NULL, y =NULL, title="Motif enrichment of non-transcribed CREs", fill="Scaled\nenrichment")+
  scale_y_discrete(position = "right", limits=rev)+
  geom_tile(data=enricha, mapping=aes(x=Set1, fill=scaled_value, y=ALT_ID), size=0.6)+
  scale_fill_gradient2(low="white",high="firebrick3")+
  theme1+theme(panel.grid.major = element_blank(),
               axis.ticks = element_blank(), axis.line = element_blank(),
               axis.text.x = element_text(color ="black",angle=90,hjust=1, vjust=0.5),
               axis.text.y = element_text(color ="black", size=5))



#===============================
# extract all CRE without separating Tx / T

enrich=read.delim("/osc-fs_home/yip/F6_meme_v2/CRE/sea.all.tsv", header=T, stringsAsFactors = F, check.names = F)
unique(enrich$Cell)
enrich$Set1=sapply(strsplit(enrich$Cell,"_"),"[",1)
enrich$Set2=sapply(strsplit(enrich$Cell,"_"),"[",2)
enrich$Set3=sapply(strsplit(enrich$Cell,"_"),"[",3)
enrich$Set4=paste0(enrich$Set3,"_",enrich$Set2,"_",enrich$Set1)

enrich1=enrich[-grep("active",enrich$Cell),] #select from non-active CREs
#enrich=enrich[c(grep("repressed",enrich$Set),grep("primed",enrich$Set)),]
need=unique(enrich1$ID[which(enrich1$ENR_RATIO >3 & enrich1$logQval >2 & enrich1$Set3=="repressed")])
enricha=enrich[which(enrich$ID %in% need),]
enricha$ID2=paste0(enricha$ALT_ID,"_",enricha$ID)
enricha=enricha%>%group_by(ID2)%>%dplyr::mutate(scaled_value=ENR_RATIO/max(ENR_RATIO))
enricha=enricha[which(enricha$ID %in% unique(expressed$motifID)),]

#select generally stronger ones if same motif names
enrichc=enricha%>%group_by(ALT_ID,ID)%>%dplyr::summarise(ENR_RATIO=max(ENR_RATIO))%>%slice_max(ENR_RATIO)
enricha=enricha[which(enricha$ID %in% enrichc$ID),]

#
enrichbscale=data.frame(spread(enricha[,c(2,10,12)], key=2, value=3))
rownames(enrichbscale)=enrichbscale$ALT_ID
enrichbscale=enrichbscale[,-1]
row_order <- rownames(enrichbscale)[hclust(dist(enrichbscale))$order]

enricha$ALT_ID=factor(enricha$ALT_ID,levels=row_order)
enricha$Set5=paste0(enricha$Set3,"_",enricha$Set2)
enricha$Set5=gsub("_en","\nenhancer",enricha$Set5)
enricha$Set5=gsub("_p","\npromoter",enricha$Set5)
enricha$Set5=factor(enricha$Set5, levels=c("active\npromoter","bivalent\npromoter","active\nenhancer","repressed\nenhancer","primed\nenhancer"))
enricha$Set1=factor(enricha$Set1, levels=c("iPSC","NSC","NRN"))
ggplot() + 
  facet_wrap(~Set5, ncol=5)+
  labs(x=NULL, y =NULL, title="Motif enrichment of non-transcribed CREs", fill="Scaled\nenrichment")+
  scale_y_discrete(position = "right", limits=rev)+
  geom_tile(data=enricha, mapping=aes(x=Set1, fill=scaled_value, y=ALT_ID), size=0.6)+
  scale_fill_gradient2(low="white",high="firebrick3")+
  theme1+theme(panel.grid.major = element_blank(),
               axis.ticks = element_blank(), axis.line = element_blank(),
               axis.text.x = element_text(color ="black",angle=90,hjust=1, vjust=0.5),
               axis.text.y = element_text(color ="black", size=5))


#tobias score in inactive enhancer?
tead_FP1=read.delim("tobias_chip_tead4_cluster1.tsv.gz", header=T, stringsAsFactors = F)
aCRE_cell=read.delim("/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/JC_merged/all_aCRE.TSS.access.markalone_cuttag.tsv", header=T, stringsAsFactors = F, check.names = F)
tead_FP1=left_join(tead_FP1, aCRE_cell[,c(1,8)], by=c("V10"="peakID"),copy=F)
tead_FP1%>%group_by(iPSC_CT,ChIP)%>%dplyr::summarise(median=median(V13))
write.table(tead_FP1,gzfile("tobias_chip_tead4_cluster1.tsv.gz"), row.names=F, col.names=T, sep="\t", quote=F)

#===============================================================================
#tobias for ZBTB6 -> go to fig2???
path1="/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/subset_BAM/BINDetect_compare_output/ZBTB6_MA1581.2/beds/"
ZBTB6_FPb1=read.delim(paste0(path1,"ZBTB6_MA1581.2_out_neuron_immature_CB_footprints_bound.bed"), header=F, stringsAsFactors = F)
ZBTB6_FPub1=read.delim(paste0(path1,"ZBTB6_MA1581.2_out_neuron_immature_CB_footprints_unbound.bed"), header=F, stringsAsFactors = F)
ZBTB6_FPb2=read.delim(paste0(path1,"ZBTB6_MA1581.2_out_neuron_mature_CB_footprints_bound.bed"), header=F, stringsAsFactors = F)
ZBTB6_FPb3=read.delim(paste0(path1,"ZBTB6_MA1581.2_out_neuron_mature_2_CB_footprints_bound.bed"), header=F, stringsAsFactors = F)
ZBTB6_FPb4=read.delim(paste0(path1,"ZBTB6_MA1581.2_out_neuron_progenitor_n_schwann_CB_footprints_bound.bed"), header=F, stringsAsFactors = F)

motif_pos=union(unique(ZBTB6_FPb1$V10),unique(ZBTB6_FPub1$V10))
need=union(union(union(unique(ZBTB6_FPb1$V10),unique(ZBTB6_FPb2$V10)),unique(ZBTB6_FPb3$V10)),unique(ZBTB6_FPb4$V10))
need1=union(union(union(unique(ZBTB6_FPb1$V10[which(ZBTB6_FPb1$V13>0.2)]),unique(ZBTB6_FPb2$V10[which(ZBTB6_FPb2$V13>0.2)])),unique(ZBTB6_FPb3$V10[which(ZBTB6_FPb3$V13>0.2)])),unique(ZBTB6_FPb4$V10[which(ZBTB6_FPb4$V13>0.2)]))

aCREe=read.delim("/analysisdata/fantom6/Interactome/KD_ATAC/BAM/Run2/edgeR/tobias_enhancer_peaks.merged.all.markalone_k16.subclass.tsv",header=T, check.names = F, stringsAsFactors = F)
aCREe$ZBTB6_motif="No"
aCREe$ZBTB6_motif[which(aCREe$V4 %in% motif_pos)]="Yes"
aCREe$ZBTB6_FootPrint="No"
aCREe$ZBTB6_FootPrint[which(aCREe$V4 %in% need)]="Yes"
aCREe$ZBTB6_FootPrint_str="No"
aCREe$ZBTB6_FootPrint_str[which(aCREe$V4 %in% need1)]="Yes"

kk=aCREe%>%group_by(k16nrn_all,ZBTB6_motif,ZBTB6_FootPrint)%>%dplyr::summarise(count=n())
kk1=aCREe%>%group_by(k16nrn_all,ZBTB6_motif,ZBTB6_FootPrint_str)%>%dplyr::summarise(count=n())

kka=aCREe%>%group_by(tCRE,k16nrn_all,ZBTB6_motif,ZBTB6_FootPrint)%>%dplyr::summarise(count=n())

#not use
#===============================================================================
#tobias of ZBTB6
cluster_info=read.delim("/analysisdata/fantom6/Interactome/single_cell_wallace/Figure/Fig1/cluster_info.tsv",header=T, stringsAsFactors = F, check.names = F)

path1="/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/subset_BAM/BINDetect_compare_output/ZBTB6_MA1581.2/beds/"
files=list.files(path=path1, pattern=".bed")
files=files[-1]
files.names=gsub("ZBTB6_MA1581.2_out_","",files)
files.names=unique(sapply(strsplit(files.names,"_CB_footprints"),"[",1))
bind0=data.frame()
for (i in 1: length(files.names)){
  files1=files[grep(files.names[i], files)]
  bind1=read.delim(paste0(path1,files1[1]),header=F, stringsAsFactors = F)
  bind2=read.delim(paste0(path1,files1[2]),header=F, stringsAsFactors = F)
  bind3=rbind(bind1[,c(5,10,13)],bind2[,c(5,10,13)])
  bind3=bind3%>%group_by(V10)%>%dplyr::slice_max(V13)
  bind3=bind3%>%group_by(V10)%>%dplyr::slice_sample(n=1)
  bind3$cell=files.names[i]
  bind0=rbind(bind0, bind3)}
write.table(bind0,gzfile(paste0(path1,"all_celltype_score.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
#replace spread with dcast of data.table
bind0=read.delim(paste0(path1,"all_celltype_score.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
bind0=bind0%>%group_by(cell)%>%dplyr::mutate(nor_score=V13/max(V13))
library(data.table)
dt <- as.data.table(bind0)
colnames(dt)[c(2,3,4)] <- c("id","value","key")
bind00 <- dcast(dt, id ~ key, value.var="value", fun.aggregate = function(x) x[1])

rownames(bind00)=bind00$id
data1=bind00[,c(2:11)]
data2=data.frame(data1[which(rowSums(data1>0.5)>0),])
data2=data2[,cluster_info$cluster]
colnames(data2)=cluster_info$label
data2[data2>3]=3

library(pheatmap)
out1=pheatmap(data2, fontsize = 6, fontsize_col = 6, fontsize_row = 4 , 
              color = colorRampPalette(c("white", "firebrick3"))(50), 
              #breaks=breaks,
              main = "ZBTB6 footprinting score" , 
              treeheight_row = 25, 
              angle_col = 90,
              cluster_cols = F,
              show_rownames = F)


#find repress enhancer in neuron
bind0=read.delim(paste0(path1,"all_celltype_score.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
aCRE_cell=read.delim("/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/JC_merged/all_aCRE.TSS.access.markalone_cuttag.tsv", header=T, stringsAsFactors = F, check.names = F)
bind0=bind0[which(bind0$V10 %in% aCRE_cell$peakID[which(aCRE_cell$Neuron_CT == "Repressed_enhancer")]),]
dt <- as.data.table(bind0)
colnames(dt)[c(2,3,4)] <- c("id","value","key")
bind00 <- dcast(dt, id ~ key, value.var="value", fun.aggregate = function(x) x[1])

rownames(bind00)=bind00$id
data1=bind00[,c(2:11)]
data2=data.frame(data1[which(rowSums(data1>0)>0),])
data2=data2[,cluster_info$cluster]
colnames(data2)=cluster_info$label
data2[data2>2]=2
library(pheatmap)
out1=pheatmap(data2, fontsize = 6, fontsize_col = 6, fontsize_row = 4 , 
              color = colorRampPalette(c("white", "firebrick3"))(50), 
              #breaks=breaks,
              main = "ZBTB6 footprinting score" , 
              treeheight_row = 25, 
              angle_col = 90,
              cluster_cols = F,
              show_rownames = F)

#====
path1="/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/subset_BAM/BINDetect_compare_output/TEAD4_MA0809.3/beds/"
files=list.files(path=path1, pattern=".bed")
files=files[-1]
files=files[-grep("hit.bed.gz",files)]
files.names=gsub("TEAD4_MA0809.3_out_","",files)
files.names=unique(sapply(strsplit(files.names,"_CB_footprints"),"[",1))
bind0=data.frame()
for (i in 1: length(files.names)){
  files1=files[grep(files.names[i], files)]
  bind1=read.delim(paste0(path1,files1[1]),header=F, stringsAsFactors = F)
  bind2=read.delim(paste0(path1,files1[2]),header=F, stringsAsFactors = F)
  bind3=rbind(bind1[,c(5,10,13)],bind2[,c(5,10,13)])
  bind3=bind3%>%group_by(V10)%>%dplyr::slice_max(V13)
  bind3=bind3%>%group_by(V10)%>%dplyr::slice_sample(n=1)
  bind3$cell=files.names[i]
  bind0=rbind(bind0, bind3)}
write.table(bind0,gzfile(paste0(path1,"all_celltype_score.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
#replace spread with dcast of data.table
bind0=read.delim(paste0(path1,"all_celltype_score.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
bind0=bind0%>%group_by(cell)%>%dplyr::mutate(nor_score=V13/max(V13))
library(data.table)
dt <- as.data.table(bind0)
colnames(dt)[c(2,3,4)] <- c("id","value","key")
bind00 <- dcast(dt, id ~ key, value.var="value", fun.aggregate = function(x) x[1])

rownames(bind00)=bind00$id
data1=bind00[,c(2:11)]
data2=data.frame(data1[which(rowSums(data1>0.5)>0),])
data2=data2[,cluster_info$cluster]
colnames(data2)=cluster_info$label
data2[data2>3]=3

library(pheatmap)
out1=pheatmap(data2, fontsize = 6, fontsize_col = 6, fontsize_row = 4 , 
              color = colorRampPalette(c("white", "firebrick3"))(50), 
              #breaks=breaks,
              main = "TEAD4 footprinting score" , 
              treeheight_row = 25, 
              angle_col = 90,
              cluster_cols = F,
              show_rownames = F)

#====
path1="/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/subset_BAM/BINDetect_compare_output/ONECUT2_MA0756.3/beds/"
files=list.files(path=path1, pattern=".bed")
files=files[-1]
files.names=gsub("ONECUT2_MA0756.3_out_","",files)
files.names=unique(sapply(strsplit(files.names,"_CB_footprints"),"[",1))
bind0=data.frame()
for (i in 1: length(files.names)){
  files1=files[grep(files.names[i], files)]
  bind1=read.delim(paste0(path1,files1[1]),header=F, stringsAsFactors = F)
  bind2=read.delim(paste0(path1,files1[2]),header=F, stringsAsFactors = F)
  bind3=rbind(bind1[,c(5,10,13)],bind2[,c(5,10,13)])
  bind3=bind3%>%group_by(V10)%>%dplyr::slice_max(V13)
  bind3=bind3%>%group_by(V10)%>%dplyr::slice_sample(n=1)
  bind3$cell=files.names[i]
  bind0=rbind(bind0, bind3)}
write.table(bind0,gzfile(paste0(path1,"all_celltype_score.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
#replace spread with dcast of data.table
bind0=read.delim(paste0(path1,"all_celltype_score.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
bind0=bind0%>%group_by(cell)%>%dplyr::mutate(nor_score=V13/max(V13))
library(data.table)
dt <- as.data.table(bind0)
colnames(dt)[c(2,3,4)] <- c("id","value","key")
bind00 <- dcast(dt, id ~ key, value.var="value", fun.aggregate = function(x) x[1])

rownames(bind00)=bind00$id
data1=bind00[,c(2:11)]
data2=data.frame(data1[which(rowSums(data1>0.5)>0),])
data2=data2[,cluster_info$cluster]
colnames(data2)=cluster_info$label
data2[data2>3]=3

library(pheatmap)
out1=pheatmap(data2, fontsize = 6, fontsize_col = 6, fontsize_row = 4 , 
              color = colorRampPalette(c("white", "firebrick3"))(50), 
              #scale = "column",
              #breaks=breaks,
              main = "ONECUT2 footprinting score" , 
              treeheight_row = 25, 
              angle_col = 90,
              cluster_cols = F,
              show_rownames = F)









