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
CRE_folder=paste0(primary_folder,"Data_and_code/Fig2.CRE_feature/")
path_fig2_data=paste0(primary_folder,"Fig2/data/")


#===============================================================================
#get the bed12 for all possible transcript, including long read data
setwd("/osc-fs_home/yip/CFC_seq_paper_fig_data/code_n_data/SALA/Neuron_THP1_full/transcript/all_gtf_file/")
system("zcat table5pENST.gtf.gz | bedparse gtf2bed | gzip > table5pENST.bed12.bed.gz")

#===============================================================================
#get the 5'end of all possible transcript, including long read data
setwd("/osc-fs_home/yip/CFC_seq_paper_fig_data/code_n_data/SALA/Neuron_THP1_full/transcript/all_gtf_file/")
options(scipen=999)
t5bed=read.delim("table5pENST.bed12.bed.gz", header=F, stringsAsFactors = F)
t5bed$V3[which(t5bed$V6 == "+")]=t5bed$V2[which(t5bed$V6 == "+")]+1
t5bed$V2[which(t5bed$V6 == "-")]=t5bed$V3[which(t5bed$V6 == "-")]-1
write.table(t5bed[order(t5bed$V1,t5bed$V2),c(1:6)], gzfile(paste0(CRE_folder,"CRE_to_gene/CFC_SALA_table5pENST.5n.bed.gz")), col.names=F, row.names = F, sep="\t", quote=F)

#===============================================================================
#bedtools
# include any GENCODE v39 transcripts 5'n and finalized SALA longread transcript 5'n
setwd(peak_folder)
system(paste0("bedtools intersect -wa -wb -a filtered_peaks.merged.all.bed.gz -b ",CRE_folder,"CRE_to_gene/CFC_SALA_table5pENST.5n.bed.gz | gzip > filtered_peaks.merged.all.table5pENST.bed.gz"))

# include +- 500 bp closest if promoter doesnt get a hit -> NOT use
system(paste0("bedtools closest -a filtered_peaks.merged.all.bed -b ",CRE_folder,"CRE_to_gene/CFC_SALA_table5pENST.5n.bed.gz | gzip > filtered_peaks.merged.all.closest_table5pENST.bed.gz"))

# identify tCRE without ATAC peak overlap
SCAFE_path=paste0(primary_folder,"Data_and_code/SCAFE/sc5end/aggregate/run_full/out/annotate/sc5end.iPSC_NSC_Neuron/")
system(paste0("bedtools intersect -wa -wb -a " ,SCAFE_path, "bed/sc5end.iPSC_NSC_Neuron.CRE.coord.bed.gz -b filtered_peaks.merged.all.bed.gz -v | gzip > sc5end.iPSC_NSC_Neuron.CRE_noATAC.bed.gz"))
system("bedtools merge -i sc5end.iPSC_NSC_Neuron.CRE_noATAC.bed.gz -c 4 -o collapse | gzip > sc5end.iPSC_NSC_Neuron.CRE_noATAC_merged.bed.gz")


#===add tCRE
tCRE=read.delim(paste0(SCAFE_path,"log/sc5end.iPSC_NSC_Neuron.CRE.info.tsv.gz"), header=T, stringsAsFactors = F)
tCRE_noATAC=read.delim("sc5end.iPSC_NSC_Neuron.CRE_noATAC_merged.bed.gz", header=F, stringsAsFactors = F)
tCRE_noATAC1=separate_rows(tCRE_noATAC,V4,sep=",")
tCRE1=tCRE[which(tCRE$CREID %in% tCRE_noATAC1$V4),]
tCRE1gene=tCRE1[which(tCRE1$typeStr == "gene_tss" & tCRE1$regionType != "intron"),]


#count ATAC back to tCRE without ATAC peak overlap

peak_a <- getPeaks("sc5end.iPSC_NSC_Neuron.CRE_noATAC_merged.bed.gz", sort_peaks = TRUE)
IDList <- c("iPS","NSC","NRN")
md.List <- list()
for(i in IDList){
  tmp <- read.table(
    file = paste0("/analysisdata/fantom6/Interactome/scATACcellranger_Kouno/",i,"/outs/singlecell.csv"), 
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  md.List[[i]] <- tmp}

## passed_filters => number of non-duplicate, usable read-pairs i.e. "fragments"
for(i in IDList){md.List[[i]] <- md.List[[i]][md.List[[i]]$passed_filters > 500, ]}

## create fragment objects
frags.List <- list()
for(i in IDList){
  tmp <- CreateFragmentObject(
    path = paste0("/analysisdata/fantom6/Interactome/scATACcellranger_Kouno/",i,"/outs/fragments.tsv.gz"),
    cells = rownames(md.List[[i]])
  )
  frags.List[[i]] <- tmp
}

## Count fragments in peaks (takes time...) x core ~30 min for 18 samples, 6 cores ~ 1h
counts.List <- list()
for(i in IDList){
  tmp <-  FeatureMatrix(
    fragments = frags.List[[i]],
    features = peak_a,
    cells = rownames(md.List[[i]])
  )
  counts.List[[i]] <- tmp
}
saveRDS(counts.List, file = "/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/counts.List_tCREnoATAC.3samples.rds")

k100=as.data.frame(counts.List[["iPS"]])
gene_sum <- data.frame(rowSums(k100))


#===============================================================================
# add transcript expression 
t_bambu=read.delim("/osc-fs_home/yip/CFC_seq_paper_fig_data/code_n_data/transcript_model_analyses_Fig3/bambu_long_t5_partialYes.ENST/transcript.count.matrix.tsv.gz",header=T, stringsAsFactors = F, check.names = F)
rownames(t_bambu)=t_bambu$TXNAME
t_bambu=t_bambu[,c(2:7)]
t_bambu <- sweep(t_bambu, 2, colSums(t_bambu), "/") * 1000000
t_bambu$exp=rowMeans(t_bambu)
t_bambu$tID=rownames(t_bambu)
#aCREanno=read.delim("/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/JC_merged/filtered_peaks.merged.all.table5transcript.bed.gz", header=F, stringsAsFactors = F)
#aCREanno=left_join(aCREanno,t_bambu[,c(6,7)], by=c("V10"="tID"))

aCREanno=read.delim(paste0(CRE_folder,"CRE_to_gene/filtered_peaks.merged.all.table5pENST.bed.gz"), header=F, stringsAsFactors = F)
aCREanno=left_join(aCREanno,t_bambu[,c(8,7)], by=c("V10"="tID"))

#connected to geneID
gene_t=read.delim(paste0(CRE_folder,"CRE_to_gene/CFC_SALA_table5pENST.transcript.gene.class.tsv.gz"), header=T, check.names = F, stringsAsFactors = F)
gene_t=gene_t[,c(1:3)]
HGNC=read.delim(paste0(primary_folder,"resources/gencode.v39.metadata.HGNC.gz"), header=F, stringsAsFactors = F)
HGNC=left_join(HGNC, gene_t, by=c("V1"="transcript"),copy=F)
HGNC1=unique(HGNC[,c(2:4)])
HGNC1=HGNC1[!is.na(HGNC1$gene),]
colnames(HGNC1)=c("HGNC","HGNCID","gene")
aCREanno=left_join(aCREanno, gene_t, by=c("V10"="transcript"),copy=F)
aCREanno=left_join(aCREanno, HGNC1, by="gene",copy=F)

aCREanno$exp[which(is.na(aCREanno$exp))]=0
aCREanno1=aCREanno%>%group_by(V4,gene)%>%dplyr::summarise(exp_gene=sum(exp))
aCREanno1=left_join(aCREanno1, unique(aCREanno[,c(12,14:17)]), by="gene",copy=F)
aCREanno1$geneName=aCREanno1$HGNC
aCREanno1$geneName[which(is.na(aCREanno1$geneName))]=aCREanno1$gene[which(is.na(aCREanno1$geneName))]

# add gene expression (use single cell expression level)
exp=read.delim(paste0(primary_folder,"Data_and_code/metacell/RNA_expression_per_cluster.matrix.tsv.gz"))
exp$scGene_exp=rowMeans(exp[,c(1:10)])
aCREanno1$simple_geneID=sapply(strsplit(aCREanno1$gene,"\\."),"[",1)
aCREanno1=left_join(aCREanno1,exp[,c(12,13)],by=c("simple_geneID"="geneID"),copy=F)

aCREanno2=aCREanno1%>%group_by(V4)%>%dplyr::slice_max(exp_gene)
write.table(aCREanno,gzfile(paste0(CRE_folder,"CRE_to_gene/filtered_peaks.merged.all.table5transcript.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(aCREanno2,gzfile(paste0(CRE_folder,"CRE_to_gene/filtered_peaks.merged.all.table5gene.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
#peakID centric collapsed gene version
aCREanno3=aCREanno1%>%group_by(V4)%>%dplyr::summarise(exp_gene=paste(exp_gene, collapse=";"),geneID=paste(gene, collapse=";"),gene_type=paste(gene_type, collapse=";"), HGNC=paste(HGNC, collapse=";"),HGNCID=paste(HGNCID, collapse=";"),geneName=paste(geneName,collapse=";"),scGene_exp=paste(scGene_exp,collapse=";"))
colnames(aCREanno3)[1]="peakID"
CREfinal=read.delim("/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/JC_merged/all_aCRE.TSS.access.markalone_cuttag.tsv", header=T, stringsAsFactors = F, check.names = F)
CREfinal=left_join(CREfinal,aCREanno3,by="peakID",copy=F)
CREfinal$gene_type2="others"
CREfinal$gene_type2[grep("ncRNA",CREfinal$gene_type)]="ncRNA"
CREfinal$gene_type2[grep("protein",CREfinal$gene_type)]="mRNA"
CREfinal$gene_type2[which(is.na(CREfinal$gene_type))]="no_gene"
write.table(unique(CREfinal[,c(1,42:48)]),gzfile(paste0(CRE_folder,"CRE_to_gene/filtered_peaks.collapsed gene.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#
aCREanno2=read.delim(paste0(CRE_folder,"CRE_to_gene/filtered_peaks.merged.all.table5gene.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
aCREanno2$peakID=gsub("-","_",aCREanno2$V4)

aCREanno2$rank=3
aCREanno2$rank[which(aCREanno2$gene_type == "protein_coding")]=1
aCREanno2$rank[grep("pseudogene",aCREanno2$gene_type)]=2
aCREanno2$rank[which(aCREanno2$gene_type %in% c("others"))]=5
aCREanno2$rank[which(aCREanno2$gene_type %in% c("lncRNA"))]=4
aCREanno2_p=aCREanno2%>%group_by(peakID)%>%dplyr::slice_min(rank)
aCREanno2_p=aCREanno2_p%>%group_by(peakID)%>%dplyr::slice_sample(n=1)
write.table(aCREanno2_p,gzfile(paste0(CRE_folder,"CRE_to_gene/filtered_peaks.merged.all.table5gene_rank_selected.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)

#update CRE file
scACREv3=read.delim(paste0(peak_folder,"all_aCRE.PE.final.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
scACREv3=left_join(scACREv3, aCREanno2_p[,c("V4","gene","geneName","simple_geneID")], by=c("peakID"="V4"), copy=F)
scACREv3=left_join(scACREv3,all_markers2, by=c("geneName"="gene"), copy=F)
write.table(scACREv3,gzfile(paste0(peak_folder,"all_aCRE.PE.final.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#update cut&tag file
CREfinal=read.delim(paste0(peak_folder,"all_aCRE.TSS.access.markalone_cuttag.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
CREfinal=left_join(CREfinal,aCREanno2_p[,c("V4","gene","geneName","simple_geneID","HGNC","exp_gene","scGene_exp","gene_type")], by=c("peakID"="V4"),copy=F)
write.table(CREfinal,gzfile(paste0(peak_folder,"all_aCRE.TSS.access.markalone_cuttag.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

> length(which(CREfinal$promoter_type_CT == "promoter-like"))
[1] 34435
> length(which(CREfinal$promoter_type_CT == "promoter-like" & CREfinal$tCRE == "Yes"))
[1] 22341
> length(which(CREfinal$promoter_type_CT == "promoter-like" & CREfinal$tCRE == "Yes" & !is.na(CREfinal$gene)))
[1] 20959
> length(which(CREfinal$promoter_type_CT == "enhancer-like" & CREfinal$tCRE == "Yes"))
[1] 20951

