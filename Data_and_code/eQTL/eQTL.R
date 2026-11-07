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

###color#####
library("ggsci")
mypal = pal_npg("nrc", alpha = 1)(9)
mypal
library("scales")
show_col(mypal)

#####################
primary_folder="/analysisdata/fantom6/Interactome/single_cell_wallace/Figure/"
eQTL_folder=paste0(primary_folder,"Data_and_code/eQTL/")
peak_folder=paste0(primary_folder,"Data_and_code/scATAC/merged_peak/")

#===================
setwd(eQTL_folder)
#data from https://www.ebi.ac.uk/eqtl/

metadata=read.delim("tabix_ftp_paths.tsv", header=T, stringsAsFactors = F, check.names = F)
metadata1=metadata[union(grep("brain",metadata$tissue_label),grep("iPSC",metadata$tissue_label)),]

path3=("/analysisdata/fantom6/Interactome/single_cell_wallace/eqtl/download/")
files=list.files(path=path3, pattern="tsv.gz")
metadata2=data.frame(files)
metadata2$dataset_id=sapply(strsplit(metadata2$files,"\\."),"[",5)
metadata2$quant_method=sapply(strsplit(metadata2$files,"\\."),"[",3)
metadata2=left_join(metadata2, metadata, by=c("dataset_id","quant_method"), copy=F)
metadata2a=metadata2
metadata2=metadata2[which(metadata2$dataset_id %in% metadata1$dataset_id),]

data=read.delim(paste0(path3,metadata2$files[1]), header=T, stringsAsFactors = F, check.names = F)
data=data[which(data$pip>0.3),]
data$dataset_id=metadata2$dataset_id[1]
data$quant_method=metadata2$quant_method[1]
for (i in 2: nrow(metadata2)){
  data1=read.delim(paste0(path3,metadata2$files[i]), header=T, stringsAsFactors = F, check.names = F)
  data1=data1[which(data1$pip>0.3),]
  data1$dataset_id=metadata2$dataset_id[i]
  data1$quant_method=metadata2$quant_method[i]
  data=rbind(data,data1)}
data$chr=sapply(strsplit(data$variant,"_"),"[",1)
data$end=sapply(strsplit(data$variant,"_"),"[",2)
data$start=as.numeric(data$end)-1
data$ID=paste0(data$chr,"_",data$start,"_",data$end)
data=left_join(data, metadata2[,c(2,6)], by="dataset_id", copy=F)
datak=data%>%group_by(ID,variant, pip, gene_id,sample_group,quant_method)%>%dplyr::summarise(rsid=paste(rsid,collapse=";"))
datab=unique(datak[,c(1,4)])%>%group_by(ID)%>%dplyr::summarise(no_gene=n())
datac=unique(datak[,c(1,5)])%>%group_by(ID)%>%dplyr::summarise(no_sample=n())
dataa=unique(datak[,c(1,2)])%>%group_by(ID)%>%dplyr::summarise(no_variant=n())
datad=unique(datak[,c(1,6)])%>%group_by(ID)%>%dplyr::summarise(no_method=n())
data2=datak%>%group_by(ID)%>%dplyr::summarise(variant=paste(variant, collapse=";"),pip=paste(pip, collapse=";"),gene_id=paste(gene_id,collapse=";"),sample_group=paste(sample_group, collapse=";"),quant_method=paste(quant_method, collapse=";"))
data2=left_join(data2, dataa,by="ID",copy=F)
data2=left_join(data2, datab,by="ID",copy=F)
data2=left_join(data2, datac,by="ID",copy=F)
data2=left_join(data2, datad,by="ID",copy=F)
data3=unique(data[,c(16,18,17,19)])

write.table(data, gzfile("eqtl_LD_pip3.tsv.gz"), col.names=T, row.names=F, quote=F, sep="\t")
write.table(data2, gzfile("eqtl_LD_pip3.collapsed.tsv.gz"), col.names=T, row.names=F, quote=F, sep="\t")
write.table(data3, gzfile("eqtl_LD_pip3.collapsed.bed.gz"), col.names=F, row.names=F, quote=F, sep="\t")

#===================================================================
#bedtools intersect with aCRE
setwd(eQTL_folder)
system(paste0("bedtools intersect -wa -wb -a ",peak_folder,"filtered_peaks.merged.all.bed.gz -b eqtl_LD_pip3.collapsed.bed.gz | gzip > aCRE.eqtl_LD_pip3.collapsed.bed.gz"))
#===========================

data2=read.delim("eqtl_LD_pip3.collapsed.tsv.gz",header=T, check.names=F, stringsAsFactors = F)
aCRE=read.delim("aCRE.eqtl_LD_pip3.collapsed.bed.gz", header=F, stringsAsFactors = F)
CREfinal=read.delim("/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/JC_merged/all_aCRE.TSS.access.markalone_cuttag.tsv", header=T, stringsAsFactors = F, check.names = F)
aCRE=left_join(aCRE, CREfinal[,c(1,8:11)], by=c("V4"="peakID"))
aCRE=left_join(aCRE, data2[,c(1:6)], by=c("V10"="ID"),copy=F)
aCRE=aCRE[,c(4,11:14,10,15:19)]
colnames(aCRE)[c(1,7,9)]=c("peakID","variantID","eqtl_gene_id")
library(tidyr)
aCRE1=separate_rows(aCRE, variantID, pip, eqtl_gene_id, sample_group,quant_method, sep=";", convert = TRUE)
write.table(aCRE1,gzfile("aCRE.eqtl_LD_pip3_all_method.separate.bed.gz"), col.names=T, row.names=F, sep="\t", quote=F)
aCRE1$pair=paste0(aCRE1$peakID,":",aCRE1$eqtl_gene_id)
aCRE2=aCRE1%>%group_by(pair, peakID, iPSC_CT, NSC_CT, Neuron_CT, promoter_type2, eqtl_gene_id)%>%dplyr::summarise(V10=paste(V10, collapse=";"), variantID=paste(variantID,collapse=";"), pip=paste(pip, collapse=";"), sample_group=paste(sample_group, collapse=";"), quant_method=paste(quant_method, collapse=";"))

length(unique(aCRE$peakID))#15209
length(unique(aCRE$variantID))#18840
length(unique(aCRE$peakID[which(aCRE$promoter_type2 == "enhancer-like")]))#9700

length(unique(aCRE2$pair[which(aCRE2$promoter_type2 == "enhancer-like")]))#10796
write.table(aCRE,gzfile("aCRE.eqtl_LD_pip3_all_method.collapsed.bed.gz"), col.names=T, row.names=F, sep="\t", quote=F)
write.table(aCRE2,"aCRE2egene_eqtl_LD_pip3_all_method_pairbase.tsv", col.names=T, row.names=F, sep="\t", quote=F)

