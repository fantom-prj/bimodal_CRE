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
TF_out=paste0(primary_folder,"Data_and_code/Fig1.metacell/output/TF_motif/")
path_SEACell=paste0(primary_folder,"Data_and_code/Fig1.coembed/SEACell/")
bigDATA=paste0(primary_folder,"Data_and_code/DATA/")
path_fig1_data=paste0(primary_folder,"Fig1/data/")
path_fig3_data=paste0(primary_folder,"Fig3/data/")

#===============================================================================
# manually curate TF function
# extract TF gene name from motif
setwd(TF_out)
motif.df=read.delim(paste0(primary_folder,"Data_and_code/Fig3.ChromVar/jaspar2024/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.ID.tsv"), header=F, stringsAsFactors = F)

motif.df1=motif.df
motif.df1$motifName1=sapply(strsplit(motif.df1$V2,"::"),"[",1)
motif.df1$motifName2=sapply(strsplit(motif.df1$V2,"::"),"[",2)
motif.df2=reshape2::melt(motif.df1, id=c(1,2))
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
# removed 3 redundant entry
function2=read.delim("idmapping_2024_04_25.tsv",header=T, stringsAsFactors = F)
function2$TF_activity2=""
function2$TF_activity2[grep("activate",function2$Function..CC.)]="activator"
function2$TF_activity2[grep("activated",function2$Function..CC.)]=""
function2$TF_activity2[grep("activator",function2$Function..CC.)]="activator"
#function2$TF_activity[which(function2$curated=="no")]=function2$TF_activity2[which(function2$curated=="no")]
function2$TF_activity[which(function2$TF_activity == "unknown")]=function2$TF_activity2[which(function2$TF_activity == "unknown")]
function2$TF_activity[which(function2$TF_activity == "unknown" | function2$TF_activity == "uknown" | function2$TF_activity=="")]="unknown"

function2%>%group_by(TF_activity)%>%dplyr::summarise(count=n())
write.table(function2, "idmapping_2024_04_25.tsv", col.names=T, row.names=F, sep="\t", quote=F)
function4=function2%>%group_by(From)%>%dplyr::filter(!(TF_activity == "unknown" & any(TF_activity != "unknown"))) %>% unique()
motif.df2=left_join(motif.df2, function4, by=c("Approved.symbol"="From"), copy=F)


combined=readRDS(paste0(bigDATA,"metacell.combined.feature.Rds"))
expression=as.data.frame(t(LayerData(combined, assay = "RNA", layer = "data")))
expression1=expression[,which(colnames(expression)%in% unique(motif.df2$Approved.symbol))]
expression1$metacell=rownames(expression1)
expression2=reshape2::melt(expression1, id=747)
write.table(expression2,gzfile("TF.expression.per.meta.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)
expression3=expression2%>%group_by(variable)%>%dplyr::summarise(max=max(value),sum=sum(value),mean=mean(value),count=sum(value!=0))
 
#==
motif.df4=left_join(motif.df2, expression3, by=c("Approved.symbol"="variable"), copy=F)
write.table(motif.df4,gzfile("motif.info.with_TF.expression.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)

# motif.df4rep=motif.df4[grep("repress",motif.df4$TF_activity),]
# motif.df4rep=motif.df4rep[which(motif.df4rep$max >=0.2),]
# motif.df4rep=motif.df4rep%>%group_by(V1)%>%dplyr::mutate(nTF=n())
# motifout=motif.df4rep$V1[intersect(grep(":",motif.df4rep$V2), which(motif.df4rep$nTF ==1))]
# motif.df4rep=motif.df4rep[-which(motif.df4rep$V1 %in% motifout),]
# length(unique(motif.df4rep$V1)) #145
# write.table(motif.df4rep,"motif.info.with_TF.expression.repress.max02.tsv", col.names=T, row.names=F, sep="\t", quote=F)
# #==

motif.df1=left_join(motif.df1,unique(motif.df2[,c(4,6)]), by=c("motifName1"="value"))
motif.df1=left_join(motif.df1,unique(motif.df2[,c(4,6)]), by=c("motifName2"="value"), suffix=c(".TF1",".TF2"))

motif.df5=left_join(motif.df1, expression2, by=c("Approved.symbol.TF1"="variable"), copy=F)
motif.df5=left_join(motif.df5, expression2, by=c("Approved.symbol.TF2"="variable", "metacell"="metacell"), copy=F, suffix=c(".TF1",".TF2"))
colnames(motif.df5)[c(1:2)]=c("motifID","motifName")
write.table(motif.df5,gzfile("motif.base.TF.expression.per.meta.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)
# -> expression of TFs from all 879 motifs

#===============================================================================
#filter motif for expression, max 0.2 for 3 metacells --> 420 motif
motif.df5=read.delim("motif.base.TF.expression.per.meta.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
motif.df6a=motif.df5[which(is.na(motif.df5$motifName2)),]
motif.df6b=motif.df5[which(!is.na(motif.df5$motifName2)),]
motif.df6a_summary=motif.df6a[which(motif.df6a$value.TF1 >= 0.2),]%>%group_by(motifID)%>%dplyr::summarise(count=n())
need1=motif.df6a_summary$motifID[which(motif.df6a_summary$count>=3)] #394
motif.df6b_summary=motif.df6b[which(motif.df6b$value.TF1 >= 0.2 & motif.df6b$value.TF2 >=0.2),]%>%group_by(motifID)%>%dplyr::summarise(count=n())
need2=motif.df6b_summary$motifID[which(motif.df6b_summary$count>=3)] #26
length(union(need1,need2)) #420
motif.df6=motif.df5[which(motif.df5$motifID %in% union(need1,need2)),]
write.table(motif.df6,gzfile("motif.base.TF.expression.per.max02.cell3.meta.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
#prepare literature review on known TF function, for FE later

motif.df6=read.delim(paste0(metacell,"output/TF_motif/motif.base.TF.expression.per.max02.cell3.meta.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
TF_involve=union(unique(motif.df6$Approved.symbol.TF1), unique(motif.df6$Approved.symbol.TF2)) #356 TF, 400 motif
motif.df6d=unique(motif.df6[,c(1:6)])
motif.df6e1=unique(motif.df6[,c(5,7,8)])
motif.df6e2=unique(motif.df6[which(!is.na(motif.df6$motifName2)),c(6,7,9)])
colnames(motif.df6e2)=colnames(motif.df6e1)
motif.df6e=unique(rbind(motif.df6e1,motif.df6e2))
motif.df6f=spread(motif.df6e, key=1, value=3)
write.table(motif.df6d,gzfile(paste0(path_fig3_data,"expressed_420motif.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(motif.df6f,gzfile(paste0(path_fig3_data,"expressed_356TFx393meta.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

# motif.df7=read.delim(paste0(metacell,"output/TF_motif/motif.info.with_TF.expression.repress.max02.tsv"),header=T, stringsAsFactors = F, check.names = F)
# motif.df7$label=paste0(motif.df7$V1,"_",motif.df7$V2)
# motif.df7$label2=gsub("motifName","TF",motif.df7$variable)
# motif.df7$label[grep("::",motif.df7$label)]=paste0(motif.df7$label[grep("::",motif.df7$label)],"_",motif.df7$label2[grep("::",motif.df7$label)])

close=c("REST","EZH2","ZNF217","KLF3","BACH1","FOXG1","MAFG","TGIF2","INSM1")

motif.df8=read.delim(paste0(metacell,"output/TF_motif/motif.info.with_TF.expression.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
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
motif.df8$TF_pioneer[grep("yes",motif.df8$`pioneer.factor`)]="Yes"
motif.df8$TF_pioneer[grep("maybe",motif.df8$`pioneer.factor`)]="Yes"
motif.df8$TF_closer="No"
motif.df8$TF_closer[which(motif.df8$Approved.symbol %in% close)]="Yes"

motif.df8%>%group_by(TF_activator,TF_repressor,TF_pioneer,TF_closer)%>%dplyr::summarise(count=n())

#TF to motif base
motif.df8a=unique(motif.df8[,c(1,2,30:33)])%>%group_by(V1,V2)%>%dplyr::summarise(
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
write.table(motif.df8b,gzfile(paste0(path_fig3_data,"motif.info.with_TF.expression_n420.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
