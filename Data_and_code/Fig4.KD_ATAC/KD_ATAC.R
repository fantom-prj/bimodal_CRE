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
KDATAC_folder=paste0(primary_folder,"Data_and_code/Fig4.KD_ATAC/")
path_fig4_data=paste0(primary_folder,"Fig4/data/")

#===============================================================================
setwd(KDATAC_folder)

#mapping by BWA
system("sh mapping.sh")
#fastq and bam are maintained in DDBJ and FANTOM webpage

#===============================================================================
#counting by subread featureCounts

#prepare SAF format file for subread
ATAC=read.delim(paste0(primary_folder,"Data_and_code/scATAC/merged_peak/filtered_peaks.merged.all.bed.gz"), header=F, stringsAsFactors = F)
ATAC$V2=ATAC$V2+1
ATAC=ATAC[,c(4,1,2,3,6)]
colnames(ATAC)=c("GeneID", "Chr", "Start", "End", "Strand")
write.table(ATAC, "all_ATAC_peak.saf", col.names=T, row.names=F, sep="\t", quote=F)

system("sh featureCounts.sh")
system("gzip snATAC_peak_counts.tsv")
#output count file "snATAC_peak_counts.tsv.gz" located in KDATAC_folder

#===============================================================================
#bam count summary
setwd(KDATAC_folder)
count_sum=read.delim("snATAC_peak_counts.summary.tsv.gz", header=T, check.names = F, stringsAsFactors = F)
count_sum_t=data.frame(t(count_sum))
colnames(count_sum_t)=count_sum_t[1,]
count_sum_t=count_sum_t[-1,]
count_sum_t$lib=rownames(count_sum_t)
v1=ggplot()+
  geom_bar(data=count_sum_t, aes(x=lib, y=as.numeric(Assigned)), stat="identity")+
  coord_cartesian(ylim=c(0,70000000))+
  labs(y="read",x=NULL)+
  theme1
jpeg("snATAC_peak_counts.summary.jpeg", width = 3.5, height = 3.5, units = "in", res=600)
print(v1)
dev.off()

#=====================================================================
setwd(KDATAC_folder)
count=read.delim("snATAC_peak_counts.tsv.gz", skip=1, header=T, stringsAsFactors = F, check.names = F)
colnames(count)=gsub(".bam","",colnames(count))
colnames(count)=gsub("-","_",colnames(count))
rownames(count)=count$Geneid
count1=count[,c(7:18)]
info=read.delim(paste0(primary_folder,"Data_and_code/scATAC/merged_peak/all_aCRE.PE.final.tsv.gz"), header=T,  stringsAsFactors = F, check.names = F)
info=info[which(info$peakID %in% rownames(count)),]
study=data.frame(rbind(c("iPSC_scr", "iPSC_TEAD24_1", "iPSC_TEAD24_2"),c("Neuron_scr","Neuron_ONECUT2_1", "Neuron_ONECUT2_2")))
colnames(study)=c("Control","Query1","Query2")
row.names(study)=c("TEAD24KD","ONECUT2KD")

#===============================================================================
#edgeR with factor of GC content
library(edgeR)
library(EDASeq)
library(GenomicAlignments)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(wesanderson)
library(Hmisc)

txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
ff = FaFile("/analysisdata/fantom6/Interactome/resources/hg38.fa")
#fasta from https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/

library(EDASeq)
count1=as.matrix(count1)

gr = GRanges(seqnames=count$Chr, ranges=IRanges(count$Start, count$End), strand="*", mcols=data.frame(peakID=count$Geneid))
peakSeqs = getSeq(x=ff, gr)
gcContentPeaks = letterFrequency(peakSeqs, "GC",as.prob=TRUE)[,1]

#divide into 20 bins by GC content
gcGroups = Hmisc::cut2(gcContentPeaks, g=20)
mcols(gr)$gc = gcContentPeaks
dataOffset = withinLaneNormalization(count1,y=gcContentPeaks,num.bins=20,which="full",offset=TRUE)
dataOffset = betweenLaneNormalization(count1,which="full",offset=TRUE)

lowListGC = list()
for(kk in 1:ncol(count1)){
  set.seed(kk)
  lowListGC[[kk]] = lowess(x=gcContentPeaks, y=log1p(count1[,kk]), f=1/10)}

names(lowListGC)=colnames(count1)
dfList = list()
for(ss in 1:length(lowListGC)){
  oox = order(lowListGC[[ss]]$x)
  dfList[[ss]] = data.frame(x=lowListGC[[ss]]$x[oox], y=lowListGC[[ss]]$y[oox], sample=names(lowListGC)[[ss]])}
dfAll = do.call(rbind, dfList)
dfAll$sample = factor(dfAll$sample)

p1.1 = ggplot(dfAll, aes(x=x, y=y, group=sample, color=sample)) +
  geom_line(size = 1) +
  xlab("GC-content") +
  ylab("log(count + 1)") +
  theme_classic()
pdf("GCcontent_peaks.pdf")
p1.1
dev.off()

path5=paste0(KDATAC_folder,"edgeR/")
for (i in 1: 2){
  name=row.names(study)[i]
  TS=c(rep(study$Control[i],2),rep(study$Query1[i],2),rep(study$Query2[i],2))
  TS <- factor(TS, levels=c(study$Control[i],study$Query1[i],study$Query2[i]))
  counttable1=count1[,c(grep(study$Control[i],colnames(count1)),grep(study$Query1[i],colnames(count1)),grep(study$Query2[i],colnames(count1)))]
  counttable1=counttable1[which(rowSums(counttable1)>=100),]
  dataOffset1=dataOffset[,c(grep(study$Control[i],colnames(dataOffset)),grep(study$Query1[i],colnames(dataOffset)),grep(study$Query2[i],colnames(dataOffset)))]
  dataOffset1=dataOffset1[rownames(counttable1),]
  design = model.matrix(~0+TS)
  colnames(design)=study[i,]
  my_data1 = DGEList(counts=counttable1, group=TS)
  keep = filterByExpr(my_data1)
  summary(keep)
  my_data1=my_data1[keep,,keep.lib.sizes=FALSE]
  my_data1$offset = -dataOffset1[keep,]
  my_data1.eda = estimateGLMCommonDisp(my_data1, design = design)
  fit = glmQLFit(my_data1.eda, design = design)
  qlf.DEASeq = glmQLFTest(fit, contrast=c(-1,0.5,0.5))
  exp.df=as.data.frame(topTags(qlf.DEASeq, nrow(qlf.DEASeq$table)))
  exp.df$peakID=row.names(exp.df)
  exp.df = right_join(info,exp.df, by="peakID", copy=F)
  write.table(exp.df[order(exp.df$FDR),], file=gzfile(paste0(path5,name,".DE.tsv.gz")), col.names=TRUE, row.names=FALSE, sep="\t")}

#===============================================================================
#put together
#the number of CRE doenst match complete because we only consider the opened CRE in a cell type
DE1=read.delim(paste0(KDATAC_folder,"edgeR/TEAD24KD.DE.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
DE2=read.delim(paste0(KDATAC_folder,"edgeR/ONECUT2KD.DE.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
DE1=DE1[,c(1:11)]
DE1$sig="no"
DE1$sig[which(DE1$FDR<0.05 & DE1$logFC < (-0.5))]="yes"
DE2=DE2[,c(1:11)]
DE2$sig="no"
DE2$sig[which(DE2$FDR<0.05 & DE2$logFC < (-0.5))]="yes"

#===============================================================================
# load TEAD4 and ONECUT2 footprinting result
aCREe=read.delim(paste0(path_fig4_data,"peaks.merged.all.markalone_k16.subclass_TEAD_ONECUT.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
aCREe$TEAD_KD="No"
aCREe$TEAD_KD[which(aCREe$V4 %in% DE1$peakID[which(DE1$sig == "yes")])]="Yes"
aCREe$ONECUT_KD="No"
aCREe$ONECUT_KD[which(aCREe$V4 %in% DE2$peakID[which(DE2$sig == "yes")])]="Yes"

aCREe$TEAD_Motif_ChIP="No"
aCREe$TEAD_Motif_ChIP[which(aCREe$TEAD_motif == "Yes" & aCREe$TEAD_ChIP == "Yes")]="Yes"
aCREe$TEAD_FootPrint_ChIP="No"
aCREe$TEAD_FootPrint_ChIP[which(aCREe$TEAD_FootPrint == "Yes" & aCREe$TEAD_ChIP == "Yes")]="Yes"
aCREe$TEAD_FootPrint_KD="No"
aCREe$TEAD_FootPrint_KD[which(aCREe$TEAD_FootPrint == "Yes" & aCREe$TEAD_KD == "Yes")]="Yes"
aCREe$TEAD_FootPrint_ChIP_KD="No"
aCREe$TEAD_FootPrint_ChIP_KD[which(aCREe$TEAD_FootPrint == "Yes" & aCREe$TEAD_KD == "Yes" & aCREe$TEAD_ChIP == "Yes")]="Yes"
aCREe$ONECUT_Motif_ChIP="No"
aCREe$ONECUT_Motif_ChIP[which(aCREe$ONECUT_motif == "Yes" & aCREe$ONECUT_ChIP == "Yes")]="Yes"
aCREe$ONECUT_FootPrint_ChIP="No"
aCREe$ONECUT_FootPrint_ChIP[which(aCREe$ONECUT_FootPrint == "Yes" & aCREe$ONECUT_ChIP == "Yes")]="Yes"
aCREe$ONECUT_FootPrint_KD="No"
aCREe$ONECUT_FootPrint_KD[which(aCREe$ONECUT_FootPrint == "Yes" & aCREe$ONECUT_KD == "Yes")]="Yes"
aCREe$ONECUT_FootPrint_ChIP_KD="No"
aCREe$ONECUT_FootPrint_ChIP_KD[which(aCREe$ONECUT_FootPrint == "Yes" & aCREe$ONECUT_KD == "Yes" & aCREe$ONECUT_ChIP == "Yes")]="Yes"

write.table(aCREe,gzfile(paste0(path_fig4_data,"peaks.merged.all.markalone_k16.subclass_TEAD_ONECUT.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)

#==
DE1$Motif="no"
DE1$Motif[which(DE1$peakID %in% aCREe$V4[which(aCREe$TEAD_motif == "Yes")])]="yes"
DE1$ChIP="no"
DE1$ChIP[which(DE1$peakID %in% aCREe$V4[which(aCREe$TEAD_ChIP == "Yes")])]="yes"
DE1$Footprint="no"
DE1$Footprint[which(DE1$peakID %in% aCREe$V4[which(aCREe$TEAD_FootPrint == "Yes")])]="yes"
DE1$sigup="no"
DE1$sigup[which(DE1$FDR<0.05 & DE1$logFC > (0.5))]="yes"

aa=DE1%>%group_by(sig, Motif)%>%dplyr::summarise(count=n())
bb=DE1%>%group_by(sig, ChIP)%>%dplyr::summarise(count=n())
cc=DE1%>%group_by(sig, Footprint)%>%dplyr::summarise(count=n())
aa1=DE1%>%group_by(sigup, Motif)%>%dplyr::summarise(count=n())
bb1=DE1%>%group_by(sigup, ChIP)%>%dplyr::summarise(count=n())
cc1=DE1%>%group_by(sigup, Footprint)%>%dplyr::summarise(count=n())

tead_summary=data.frame(rbind(aa$count,bb$count,cc$count,aa1$count,c(bb1$count,0),cc1$count))
colnames(tead_summary)=c("nono","noyes","yesno","yesyes")
tead_summary$feature1=c(rep("KD_sig_down",3),rep("KD_sig_up",3))
tead_summary$feature2=c("motif","ChIP","Footprint")
for(i in 1:6){
  GSEATasting <- matrix(c(tead_summary$yesyes[i], tead_summary$yesno[i], tead_summary$noyes[i], tead_summary$nono[i]), nrow = 2, dimnames = list(oligo1 = c("yes", "no"), oligo2 = c("yes", "no")))
  tead_summary$p.val[i] = fisher.test(GSEATasting, alternative = "two.sided")$p.value
  tead_summary$OR[i] = fisher.test(GSEATasting, alternative = "two.sided")$estimate}
tead_summary$logOdds=log(tead_summary$OR)

DE2$Motif="no"
DE2$Motif[which(DE2$peakID %in% aCREe$V4[which(aCREe$ONECUT_motif == "Yes")])]="yes"
DE2$ChIP="no"
DE2$ChIP[which(DE2$peakID %in% aCREe$V4[which(aCREe$ONECUT_ChIP == "Yes")])]="yes"
DE2$Footprint="no"
DE2$Footprint[which(DE2$peakID %in% aCREe$V4[which(aCREe$ONECUT_FootPrint == "Yes")])]="yes"
DE2$sigup="no"
DE2$sigup[which(DE2$FDR<0.05 & DE2$logFC > (0.5))]="yes"

aa=DE2%>%group_by(sig, Motif)%>%dplyr::summarise(count=n())
bb=DE2%>%group_by(sig, ChIP)%>%dplyr::summarise(count=n())
cc=DE2%>%group_by(sig, Footprint)%>%dplyr::summarise(count=n())
aa1=DE2%>%group_by(sigup, Motif)%>%dplyr::summarise(count=n())
bb1=DE2%>%group_by(sigup, ChIP)%>%dplyr::summarise(count=n())
cc1=DE2%>%group_by(sigup, Footprint)%>%dplyr::summarise(count=n())

onecut_summary=data.frame(rbind(aa$count,bb$count,cc$count,aa1$count,c(bb1$count,0),cc1$count))
colnames(onecut_summary)=c("nono","noyes","yesno","yesyes")
onecut_summary$feature1=c(rep("KD_sig_down",3),rep("KD_sig_up",3))
onecut_summary$feature2=c("motif","ChIP","Footprint")
for(i in 1:6){
  GSEATasting <- matrix(c(onecut_summary$yesyes[i], onecut_summary$yesno[i], onecut_summary$noyes[i], onecut_summary$nono[i]), nrow = 2, dimnames = list(oligo1 = c("yes", "no"), oligo2 = c("yes", "no")))
  onecut_summary$p.val[i] = fisher.test(GSEATasting, alternative = "two.sided")$p.value
  onecut_summary$OR[i] = fisher.test(GSEATasting, alternative = "two.sided")$estimate}
onecut_summary$logOdds=log(onecut_summary$OR)

write.table(tead_summary, gzfile(paste0(path_fig4_data,"TEAD24KD.DE.summary.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(onecut_summary, gzfile(paste0(path_fig4_data,"ONECUT2KD.DE.summary.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
# -> for fig 4b

aCRE_cell=read.delim(paste0(peak_folder,"all_aCRE.TSS.access.markalone_cuttag.tsv.gz"), header=T, stringsAsFactors = F)
DE1=left_join(DE1, aCRE_cell[,c(1,8:14)], by="peakID", copy=F)
DE2=left_join(DE2, aCRE_cell[,c(1,8:14)], by="peakID", copy=F)

write.table(DE1, gzfile(paste0(KDATAC_folder,"edgeR/TEAD24KD.DE.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(DE2, gzfile(paste0(KDATAC_folder,"edgeR/ONECUT2KD.DE.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(DE1[which(DE1$sig == "yes" & DE1$Footprint == "yes"),], gzfile(paste0(KDATAC_folder,"edgeR/TEAD24KD.DE.final.list.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(DE2[which(DE2$sig == "yes" & DE2$Footprint == "yes"),], gzfile(paste0(KDATAC_folder,"edgeR/ONECUT2KD.DE.final.list.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#table for plotting
write.table(DE1, gzfile(paste0(path_fig4_data,"TEAD24KD.DE.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(DE2, gzfile(paste0(path_fig4_data,"ONECUT2KD.DE.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
# prepare table for number of features
aCREe=read.delim(paste0(path_fig4_data,"peaks.merged.all.markalone_k16.subclass_TEAD_ONECUT.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
maCREe1=reshape2::melt(aCREe[,c(1,3,9,grep("TEAD",colnames(aCREe)))],id=c(1:3))
maCREe2=reshape2::melt(aCREe[,c(1,7,9,grep("ONECUT",colnames(aCREe)))],id=c(1:3))
colnames(maCREe2)=colnames(maCREe1)
maCREe=rbind(maCREe1,maCREe2)
maCREe=maCREe[which(maCREe$value == "Yes"),]
maCREe$k16ips_all[which(is.na(maCREe$k16ips_all))]="Others"
maCREe$k16ips_all[which(maCREe$k16ips_all=="CTCF")]="Others"
maCREe$group=sapply(strsplit(as.character(maCREe$variable),"_"),"[",1)
maCREe$group=gsub("TEAD","TEAD in iPSC",maCREe$group)
maCREe$group=gsub("ONECUT","ONECUT in Neuron",maCREe$group)
maCREe$variable=gsub("TEAD_","",as.character(maCREe$variable))
maCREe$variable=gsub("ONECUT_","",as.character(maCREe$variable))
colnames(maCREe)[2]="chromatinState"
write.table(maCREe,gzfile(paste0(path_fig4_data,"proportion_chromatin_state.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)
# -> for fig 4c

#===============================================================================
# prepare table for motif activity of TEAD and ONECUT
combined=readRDS(paste0(primary_folder,"Data_and_code/DATA/metacell.combined.feature.Rds"))
# only included CRE that i needed
acre=unique(maCREe$V4[which(maCREe$variable %in% c("FootPrint_ChIP","FootPrint_KD"))])
need.atac=which(rownames(combined[["aCRE"]]) %in% acre)
need.rna=which(rownames(combined[["tCRE"]]) %in% acre)
combined[["aCRE"]] <- subset(combined[["aCRE"]], features = rownames(combined[["aCRE"]])[need.atac])
combined[["tCRE"]] <- subset(combined[["tCRE"]], features = rownames(combined[["tCRE"]])[need.rna])
zt=as.data.frame(combined@assays[["tCRE"]]@data)
za=as.data.frame(combined@assays[["aCRE"]]@data)
rm(combined)
zt$peakID=rownames(zt)
za$peakID=rownames(za)
zt = left_join(zt, aCREe[,c(1,18:23)], by=c("peakID"="V4"),copy=F)
za = left_join(za, aCREe[,c(1,18:23)], by=c("peakID"="V4"),copy=F)
write.table(zt,gzfile(paste0(path_fig4_data,"tCRE.matrix.391.TEAD.ONECUT.ChIP_KD.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(za,gzfile(paste0(path_fig4_data,"aCRE.matrix.391.TEAD.ONECUT.ChIP_KD.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===========
pseudotime1=read.delim(paste0(primary_folder,"Data_and_code/Fig1.metacell/output/pseudotime.scRNA.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
pseudotime1=pseudotime1[which(pseudotime1$SEACell %in% colnames(zt)[c(1:391)]),]
pseudotime1$cluster=sapply(strsplit(pseudotime1$SEACell,"_meta"),"[",1)
cluster_info=read.delim(paste0(primary_folder,"Fig1/cluster_info.tsv"),header=T, stringsAsFactors = F, check.names = F)
pseudotime1=left_join(pseudotime1,cluster_info[,c(1,3)],by="cluster",copy=F)

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

pseudotime4=left_join(pseudotime1[,c(1,2)],pseudotime2[,c(1,7)],by="SEACell")
pseudotime4=left_join(pseudotime4,pseudotime3[,c(1,7)],by="SEACell",suffix=c("_1","_2"))

#for TEAD in iPSC
zt=read.delim("tCRE.matrix.391.TEAD.ONECUT.ChIP_KD.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
za=read.delim("aCRE.matrix.391.TEAD.ONECUT.ChIP_KD.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
zt=left_join(zt, aCREe[,c(1,3)], by=c("peakID"="V4"),copy=F)
za=left_join(za, aCREe[,c(1,3)], by=c("peakID"="V4"),copy=F)
TEADt=reshape2::melt(zt[,c(1:395,399)], id=c(1:392,396))
TEADt=TEADt[which(TEADt$value == "Yes"),]
TEADt1=data.frame(TEADt%>%group_by(variable,k16ips_all)%>%dplyr::summarise(across(1:391,mean)))
TEADt2=reshape2::melt(TEADt1, id=c(1,2))
TEADa=reshape2::melt(za[,c(1:395,399)], id=c(1:392,396))
TEADa=TEADa[which(TEADa$value == "Yes"),]
TEADa1=data.frame(TEADa%>%group_by(variable,k16ips_all)%>%dplyr::summarise(across(1:391,mean)))
TEADa2=reshape2::melt(TEADa1, id=c(1,2))

TEADt2$bimodal="RNA"
TEADa2$bimodal="ATAC"
TEAD=rbind(TEADt2,TEADa2)
colnames(TEAD)=c("group","chromatin_state","metacell","value","bimodal")

TEAD=left_join(TEAD,pseudotime4[,c(1,3,4)], by=c("metacell"="SEACell"), copy=F)
TEAD$chromatin_state[grep("romoter",TEAD$chromatin_state)]="All_promoter"
TEAD1=TEAD%>%group_by(bimodal,chromatin_state,group,T2_1)%>%dplyr::summarise(value=mean(value))
TEAD1$trajectory="Trajectory 1"
colnames(TEAD1)[4]="pseudotime"
TEAD2=TEAD%>%group_by(bimodal,chromatin_state,group,T2_2)%>%dplyr::summarise(value=mean(value))
TEAD2$trajectory="Trajectory 2"
colnames(TEAD2)[4]="pseudotime"
TEAD=rbind(TEAD1,TEAD2)
TEAD=TEAD[which(!is.na(TEAD$pseudotime)),]
TEAD=TEAD%>%group_by(trajectory,bimodal,chromatin_state)%>%dplyr::mutate(scale_value=value/max(value))
TEAD$group=gsub("TEAD_FootPrint_","",as.character(TEAD$group))
TEAD$group=factor(TEAD$group, levels=c("ChIP","KD","ChIP_KD"))
TEAD$chromatin_state=factor(TEAD$chromatin_state, levels=c("Active_enhancer","Primed_enhancer","All_promoter"))

#==
# ONECUT in Neuron
zt=read.delim(paste0(path_fig4_data,"tCRE.matrix.391.TEAD.ONECUT.ChIP_KD.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
za=read.delim(paste0(path_fig4_data,"aCRE.matrix.391.TEAD.ONECUT.ChIP_KD.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
zt=left_join(zt, aCREe[,c("V4","k16nrn_all")], by=c("peakID"="V4"),copy=F)
za=left_join(za, aCREe[,c("V4","k16nrn_all")], by=c("peakID"="V4"),copy=F)

ONECUTt=reshape2::melt(zt[,c(1:392,396:399)], id=c(1:392,396))
ONECUTt=ONECUTt[which(ONECUTt$value == "Yes"),]
ONECUTt1=data.frame(ONECUTt%>%group_by(variable,k16nrn_all)%>%dplyr::summarise(across(1:391,mean)))
ONECUTt2=reshape2::melt(ONECUTt1, id=c(1,2))
ONECUTa=reshape2::melt(za[,c(1:392,396:399)], id=c(1:392,396))
ONECUTa=ONECUTa[which(ONECUTa$value == "Yes"),]
ONECUTa1=data.frame(ONECUTa%>%group_by(variable,k16nrn_all)%>%dplyr::summarise(across(1:391,mean)))
ONECUTa2=reshape2::melt(ONECUTa1, id=c(1,2))

ONECUTt2$bimodal="RNA"
ONECUTa2$bimodal="ATAC"
ONECUT=rbind(ONECUTt2,ONECUTa2)
colnames(ONECUT)=c("group","chromatin_state","metacell","value","bimodal")

ONECUT=left_join(ONECUT,pseudotime4[,c(1,3,4)], by=c("metacell"="SEACell"), copy=F)
ONECUT$chromatin_state[grep("romoter",ONECUT$chromatin_state)]="All_promoter"
ONECUT1=ONECUT%>%group_by(bimodal,chromatin_state, group,T2_1)%>%dplyr::summarise(value=mean(value))
ONECUT1$trajectory="Trajectory 1"
colnames(ONECUT1)[4]="pseudotime"
ONECUT2=ONECUT%>%group_by(bimodal,chromatin_state, group,T2_2)%>%dplyr::summarise(value=mean(value))
ONECUT2$trajectory="Trajectory 2"
colnames(ONECUT2)[4]="pseudotime"
ONECUT=rbind(ONECUT1,ONECUT2)
ONECUT=ONECUT[which(!is.na(ONECUT$pseudotime)),]
ONECUT=ONECUT%>%group_by(trajectory,bimodal,chromatin_state)%>%dplyr::mutate(scale_value=value/max(value))
ONECUT$group=gsub("ONECUT_FootPrint_","",as.character(ONECUT$group))
ONECUT$group=factor(ONECUT$group, levels=c("ChIP","KD","ChIP_KD"))
ONECUT$chromatin_state=factor(ONECUT$chromatin_state, levels=c("Active_enhancer","Primed_enhancer","All_promoter"))

TEAD$TF="TEAD"
ONECUT$TF="ONECUT"
both=rbind(TEAD, ONECUT)
write.table(both,gzfile(paste0(path_fig4_data,"TEAD_ONECUT.neuronal.trajectory.chromatin_state.aCRE.vs.tCRE.activity.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)
# -> for fig 4d

#==============
