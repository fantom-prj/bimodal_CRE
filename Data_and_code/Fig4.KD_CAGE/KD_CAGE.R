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
KDCAGE_folder=paste0(primary_folder,"Data_and_code/Fig4.KD_CAGE/")
path_fig4_data=paste0(primary_folder,"Fig4/data/")

#===============================================================================
setwd(KDCAGE_folder)
files=list.files(path=paste0(KDCAGE_folder,"bam_to_ctss"), pattern="unencoded_G.collapse.ctss.bed.gz", recursive = T)
files=files[-grep("tbi", files)]
files.names=sapply(strsplit(files, "\\/"), "[", 3)
files.names=gsub(".unencoded_G.collapse.ctss.bed.gz","",files.names)
summary=data.frame(matrix(nrow=12, ncol=3))
colnames(summary)=c("libID","CTSS_count","unG_count")

for (i in 1:6){
  data=read.delim(paste0(KDCAGE_folder,"bam_to_ctss/",files[i]), header=F, stringsAsFactors = F)
  summary$libID[i]=files.names[i]
  summary$CTSS_count[i]= sum(data$V4)
  summary$unG_count[i]=sum(data$V5)}

# intersect with sc5nRNA-seq SCAFE defined TSS cluster, and then to ATAC peak defined CRE
setwd(KDCAGE_folder)
files=list.files(path=paste0(KDCAGE_folder,"bam_to_ctss"), pattern="unencoded_G.collapse.ctss.bed.gz", recursive = T)

for (i in 1:6){
  system(paste0("bedtools intersect -wa -a ",KDCAGE_folder,"bam_to_ctss/",files[i]," -b ",primary_folder,"Data_and_code/SCAFE/sc5end/aggregate/run_full/out/cluster/sc5end.iPSC_NSC_Neuron/bed/sc5end.iPSC_NSC_Neuron.tssCluster.bed.gz -s | gzip > ./CTSS_in_cluster/", files.names[i],".TSScluster.bed.gz"))}

files=list.files(path=paste0(KDCAGE_folder,"CTSS_in_cluster"), pattern="TSScluster.bed.gz", recursive = T)
files.names=gsub(".TSScluster.bed.gz","",files)

for (i in 1:6){
  system(paste0("bedtools intersect -wa -wb -a ",KDCAGE_folder,"CTSS_in_cluster/",files[i]," -b ",primary_folder,"Data_and_code/scATAC/merged_peak/filtered_peaks.merged.all.bed.gz | gzip > ./CTSS_in_aCRE/", files.names[i],".aCRE.bed.gz"))}

#===============================================================================
# collapse into CRE count
files=list.files(path=paste0(KDCAGE_folder,"CTSS_in_aCRE"),pattern=".aCRE.bed.gz", recursive = T)
files.names=gsub(".aCRE.bed.gz","",files)

summary1=data.frame(matrix(nrow=12, ncol=3))
colnames(summary1)=c("libID","CTSS_count_in_aCRE","unG_count_in_aCRE")
#aCRE=read.delim(paste0(primary_folder,"Data_and_code/scATAC/merged_peak/filtered_peaks.merged.all.bed.gz"), header=F, stringsAsFactors = F)
data0=data.frame()

for (i in 1:6){
  data=read.delim(paste0(KDCAGE_folder,"CTSS_in_aCRE/",files[i]), header=F, stringsAsFactors = F)
  data1=data%>%group_by(V10)%>%dplyr::summarise(CTSS=sum(V4), unG=sum(V5), libID=files.names[i])
  summary1$libID[i]=files.names[i]
  summary1$CTSS_count_in_aCRE[i]= sum(data1$CTSS)
  summary1$unG_count_in_aCRE[i]=sum(data1$unG)
  data0=rbind(data0, data1)}

summary=left_join(summary, summary1, by="libID", copy=F)

CRE_count=spread(data0[,c(1,2,4)], key=3, value=2)
CRE_count[is.na(CRE_count)]=0
colnames(CRE_count)[1]="peakID"
colnames(CRE_count)[2:7]=c("iPSC_scr_r1","iPSC_scr_r2", "iPSC_TEAD24_1_r1","iPSC_TEAD24_1_r2", "iPSC_TEAD24_2_r1", "iPSC_TEAD24_2_r2")

write.table(CRE_count[,c(1:7)],gzfile("aCRE_count_TEAD24KD.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
# edgeR
setwd(KDCAGE_folder)
CRE_count=read.delim("aCRE_count_TEAD24KD.tsv.gz", header=T, row.names=1)

library(edgeR)
TS=c(rep("iPSC_scr",2),rep("iPSC_TEAD24_1",2),rep("iPSC_TEAD24_2",2))
TS <- factor(TS, levels=c("iPSC_scr", "iPSC_TEAD24_1", "iPSC_TEAD24_2"))
CRE_count=CRE_count[which(rowSums(CRE_count)>=0),]
design = model.matrix(~0+TS)
colnames(design)=c("iPSC_scr", "iPSC_TEAD24_1", "iPSC_TEAD24_2")
my_data1 = DGEList(counts=CRE_count, group=TS)
keep = filterByExpr(my_data1)
summary(keep)
my_data1=my_data1[keep,,keep.lib.sizes=FALSE]
my_data1 = calcNormFactors(my_data1, method = "TMM")
my_data1.eda = estimateGLMCommonDisp(my_data1, design = design)
fit = glmQLFit(my_data1.eda, design = design)
qlf.DEASeq = glmQLFTest(fit, contrast=c(-1,0.5,0.5))
exp.df=as.data.frame(topTags(qlf.DEASeq, nrow(qlf.DEASeq$table)))
exp.df$peakID=row.names(exp.df)
write.table(exp.df[order(exp.df$FDR),], file=gzfile("TEAD_CAGE.DE.tsv.gz"), col.names=TRUE, row.names=FALSE, sep="\t")

aCREe=read.delim(paste0(path_fig4_data,"peaks.merged.all.markalone_k16.subclass_TEAD_ONECUT.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
DE1=read.delim("TEAD_CAGE.DE.tsv.gz", header=T, stringsAsFactors = F, check.names = F)

DE1=left_join(DE1,aCREe[,c(1:4,9,11,13,15,17,19)], by=c("peakID"="V4"), copy=F)
DE1$TSSdown="No"
DE1$TSSdown[which(DE1$FDR<0.05 & DE1$logFC < (-0.5))]="Yes"
write.table(DE1,file=gzfile("TEAD_CAGE.DE.tsv.gz"), col.names=TRUE, row.names=FALSE, sep="\t")
write.table(DE1,file=gzfile(paste0(path_fig4_data,"TEAD_CAGE.DE.tsv.gz")), col.names=TRUE, row.names=FALSE, sep="\t")

#==
DE1=read.delim("TEAD_CAGE.DE.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
DE1$sig="n.s."
DE1$sig[which(DE1$FDR < 0.05 & DE1$logFC>0.5)]="Up"
DE1$sig[which(DE1$FDR < 0.05 & DE1$logFC<(-0.5))]="Down"
DE1%>%group_by(sig)%>%dplyr::summarise(count=n())
ggplot()+
  geom_point_rast(data=DE1, aes(shape=TEAD_motif, x=logFC, y=-log10(FDR), color=sig), stroke=0.1, fill="white", size=0.9, alpha=0.6)+
  scale_color_manual(values = c("Down"="blue", "n.s."="black","Up"="red"), labels=c("Down"="Down\nn=244","Up"="Up\nn=56"), name="Sig.DE")+
  scale_shape_manual(values=c("Yes"=19, "No"=22), name="contain\nmotif")+
  coord_cartesian(ylim=c(0,3.5))+
  labs(color = NULL, title="Differential expression\nfrom TEAD2/4 KD")+
  theme1+theme(legend.position = c(0.65,0.8), legend.box = "horizontal")
ggplot()+
  geom_point_rast(data=DE1[which(DE1$TEAD_KD == "Yes"),], aes(shape=TEAD_motif, x=logFC, y=-log10(FDR), color=sig), stroke=0.1, fill="white", size=0.9, alpha=0.6)+
  scale_color_manual(values = c("Down"="blue", "n.s."="black","Up"="red"),  name="Sig.DE")+
  scale_shape_manual(values=c("Yes"=19, "No"=22), name="contain\nmotif")+
  coord_cartesian(ylim=c(0,3.5))+
  labs(color = NULL, title="Differential expression\nfrom TEAD2/4 KD")+
  theme1+theme(legend.position = c(0.65,0.8), legend.box = "horizontal")

DE1%>%group_by(sig,TEAD_KD)%>%dplyr::summarise(count=n())
# sig   TEAD_KD count
# <chr> <chr>   <int>
# 1 Down  No        204
# 2 Down  Yes        40
# 3 Up    No         56
# 4 n.s.  No      21211
# 5 n.s.  Yes        88
DE1%>%group_by(sig,TEAD_motif)%>%dplyr::summarise(count=n())
# sig   TEAD_motif count
# <chr> <chr>      <int>
# 1 Down  No           143
# 2 Down  Yes          101
# 3 Up    No            40
# 4 Up    Yes           16
# 5 n.s.  No         16754
# 6 n.s.  Yes         4545
DE1%>%group_by(sig,TEAD_motif,TEAD_KD,k16ips)%>%dplyr::summarise(count=n())
#==================
# ATAC-CAGE DE xy plot
DE1_ATAC=read.delim(paste0(path_fig4_data,"TEAD24KD.DE.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
DE1_TSS=read.delim(paste0(path_fig4_data,"TEAD_CAGE.DE.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

DE1_TSS=left_join(DE1_TSS,DE1_ATAC[,c(1,2,6)], by="peakID", copy=F, suffix=c("_TSS","_ATAC"))
DE1_TSS1=DE1_TSS[which(DE1_TSS$TEAD_motif=="Yes" & DE1_TSS$k16ips %in% c("enhancer","promoter")),]

ggplot()+
  geom_point_rast(data=DE1_TSS1, aes(shape=TEAD_ChIP ,x=logFC_TSS, y=logFC_ATAC, color=k16ips), stroke=0.1, fill="white", size=0.9, alpha=0.7)+
  scale_color_manual(values = c("enhancer"="red", "promoter"="#B09C85FF"), labels=c("enhancer"="Enhancer(1492)","promoter"="Promoter(3146)"), name="Promoter-type")+
  scale_shape_manual(values=c("Yes"=19, "No"=22), name="ChIP")+
  labs(color = NULL, title="Change of CRE after TEAD2/4\nKD (with TEAD motif)", x="log2FC (Transcription)", y="log2FC (Chromatin accessibility)")+
  theme1+theme(legend.position = c(0.15,0.8), legend.box = "horizontal")

DE1_TSS1%>%group_by(k16ips)%>%dplyr::summarise(n())
cor.test(DE1_TSS1$logFC_TSS, DE1_TSS1$logFC_ATAC, method="pearson") #r 0.4359762
cor.test(DE1_TSS1$logFC_TSS[which(DE1_TSS1$TEAD_ChIP=="Yes")], DE1_TSS1$logFC_ATAC[which(DE1_TSS1$TEAD_ChIP=="Yes")], method="pearson") #r 0.5062145 


#====
aCREe=read.delim(paste0(path_fig4_data,"peaks.merged.all.markalone_k16.subclass_TEAD_ONECUT.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
DE1=read.delim("TEAD_CAGE.DE.tsv.gz", header=T, stringsAsFactors = F, check.names = F)

aCREe$TEAD_TSS_KD="No"
aCREe$TEAD_TSS_KD[which(aCREe$V4 %in% DE1$peakID[which(DE1$TSSdown == "Yes")])]="Yes"
aCREe$TEAD_TSS_KD_detected="No"
aCREe$TEAD_TSS_KD_detected[which(aCREe$V4 %in% DE1$peakID)]="Yes"
aCREe$TEAD_Motif_KD_TSS="No"
aCREe$TEAD_Motif_KD_TSS[which(aCREe$TEAD_motif == "Yes" & aCREe$TEAD_TSS_KD == "Yes")]="Yes"
write.table(aCREe,gzfile(paste0(path_fig4_data,"peaks.merged.all.markalone_k16.subclass_TEAD_ONECUT.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

aCREe1=aCREe[which(aCREe$TEAD_KD == "Yes" & aCREe$TEAD_motif == "Yes"),]
length(which(aCREe1$TEAD_TSS_KD_detected == "No")) #5066
need2=aCREe1$V4[which(aCREe1$TEAD_TSS_KD_detected == "Yes" & aCREe1$TEAD_TSS_KD == "No")] #48
length(which(aCREe1$TEAD_TSS_KD_detected == "Yes" & aCREe1$TEAD_TSS_KD == "Yes")) #27
need=aCREe$V4[which(aCREe$TEAD_KD == "No" & aCREe$TEAD_TSS_KD == "Yes" & aCREe$TEAD_motif == "Yes")] #74
                                 
aCREe_detected <- aCREe[aCREe$TEAD_motif == "Yes" & aCREe$TEAD_TSS_KD_detected == "Yes",]
tab <- table(ATAC_down = aCREe_detected$TEAD_KD, CAGE_down = aCREe_detected$TEAD_TSS_KD)
tab
fisher.test(tab)    #OR=34.17526                  
                                 
kk=DE1_TSS1[which(DE1_TSS1$peakID %in% need),]
kk2=DE1_TSS1[which(DE1_TSS1$peakID %in% need2),]

#===============================================================================
# zenbu visualize
#DE1=read.delim("TEAD_CAGE.DE.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
#DE1a=DE1[which(DE1$FDR < 0.05)]




