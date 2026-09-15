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
chromHMM_folder=paste0(primary_folder,"Data_and_code/Fig2.ChromHMM/")
path_fig2_data=paste0(primary_folder,"Fig2/data/")

#===============================================================================
setwd(chromHMM_folder)
#=====
#generate binary file from bam files
system("sh bin.20231108.sh")
#-> output files stored in [chromHMM_folder]/data_files/all_HG38
#-> only primary chromosomes were kept

#===============================================================================
#include the CTSS as annotation
#CTSS summit inside SCAFE genuine TSS cluster from sc5nRNAseq was used

clusterbed=read.delim(paste0(primary_folder,"Data_and_code/SCAFE/sc5end/merge/annotate/sc5end.iPSC_NSC_Neuron/bed/sc5end.iPSC_NSC_Neuron.cluster.annot.bed.gz"), header=F, stringsAsFactors = F)
write.table(clusterbed[,c(1,7,8)], gzfile(paste0(chromHMM_folder,"data_files/COORDS/hg38/CTSS.hg38.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#===============================================================================
# run learnmodel from chromHMM
setwd(paste0(chromHMM_folder,"data_files/"))
system("sh learnmodel.20240305.sh")
# -> output located in [chromHMM_folder]/data_files/chromatin_state/

#===============================================================================
# plot overall figure to visualize the result of all clusters
pathl=paste0(chromHMM_folder,"data_files/chromatin_state/")

files=list.files(path=pathl, pattern ="model", recursive=T)
files_key=sapply(strsplit(files,"_"),"[",2)
files_key=as.numeric(gsub(".txt","",files_key))
beds=list.files(path=pathl, pattern="dense.bed",recursive=T)
overlaps=list.files(path=pathl, pattern="overlap.txt",recursive=T)
cell=sapply(strsplit(beds,"\\/"),"[",2)
cell=unique(sapply(strsplit(cell,"_"),"[",1))
region=list.files(path=pathl, pattern="RefSeqTSS_neighborhood.txt",recursive=T)
for (i in 1:length(files)){
  lineskip=files_key[i]*files_key[i]+files_key[i]+1
  emission=read.delim(paste0(pathl, files[i]), header=F, stringsAsFactors = F, check.names = F, skip=lineskip)
  emission1=emission[which(emission$V5 == 1),]
  emission1$label=signif(emission1$V6,2)
  emission1$label[which(emission1$V6<0.2)]=NA
  emission1$V4=factor(emission1$V4, levels=c("ATAC","CTCF","K27ac","K27m3","K4m1","K4m3"))
  emission1$V2=factor(emission1$V2)
  
  z21=ggplot(emission1, aes(y=V2, x=V4, fill=V6)) + 
    labs(y="Chromatin state cluster", x =NULL, title="Emission parameter", color=NULL)+
    geom_tile()+
    geom_text(mapping=aes(label=label), size=2.2)+
    scale_fill_gradient(low = "white", high = "#3C5488FF")+
    theme(  panel.background = element_rect(fill = "white", color = "white", size = 0.5, linetype = "solid"),
            panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = "white"), 
            axis.text.x = element_text(color ="black", angle = 90, hjust=1, vjust=0.5),
            axis.text.y = element_blank(),
            plot.title = element_text(hjust = 0.5,  color ="black", size=11),
            plot.subtitle = element_text(color="black", hjust=0.5),
            text = element_text(size=10, family="Arial"),
            strip.text = element_text(size=6, family="Arial"),
            legend.background = element_rect(fill="white", color="black"),
            legend.title = element_text(hjust=0.5),
            legend.text = element_text(lineheight = 0.8),
            axis.line = element_line(size = 0.25, colour = "black"),
            axis.ticks = element_line(size = 0.25,colour = "black"),
            strip.background =element_rect(fill="grey90"),
            legend.key = element_blank(),
            legend.position = "none",
            plot.margin = unit(c(0.1, 0.1, 1.8, 0.25), "cm"))
  
  mregion2=data.frame(matrix(nrow=0,ncol=4))
  colnames(mregion2)=c("coord","variable","value","cell")
  for (j in 1:3){
    region1=region[grep(cell[j],region)]
    region3=read.delim(paste0(pathl,region1[i]),header=F, stringsAsFactors = F, check.names = F, row.names=1)
    region3=data.frame(t(region3))
    colnames(region3)[1]="coord"
    mregion1=melt(region3, id=1)
    mregion1$cell=cell[j]
    mregion2=rbind(mregion2,mregion1)}
  mregion2$cell=factor(mregion2$cell, levels=c("iPS","NSC","NRN"))
  z4=ggplot(mregion2, aes(y=variable, x=coord, fill=value)) + 
    labs(y=NULL, x =NULL, title="RefSeq TSS", fill=NULL)+
    geom_tile()+
    #geom_text(mapping=aes(label=label), size=2.2)+
    facet_wrap(vars(cell), ncol=3, strip.position="right")+
    scale_fill_gradient(low = "white", high = "#E64B35FF")+
    theme(  panel.background = element_rect(fill = "white", color = "white", size = 0.5, linetype = "solid"),
            panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = "white"), 
            axis.text.x = element_text(color ="black", angle = 90, hjust=1, vjust=0.5),
            axis.text.y = element_blank(),
            plot.title = element_text(hjust = 0.5,  color ="black", size=11),
            plot.subtitle = element_text(color="black", hjust=0.5),
            text = element_text(size=10, family="Arial"),
            strip.text = element_text(size=6, family="Arial"),
            legend.background = element_rect(fill="white", color="black"),
            legend.title = element_text(hjust=0.5),
            legend.text = element_text(lineheight = 0.8),
            axis.line = element_line(size = 0.25, colour = "black"),
            axis.ticks = element_line(size = 0.25,colour = "black"),
            strip.background = element_blank(),
            strip.text.y = element_blank(),
            legend.key = element_blank(),
            legend.position = "none",
            plot.margin = unit(c(0.1, 0.1, 2.0, 0.25), "cm"))
  
  content2=data.frame(matrix(nrow=0,ncol=4))
  colnames(content2)=c("V4","cover","label","cell")
  for (j in 1:3){
    beds1=beds[grep(cell[j],beds)]
    content=read.delim(paste0(pathl,beds1[i]),header=F, stringsAsFactors = F, check.names = F, skip=1)
    content$cover=content$V3-content$V2
    content1=content%>%group_by(V4)%>%dplyr::summarise(cover=sum(cover)/sum(content$cover))
    content1$label=paste0(signif(content1$cover,2)*100,"%")
    content1$cell=cell[j]
    content2=rbind(content2,content1)}
  
  content2$cell=factor(content2$cell, levels=c("iPS","NSC","NRN"))
  z22=ggplot(content2, aes(y=as.factor(V4), x=cell, fill=cover)) + 
    labs(y=NULL, x =NULL, title="Genome", color=NULL)+
    geom_tile()+
    geom_text(mapping=aes(label=label), size=2.2)+
    scale_fill_gradient(low = "white", high = "#00A087FF", limits = c(0, 0.051))+
    theme(  panel.background = element_rect(fill = "white", color = "white", size = 0.5, linetype = "solid"),
            panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = "white"), 
            axis.text.x = element_text(color ="black", angle = 90, hjust=1, vjust=0.5),
            axis.text.y = element_blank(),
            plot.title = element_text(hjust = 0.5,  color ="black", size=11),
            plot.subtitle = element_text(color="black", hjust=0.5),
            text = element_text(size=10, family="Arial"),
            strip.text = element_text(size=6, family="Arial"),
            legend.background = element_rect(fill="white", color="black"),
            legend.title = element_text(hjust=0.5),
            legend.text = element_text(lineheight = 0.8),
            axis.line = element_line(size = 0.25, colour = "black"),
            axis.ticks.x = element_line(size = 0.25,colour = "black"),
            axis.ticks.y = element_line(size = 0.25,colour = "black"),
            strip.background = element_blank(),
            strip.text.y = element_blank(),
            legend.key = element_blank(),
            legend.position = "none",
            plot.margin = unit(c(0.1, 0.1, 2.1, 0.25), "cm"))
  
  feature2=data.frame(matrix(nrow=0,ncol=7))
  colnames(feature2)=c("V4","variable","value","bg_rate","value2","label","cell")
  for (j in 1:3){
    overlaps1=overlaps[grep(cell[j],overlaps)]
    feature=read.delim(paste0(pathl,overlaps1[i]), header=T, stringsAsFactors = F, check.names = F)
    colnames(feature)=sapply(strsplit(colnames(feature), "\\."),"[",1)
    feature=melt(feature[,-2], id=1)
    colnames(feature)[1]="V4"
    feature1=feature[which(feature$V4 == "Base"),]
    colnames(feature1)[3]="bg_rate"
    feature=feature[which(feature$V4 != "Base"),]
    feature=left_join(feature,feature1[,c(2,3)], by="variable", copy=F)
    feature$V4=factor(as.numeric(feature$V4), levels=c(1:files_key[i]))
    feature$value2=feature$value-feature$bg_rate
    feature$label=feature$value
    feature$label=signif(feature$label,2)
    feature$label[which(feature$value<2)]=NA
    feature$value[which(feature$value>10)]=10
    feature$cell=cell[j]
    feature2=rbind(feature2, feature)
  }
  
  feature2$cell=factor(feature2$cell, levels=c("iPS","NSC","NRN"))  
  z23=ggplot(feature2, aes(y=V4, x=variable, fill=value)) + 
    labs(y=NULL, x =NULL, title="Feature enrichment", color=NULL)+
    geom_tile()+
    #geom_text(mapping=aes(label=label), size=2.2)+
    facet_wrap(vars(cell), ncol=3, strip.position="right")+
    scale_fill_gradient(low = "white", high = "#7E6148FF")+
    theme(  panel.background = element_rect(fill = "white", color = "white", size = 0.5, linetype = "solid"),
            panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = "white"), 
            axis.text.x = element_text(color ="black", angle = 90, hjust=1, vjust=0.5),
            axis.text.y = element_blank(),
            plot.title = element_text(hjust = 0.5,  color ="black", size=11),
            plot.subtitle = element_text(color="black", hjust=0.5),
            text = element_text(size=10, family="Arial"),
            strip.text = element_text(size=6, family="Arial"),
            legend.background = element_rect(fill="white", color="black"),
            legend.title = element_text(hjust=0.5),
            legend.text = element_text(lineheight = 0.8),
            axis.line = element_line(size = 0.25, colour = "black"),
            axis.ticks.x = element_line(size = 0.25,colour = "black"),
            axis.ticks.y = element_line(size = 0.25,colour = "black"),
            strip.background = element_blank(),
            strip.text.y = element_blank(),
            legend.key = element_blank(),
            legend.position = "none",
            plot.margin = unit(c(0.1, 0.1, 0.2, 0.25), "cm"))
  jpeg(paste0(pathl,"/k", files_key[i], "/overall.jpg"), width = 9, height = 3.5, units = "in", res=600)
  grid.arrange(z21,z4,z22,z23,ncol=4, widths=c(2,2,1.2,4.5))
  dev.off() }

#===============================================================================
#use clustering k16
pathm=paste0(chromHMM_folder,"data_files/chromatin_state/k16/")

p_bed=list.files(path=pathm, pattern="dense.bed", recursive=T)
p_bedi=read.delim(paste0(pathm,p_bed[1]), header=F, stringsAsFactors = F, skip=1)
p_bedn=read.delim(paste0(pathm,p_bed[2]), header=F, stringsAsFactors = F, skip=1)
p_beds=read.delim(paste0(pathm,p_bed[3]), header=F, stringsAsFactors = F, skip=1)

p_bedi=p_bedi[which(p_bedi$V4 %in% c(2:4,7,10:16)),]
p_bedn=p_bedn[which(p_bedn$V4 %in% c(2:4,7,10:16)),]
p_beds=p_beds[which(p_beds$V4 %in% c(2:4,7,10:16)),]

labeling=data.frame(cbind(c(2:4,7,10:16),c("Active_enhancer","Active_enhancer","Primed_enhancer","Repressed_enhancer","Bivalent_promoter","Primed_enhancer","Promoter","Flanking_promoter","Promoter","CTCF","Active_enhancer")))
colnames(labeling)=c("V4","V10")
labeling$V4=as.numeric(labeling$V4)
p_bedi=left_join(p_bedi, labeling, by="V4",copy=F)
p_beds=left_join(p_beds, labeling, by="V4",copy=F)
p_bedn=left_join(p_bedn, labeling, by="V4",copy=F)

write.table(p_bedi[order(p_bedi$V1,p_bedi$V2),c(1:3,10,5:9)], gzfile(paste0(chromHMM_folder,"zenbu/markalone_iPS_k16.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(p_beds[order(p_beds$V1,p_beds$V2),c(1:3,10,5:9)], gzfile(paste0(chromHMM_folder,"zenbu/markalone_NSC_k16.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(p_bedn[order(p_bedn$V1,p_bedn$V2),c(1:3,10,5:9)], gzfile(paste0(chromHMM_folder,"zenbu/markalone_Neuron_k16.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#===============================================================================
# upload chromatin state region to zenbu
labeling1=unique(data.frame(cbind(c(3,3,5,4,2,5,1,1,1,6,3),c("Active_enhancer","Active_enhancer","Primed_enhancer","Repressed_enhancer","Bivalent_promoter","Primed_enhancer","Promoter","Flanking_promoter","Promoter","CTCF","Active_enhancer"))))
p_bedi=read.delim(paste0(chromHMM_folder,"zenbu/markalone_iPS_k16.bed.gz"), header=F, stringsAsFactors = F)
p_bedi=left_join(p_bedi, labeling1, by=c("V4"="X2"), copy=F)
write.table(p_bedi[,c(1:4,10,6:9)],gzfile(paste0(chromHMM_folder,"zenbu/markalone_iPS_k16.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)

p_beds=read.delim(paste0(chromHMM_folder,"zenbu/markalone_NSC_k16.bed.gz"), header=F, stringsAsFactors = F)
p_beds=left_join(p_beds, labeling1, by=c("V4"="X2"), copy=F)
write.table(p_beds[,c(1:4,10,6:9)],gzfile(paste0(chromHMM_folder,"zenbu/markalone_NSC_k16.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)

p_bedn=read.delim(paste0(chromHMM_folder,"zenbu/markalone_Neuron_k16.bed.gz"), header=F, stringsAsFactors = F)
p_bedn=left_join(p_bedn, labeling1, by=c("V4"="X2"), copy=F)
write.table(p_bedn[,c(1:4,10,6:9)],gzfile(paste0(chromHMM_folder,"zenbu/markalone_Neuron_k16.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)


#===============================================================================
#extract chromHMM 200bp binary files
path1=paste0(chromHMM_folder,"data_files/all_HG38/")
path2=paste0(chromHMM_folder,"data_files/binary_bed/")
files=list.files(path=path1, pattern="txt.gz")
cell.names=unique(sapply(strsplit(files, "_"),"[",1))

data_27ac=data.frame()
data_27m=data.frame()
data_m1=data.frame()
data_m3=data.frame()
data_ctcf=data.frame()

for (j in 1: length(cell.names)){
  files1=files[grep(cell.names[j], files)]
  chr.names=sapply(strsplit(files1, "_"),"[",2)
  for (i in 1: length(files1)){
    see=read.delim(paste0(path1,files1[i]), header=T, skip=1, stringsAsFactors = F, check.names = F)
    see$V1=chr.names[i]
    see$V2=seq(0,(nrow(see)-1)*200,200)
    see$V3=seq(200,(nrow(see))*200,200)
    see$V4=paste0(see$V1,"_",see$V2,"_",see$V3)
    data_27ac1=see[which(see$K27ac == 1),c(7:10)]
    data_27m1=see[which(see$K27m3 == 1),c(7:10)]
    data_m11=see[which(see$K4m1 == 1),c(7:10)]
    data_m31=see[which(see$K4m3 == 1),c(7:10)]
    data_ctcf1=see[which(see$CTCF == 1),c(7:10)]
    data_27ac=rbind(data_27ac, data_27ac1)
    data_27m=rbind(data_27m, data_27m1)
    data_m1=rbind(data_m1,data_m11)
    data_m3=rbind(data_m3,data_m31)
    data_ctcf=rbind(data_ctcf,data_ctcf1)}
  write.table(data_27ac[order(data_27ac$V1,data_27ac$V2),],gzfile(paste0(path2,cell.names[j],".27ac.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
  write.table(data_27m[order(data_27m$V1,data_27m$V2),],gzfile(paste0(path2,cell.names[j],".27m3.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
  write.table(data_m1[order(data_m1$V1,data_m1$V2),],gzfile(paste0(path2,cell.names[j],".4m1.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
  write.table(data_m3[order(data_m3$V1,data_m3$V2),],gzfile(paste0(path2,cell.names[j],".4m3.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
  write.table(data_ctcf[order(data_ctcf$V1,data_ctcf$V2),],gzfile(paste0(path2,cell.names[j],".ctcf.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)}

#===============================================================================
# prepare figure files for k16 chromatin state
state=c("Active enhancer","CTCF","Promoter","Promoter flank","Promoter","Primed enhancer","Bivalent promoter","Repressed region",
        "Genomic region","Repressed enhancer","Genomic region","Genomic region","Primed enhancer","Active enhancer","Active enhancer","Genomic region")
state_number=16:1
state_df=data.frame(state_number, state)
state_df$group="Others"
state_df$group[grep("nhancer",state_df$state)]="Enhancers"
state_df$group[grep("romoter",state_df$state)]="fPromoters"
state_df=state_df[order(state_df$group,state_df$state),]
state_df$group[grep("romoter",state_df$state)]="Promoters"

#=====
# plot overall figure #all cell types
pathl=paste0(chromHMM_folder,"data_files/chromatin_state/k16/")
beds=list.files(path=pathl, pattern="dense.bed")
overlaps=list.files(path=pathl, pattern="overlap.txt")
cell=sapply(strsplit(beds,"_"),"[",1)
region=list.files(path=pathl, pattern="RefSeqTSS_neighborhood.txt")

lineskip=16*16+16+1
emission=read.delim(paste0(pathl, "model_16.txt"), header=F, stringsAsFactors = F, check.names = F, skip=lineskip)
emission1=emission[which(emission$V5 == 1),]
emission1$label=signif(emission1$V6,2)
emission1$label[which(emission1$V6<0.2)]=NA
write.table(emission1,gzfile(paste0(path_fig2_data,"chromatin_state_emission.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)

mregion2=data.frame(matrix(nrow=0,ncol=4))
colnames(mregion2)=c("coord","variable","value","cell")
for (j in 1:3){
  region1=region[grep(cell[j],region)]
  region3=read.delim(paste0(pathl,region1),header=F, stringsAsFactors = F, check.names = F, row.names=1)
  region3=data.frame(t(region3))
  colnames(region3)[1]="coord"
  mregion1=reshape2::melt(region3, id=1)
  mregion1$cell=cell[j]
  mregion2=rbind(mregion2,mregion1)}
mregion2$variable=gsub("X","",mregion2$variable)
write.table(mregion2,gzfile(paste0(path_fig2_data,"chromatin_state_TSS.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)

content2=data.frame(matrix(nrow=0,ncol=4))
colnames(content2)=c("V4","cover","label","cell")
for (j in 1:3){
  beds1=beds[grep(cell[j],beds)]
  content=read.delim(paste0(pathl,beds1),header=F, stringsAsFactors = F, check.names = F, skip=1)
  content$cover=content$V3-content$V2
  content1=content%>%group_by(V4)%>%dplyr::summarise(cover=sum(cover)/sum(content$cover))
  content1$label=paste0(signif(content1$cover,2)*100,"%")
  content1$cell=cell[j]
  content2=rbind(content2,content1)}
write.table(content2,gzfile(paste0(path_fig2_data,"chromatin_state_genome.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)

feature2=data.frame(matrix(nrow=0,ncol=7))
colnames(feature2)=c("V4","variable","value","bg_rate","value2","label","cell")
for (j in 1:3){
  overlaps1=overlaps[grep(cell[j],overlaps)]
  feature=read.delim(paste0(pathl,overlaps1), header=T, stringsAsFactors = F, check.names = F)
  colnames(feature)=sapply(strsplit(colnames(feature), "\\."),"[",1)
  feature=reshape2::melt(feature[,-2], id=1)
  colnames(feature)[1]="V4"
  feature1=feature[which(feature$V4 == "Base"),]
  colnames(feature1)[3]="bg_rate"
  feature=feature[which(feature$V4 != "Base"),]
  feature=left_join(feature,feature1[,c(2,3)], by="variable", copy=F)
  feature$V4=factor(as.numeric(feature$V4), levels=c(1:16))
  feature$value2=feature$value-feature$bg_rate
  feature$label=feature$value
  feature$label=signif(feature$label,2)
  feature$label[which(feature$value<2)]=NA
  feature$value[which(feature$value>10)]=10
  feature$cell=cell[j]
  feature2=rbind(feature2, feature)
}
write.table(feature2,gzfile(paste0(path_fig2_data,"chromatin_state_feature.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)

#=====
setwd(paste0(chromHMM_folder,"data_files/chromatin_state/"))
system("find . -type f -name '*.bed' -exec gzip {} +")





