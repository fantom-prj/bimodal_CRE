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

primary_folder=[primary_folder]
#primary_folder="/analysisdata/fantom6/Interactome/single_cell_wallace/Figure/"
FP_output=paste0(primary_folder,"Data_and_code/Fig4.Foot_Printing/output/")
primary_FP_path="/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/subset_BAM/BINDetect_compare_output/"

ChIP_path=paste0(primary_folder,"Data_and_code/Fig4.ChIP_seq/")
path_fig4_data=paste0(primary_folder,"Fig4/data/")
path_fig4=paste0(primary_folder,"Fig4/out/")

setwd(path_fig4)

#=====================
theme1=theme(  panel.background = element_rect(fill = "white", color = "white", linewidth = 0.25, linetype = "solid"),
               panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid', color = "grey90"), 
               axis.text.x = element_text(color ="black"),
               axis.text.y = element_text(color ="black"),
               plot.title = element_text(hjust = 0.5,  color ="black", margin = margin(0.1,0.1,0.1,0.1, "cm")),
               plot.subtitle = element_text(color="black", hjust=0.5),
               text = element_text(size=6),
               strip.text = element_text(size=6),
               legend.background = element_rect(fill="white", color="black", linewidth=0.25),
               legend.title = element_text(hjust=0.5),
               legend.text = element_text(lineheight = 0.6, margin = margin(l = 1, unit = "pt")),
               legend.key.size = unit(0.2, 'cm'),
               legend.margin = margin(0.02,0.02,0.02,0.02, "cm"),
               axis.line = element_line(linewidth = 0.25, colour = "black"),
               axis.ticks = element_line(linewidth = 0.25,colour = "black"),
               strip.background =element_rect(fill="grey90"),
               strip.text.x = element_text(margin = margin(0.02,0.02,0.02,0.02, "cm")),
               strip.text.y = element_text(margin = margin(0.02,0.02,0.02,0.02, "cm")),
               legend.key = element_blank(),
               legend.position = "right",
               legend.box.margin = margin(t=-10, r=0, b=0, l=-10),
               plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
#===============================================================================
# f4a
aCREe=read.delim(paste0(path_fig4_data,"peaks.merged.all.markalone_k16.subclass_TEAD_ONECUT.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
length(which(aCREe$TEAD_motif == "Yes")) #176653
length(which(aCREe$TEAD_FootPrint == "Yes")) #93629
length(which(aCREe$TEAD_motif == "Yes" & aCREe$TEAD_ChIP == "Yes")) #4459
length(which(aCREe$TEAD_motif == "Yes" & aCREe$TEAD_KD == "Yes")) #5141
length(which(aCREe$TEAD_motif == "Yes" & aCREe$TEAD_ChIP == "Yes" & aCREe$TEAD_KD == "Yes")) #1024
length(which(aCREe$TEAD_motif == "Yes" & aCREe$TEAD_TSS_KD == "Yes")) #104
length(which(aCREe$TEAD_motif == "Yes" & aCREe$TEAD_TSS_KD == "Yes" & aCREe$TEAD_KD == "Yes")) #28

length(which(aCREe$ONECUT_motif == "Yes")) #125906
length(which(aCREe$ONECUT_FootPrint == "Yes")) #60483
length(which(aCREe$ONECUT_motif == "Yes" & aCREe$ONECUT_ChIP == "Yes")) #2050
length(which(aCREe$ONECUT_motif == "Yes" & aCREe$ONECUT_KD == "Yes")) #219
length(which(aCREe$ONECUT_motif == "Yes" & aCREe$ONECUT_ChIP == "Yes" & aCREe$ONECUT_KD == "Yes")) #28
#===============================================================================
# f4b
#FE
tead_summary=read.delim(paste0(path_fig4_data,"TEAD24KD.DE.summary.tsv.gz"), header=T, check.names = F, stringsAsFactors = F)
onecut_summary=read.delim(paste0(path_fig4_data,"ONECUT2KD.DE.summary.tsv.gz"), header=T, check.names = F, stringsAsFactors = F)
tead_summary$group="TEAD2/4 KD"
onecut_summary$group="ONECUT2 KD"
both_summary=rbind(tead_summary,onecut_summary)

both_summary$FE_sig="ns"
both_summary$FE_sig[which(both_summary$p.val<0.05)]="*"
both_summary$FE_sig[which(both_summary$p.val<0.01)]="**"
both_summary$FE_sig[which(both_summary$p.val<0.001)]="***"

both_summary$group=factor(both_summary$group, levels=c("TEAD2/4 KD","ONECUT2 KD"))
both_summary$feature1[which(both_summary$feature1 == "KD_sig_up")]="up-regulated"
both_summary$feature1[which(both_summary$feature1 == "KD_sig_down")]="down-regulated"
both_summary$feature1=factor(both_summary$feature1, levels=c("up-regulated","down-regulated"))
#both_summary=both_summary[which(both_summary$feature2 %in% c("motif","Tobias","ChIP")),]
#both_summary$feature2=gsub("Tobias","Footprint",both_summary$feature2)
both_summary$feature2=gsub("motif","Motif",both_summary$feature2)
both_summary$feature2=factor(both_summary$feature2, levels=c("ChIP","Footprint","Motif"))
both_summary$FE_sig=factor(both_summary$FE_sig, levels=c("ns","*","**","***"))

f4b=ggplot() + 
  labs(x=NULL, y = NULL, title= "Enrichment of differential CREs") +
  facet_grid(cols=vars(group))+
  geom_point(data=both_summary, mapping=aes(y=feature2, x=feature1, fill=logOdds, size=FE_sig),shape=21)+
  scale_fill_gradient2(low = "#3C5488FF", mid = "white", high = "#DC0000FF", midpoint = 0)+ 
  scale_size_manual(values=c(0.2,1.8))+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig4,"f4b.FE_KD_ChIP.pdf"), width = 1.6, height = 1.6)
print(f4b)
dev.off() 

#===============
# fig4c
maCREe=read.delim(paste0(path_fig4_data,"proportion_chromatin_state.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)

maCREe1=maCREe[-which(maCREe$variable %in% c("KD","TSS_KD_detected","TSS_KD","Motif_KD_TSS")),]%>%group_by(group, variable, chromatinState)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
maCREe1$chromatinState=gsub("_","\n",maCREe1$chromatinState)
maCREe1$chromatinState=factor(maCREe1$chromatinState, levels=c("Promoter","Bivalent\npromoter", "Active\nenhancer", "Primed\nenhancer", "Repressed\nenhancer", "Others") )
maCREe1$variable=factor(maCREe1$variable, levels=c("motif","FootPrint","ChIP","Motif_ChIP","FootPrint_ChIP","Motif_KD","FootPrint_KD","Motif_ChIP_KD","FootPrint_ChIP_KD"))
maCREe1$group=factor(maCREe1$group, levels=c("TEAD in iPSC","ONECUT in Neuron"))
maCREe2=maCREe1%>%group_by(group, variable)%>%dplyr::summarise(total=paste0("n= ", sum(count)))

f4c=ggplot() + 
  labs(x=NULL, y = "% of CRE",title= "Chromatin state of interacting CREs", fill=NULL) +
  scale_fill_npg()+
  facet_grid(cols=vars(group))+
  scale_y_continuous(labels = scales::percent, limits=c(0, 1.2), breaks=c(0, 0.25,0.5,0.75,1))+
  geom_bar(data=maCREe1, mapping=aes(y=percent, x=variable, fill=chromatinState), linewidth=0.25, color = "black", stat="identity") + 
  geom_text(data=maCREe2, mapping=aes(y=1, x=variable, label=total), size=1.8, angle=25, hjust=0, vjust=0)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig4,"f4c.proportion_chromatin_state.pdf"), width = 3.6, height = 1.7)
print(f4c)
dev.off() 

#==================
# differential ATAC with ChIP support
# f4d
DE1=read.delim(paste0(path_fig4_data,"TEAD24KD.DE.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
DE2=read.delim(paste0(path_fig4_data,"ONECUT2KD.DE.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
DE1$group="TEAD2/4 KD"
DE2$group="ONECUT2 KD"
colnames(DE1)[13]="Chromatin_state"
colnames(DE2)[15]="Chromatin_state"
both =rbind(DE1[,c("logFC","ChIP","peakID","Chromatin_state","group")],DE2[,c("logFC","ChIP","peakID","Chromatin_state","group")])

both$Chromatin_state[grep("romoter",both$Chromatin_state)]="Promoter"
both=both[which(both$Chromatin_state %in% c("Active_enhancer","Primed_enhancer","Promoter")),]
both$group=factor(both$group, levels=c("TEAD2/4 KD","ONECUT2 KD"))
both1=both%>%group_by(group,Chromatin_state)%>%dplyr::summarise(p=wilcox.test(logFC ~ChIP, alternative = "two.sided")$p.value)

f4d=ggplot() + 
  labs(x=NULL, y = "log2FC",title= "Differential accessibility with\nChIP-seq support") +
  geom_boxplot(data=both, mapping=aes(y=logFC, x=Chromatin_state, fill=ChIP), alpha=0.5, size=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.9, outlier.shape = NA) + 
  facet_grid(cols=vars(group))+
  coord_cartesian(ylim=c(-2,1))+
  scale_fill_manual(values=c("yes"="#E64B35FF", "no"="#8491B4FF"))+
  theme1+theme(axis.text.x = element_text(angle=25, hjust=1, vjust=1), legend.position = c(0.8,0.2))
pdf(paste0(path_fig4,"f4d.diff_ATAC_ChIP.pdf"), width = 1.5, height = 1.9)
print(f4d)
dev.off() 


#================
# fig 4e
# plot data from related CRE as subset ATAC activity
# Prime /active /promoter
both=read.delim(paste0(path_fig4_data,"TEAD_ONECUT.neuronal.trajectory.chromatin_state.aCRE.vs.tCRE.activity.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
both$TF=factor(both$TF, levels=c("TEAD","ONECUT"))
both=both[which(both$chromatin_state %in% c("Active_enhancer","Primed_enhancer")),]
both$chromatin_state=factor(both$chromatin_state, levels=c("Active_enhancer","Primed_enhancer"))

f4e= ggplot(both[which(both$trajectory == "Trajectory 1" & both$chromatin_state %in% c("Active_enhancer","Primed_enhancer") & both$bimodal == "ATAC"),]) +
  stat_smooth(mapping=aes(x =pseudotime, y = value, color=group), method = "lm", formula = y ~ poly(x, 21), se = FALSE, alpha=1, linewidth=0.2) +
  scale_color_manual(values=c("#E64B35FF","#3C5488FF","#91D1C2FF"), labels=c("ChIP_KD"="ChIP & KD"))+
  facet_wrap(~TF+chromatin_state, nrow=2, ncol=2, scales = "free")+
  labs(color=NULL, y="Average activity", x = "Ranked by pseudotime", title = "Chromatin activity of interacting CREs") +
  theme1+theme(axis.text.x = element_blank(), legend.position = "bottom")
pdf(paste0(path_fig4,"f4e.both.neuronal.trajectory.chromatin_state.aCRE.vs.tCRE.activity.pdf"), width = 2, height = 1.9)
f4e
dev.off()

ggplot(both[which(both$trajectory == "Trajectory 1" & both$chromatin_state %in% c("Active_enhancer","Primed_enhancer") & both$bimodal == "TSS"),]) +
  stat_smooth(mapping=aes(x =pseudotime, y = value, color=group), method = "lm", formula = y ~ poly(x, 21), se = FALSE, alpha=1, linewidth=0.2) +
  scale_color_manual(values=c("#E64B35FF","#3C5488FF","#91D1C2FF"), labels=c("ChIP_KD"="ChIP & KD"))+
  facet_wrap(~TF+chromatin_state, nrow=2, ncol=2, scales = "free")+
  labs(color=NULL, y="Average activity", x = "Ranked by pseudotime", title = "Transcription activity of interacting CREs") +
  theme1+theme(axis.text.x = element_blank())

#===========
# fig 4f
both=read.delim(paste0(path_fig4_data,"TEAD_ONECUT_ChIP_peak_merged_withmotif_overlapATAC.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
mboth=reshape2::melt(both[,c(5,2:4)],id=1)
mboth2=mboth%>%group_by(TF, variable, value)%>%dplyr::summarise(count=n())%>%mutate(percentage=count/sum(count))
mboth2$variable=gsub("overlap_ATAC","",as.character(mboth2$variable))
mboth2$variable=factor(mboth2$variable, levels=c("iPSC","NSC","Neuron"))
mboth2$TF=factor(mboth2$TF, levels=c("TEAD4","ONECUT2"))
mboth2$ChIP="No"
mboth2$ChIP[which(mboth2$TF == "TEAD4" & mboth2$variable == "iPSC")]="Yes"
mboth2$ChIP[which(mboth2$TF == "ONECUT2" & mboth2$variable == "Neuron")]="Yes"
mboth2$label=paste0(signif(mboth2$percentage*100, 2),"%")
mboth2$label[which(mboth2$value == "Yes")]=NA

f4f1= ggplot() +
  geom_bar(data=mboth2, aes(x=variable, y=percentage, fill=value, alpha=ChIP), linewidth=0.25, color ="black", stat="identity") + 
  scale_y_continuous(labels = scales::percent)+
  facet_grid(cols=vars(TF))+
  labs(title="TF binds to closed chromatin", x=NULL, y="% of ChIP-seq peaks", fill="Opened chromatin")+
  scale_fill_manual(values=c("Yes"="#4DBBD5FF","No"="white"))+
  scale_alpha_manual(values=c("Yes"=1,"No"=0.5),guide=NULL)+
  geom_text(data=mboth2, aes(x=variable, label = label, y = 0.95), color = "black", size = 1.8)+
  theme1+theme(axis.text.x = element_text(angle=25, hjust=1, vjust=1), legend.position = "bottom")
pdf(paste0(path_fig4,"f4f.ChIP_ATAC.pdf"), width = 1.7, height = 1.5)
print(f4f1)
dev.off()

#expression level of the TF
overall=read.delim(paste0(primary_folder,"Fig3/data/plot.table.all.noFP.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
overall$cell=sapply(strsplit(overall$metacell,"_"),"[",1)
overall1=overall%>%group_by(motifName, cell)%>%dplyr::summarise(expression=mean(expression_TF1))%>%dplyr::mutate(Rel.Exp=expression/max(expression))
overall1=overall1[which(overall1$motifName %in% c("TEAD4","ONECUT2")),]

overall1$motifName=factor(overall1$motifName, levels=c("TEAD4","ONECUT2"))
overall1$cell[which(overall1$cell == "iPS")]="iPSC"
overall1$cell[which(overall1$cell == "neuron")]="Neuron"
overall1$cell=factor(overall1$cell, levels=c("iPSC","NSC","Neuron"))

f4f2= ggplot() +
  geom_line(data=overall1, aes(x=cell, y=Rel.Exp, group=1), linewidth=0.25, color ="black") + 
  facet_grid(cols=vars(motifName))+
  scale_y_continuous(breaks=c(0,0.5,1))+
  labs(title="TF expression", x=NULL)+
  theme1+theme(axis.text.x = element_text(angle=25, hjust=1, vjust=1), legend.position = c(0.6,0.6))
pdf(paste0(path_fig4,"f4f.ChIP_ATAC2.pdf"), width = 1.6, height = 0.85)
print(f4f2)
dev.off()

#==================
# f4g
# ATAC-CAGE DE xy plot
DE1_ATAC=read.delim(paste0(path_fig4_data,"TEAD24KD.DE.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
DE1_TSS=read.delim(paste0(path_fig4_data,"TEAD_CAGE.DE.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

DE1_TSS=inner_join(DE1_TSS,DE1_ATAC[,c(1,2,6)], by="peakID", copy=F, suffix=c("_TSS","_ATAC"))
DE1_TSS1=DE1_TSS[which(DE1_TSS$TEAD_motif=="Yes" & DE1_TSS$k16ips %in% c("enhancer","promoter")),]
DE1_TSS1%>%group_by(k16ips)%>%dplyr::summarise(count=n())
DE1_TSS1[which(DE1_TSS1$TEAD_ChIP == "Yes"),]%>%group_by(k16ips)%>%dplyr::summarise(count=n())
cor.test(DE1_TSS1$logFC_TSS, DE1_TSS1$logFC_ATAC)# 0.4275399
cor.test(DE1_TSS1$logFC_TSS[which(DE1_TSS1$k16ips=="enhancer")], DE1_TSS1$logFC_ATAC[which(DE1_TSS1$k16ips=="enhancer")]) #0.4625973
cor.test(DE1_TSS1$logFC_TSS[which(DE1_TSS1$k16ips=="promoter")], DE1_TSS1$logFC_ATAC[which(DE1_TSS1$k16ips=="promoter")]) #0.2285357
DE1_TSS1a=DE1_TSS1[which(DE1_TSS1$TEAD_ChIP == "Yes"),]
cor.test(DE1_TSS1a$logFC_TSS[which(DE1_TSS1a$k16ips=="enhancer")], DE1_TSS1a$logFC_ATAC[which(DE1_TSS1a$k16ips=="enhancer")]) #0.5693516
cor.test(DE1_TSS1a$logFC_TSS[which(DE1_TSS1a$k16ips=="promoter")], DE1_TSS1a$logFC_ATAC[which(DE1_TSS1a$k16ips=="promoter")]) #0.2741366 

f4g=ggplot()+
  geom_point_rast(data=DE1_TSS1a, aes(x=logFC_TSS, y=logFC_ATAC, color=k16ips), shape=19, stroke=0.1, fill="white", size=0.8, alpha=0.6)+
  facet_grid(rows=vars(k16ips))+
  scale_color_manual(values = c("enhancer"="red", "promoter"="#B09C85FF"), labels=c("enhancer"="En(301)\nr = 0.57","promoter"="P(841)\nr = 0.27"), name="CRE-type")+
  #scale_shape_manual(values=c("Yes"=19, "No"=22), name="ChIP", labels=c("Yes"="Yes (r=0.50)", "No"="No (r=0.41)"))+
  labs(color = NULL, title="Change of CRE after TEAD2/4\nKD (with ChIP support)", x="log2FC (Transcription)", y="log2FC (Chromatin accessibility)")+
  theme1+theme(legend.position = "right")
pdf(paste0(path_fig4,"f4g.KD_ATAC_CAGE_TEAD_compare.pdf"), width = 1.6, height = 1.7)
print(f4g)
dev.off()

#==================
# f4g
aCREe=read.delim(paste0(path_fig4_data,"peaks.merged.all.markalone_k16.subclass_TEAD_ONECUT.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
aCREe1=aCREe[which(aCREe$TEAD_motif == "Yes"),]
aCREe1=aCREe1[which(aCREe1$TEAD_KD=="Yes" | aCREe1$TEAD_TSS_KD == "Yes"),]
aCREe1$group=""
aCREe1$group[which(aCREe1$TEAD_TSS_KD=="Yes")]="KD_TSS"
aCREe1$group[which(aCREe1$TEAD_KD=="Yes")]="KD_ATAC"
aCREe1$group[which(aCREe1$TEAD_KD=="Yes" & aCREe1$TEAD_TSS_KD == "Yes")]="KD_both"
aCREe1$group=factor(aCREe1$group, levels=c("KD_ATAC","KD_both","KD_TSS"))
aCREe1$TEAD_TSS_KD_detected[which(aCREe1$TEAD_TSS_KD_detected=="Yes")]="KD_TSS_evaluable"
aCREe1$TEAD_TSS_KD_detected[which(aCREe1$TEAD_TSS_KD_detected=="No")]="KD_TSS_non_evaluable"

aCREe3=aCREe1%>%group_by(TEAD_TSS_KD_detected,group)%>%dplyr::summarise(count=n())%>%mutate(percent=count/sum(count))
kk=aCREe1%>%group_by(TEAD_TSS_KD_detected,group,k16ips)%>%dplyr::summarise(count=n())%>%mutate(percent=count/sum(count))
kk1=aCREe1%>%group_by(TEAD_TSS_KD_detected,k16ips)%>%dplyr::summarise(count=n())%>%mutate(percent=count/sum(count))

f4h1=ggplot() + 
  labs(x=NULL, y = NULL, title= "ATAC and TSS profile from TEAD2/4 KD", fill=NULL) +
  scale_fill_manual(values=c("KD_ATAC"="#91D1C2FF", "KD_both"="purple", "KD_TSS"="firebrick3"))+
  scale_x_continuous(labels = scales::percent, breaks=c(0, 0.25,0.5,0.75,1))+
  geom_bar(data=aCREe3, mapping=aes(x=percent, y=TEAD_TSS_KD_detected, fill=group), position = position_stack(reverse = TRUE), linewidth=0.25, width=0.7, alpha=0.5, color = "black", stat="identity") + 
  geom_text(data=aCREe3, mapping=aes(y=TEAD_TSS_KD_detected, x=percent, label=paste0("n = ",count)), position = position_fill(vjust = 0.5), size=1.8)+
  theme1+theme(legend.position = "bottom")
pdf(paste0(path_fig4,"f4h1.TEAD24_KD_ATAC_TSS_summary.pdf"), width = 2.4, height = 0.7)
print(f4h1)
dev.off() 

DE1_ATAC=read.delim(paste0(path_fig4_data,"TEAD24KD.DE.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
DE1_TSS=read.delim(paste0(path_fig4_data,"TEAD_CAGE.DE.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
aCREe1=left_join(aCREe1,DE1_ATAC[,c(6,1)], by=c("V4"="peakID"),copy=F)
aCREe1=left_join(aCREe1,DE1_TSS[,c(6,1)], by=c("V4"="peakID"),copy=F, suffix=c("_ATAC","_TSS"))                         

aCREe2=reshape2::melt(aCREe1[,c(1:4,32,34:36)], id=c(1:6))
aCREe2$TEAD_TSS_KD_detected=gsub("KD_TSS_","",aCREe2$TEAD_TSS_KD_detected)
aCREe2$group2=aCREe2$group
aCREe2$group2[which(aCREe2$TEAD_TSS_KD_detected=="non_evaluable")]=paste0(aCREe2$TEAD_TSS_KD_detected[which(aCREe2$TEAD_TSS_KD_detected=="non_evaluable")],"\n",aCREe2$group[which(aCREe2$TEAD_TSS_KD_detected=="non_evaluable")])
aCREe2$group2=factor(aCREe2$group2, levels =c("non_evaluable\nKD_ATAC","KD_ATAC","KD_both","KD_TSS"))
aCREe2$k16ips[which(aCREe2$k16ips %in% c("CTCF",NA))]="Others"
aCREe2$variable=gsub("logFC_","",aCREe2$variable)

aCREe4=unique(aCREe2[,c(1:6,9)])%>%group_by(group2,k16ips)%>%dplyr::summarise(count=n())%>%mutate(percent=count/sum(count))
aCREe4$variable="CRE-type"
aCREe4$k16ips=factor(aCREe4$k16ips, levels=c("Others","promoter","enhancer"))
panel_order = c("ATAC", "TSS", "CRE-type")
aCREe2$variable = factor(aCREe2$variable, levels = panel_order)
aCREe4$variable = factor(aCREe4$variable, levels = panel_order)

library(ggh4x)
f4h2=ggplot() + 
  labs(x=NULL, y=NULL, title= NULL) +
  geom_boxplot(data=aCREe2, mapping=aes(y=group2, x=value, fill=group), alpha=0.5, linewidth=0.25, color = "black", outlier.colour="grey25", notch = FALSE, width = 0.7, outlier.shape = NA) + 
  geom_col(data=aCREe4, mapping=aes(y=group2, x =percent, fill=k16ips), alpha=0.6, linewidth=0.25, color = "black", width = 0.7)+
  scale_y_discrete(limits = rev)+
  facet_grid(cols=vars(variable), scales = "free_x")+
  ggh4x::facetted_pos_scales(x = list(variable == "CRE-type" ~ scale_x_continuous(breaks = seq(0, 1, 0.25), labels = scales::label_percent())))+
  scale_fill_manual(values=c("KD_ATAC"="#91D1C2FF", "KD_both"="purple", "KD_TSS"="firebrick3", "promoter"="#B09C85FF", "enhancer"="red", "Others"="white"), labels=c("enhancer"="Enhancer","promoter"="Promoter"), name="CRE-type")+
  theme1+theme(axis.text.x = element_text(angle=25, hjust=1, vjust=1))

pdf(paste0(path_fig4,"f4h2.TEAD24_KD_ATAC_TSS_summary.pdf"), width = 3, height = 1)
print(f4h2)
dev.off() 


#=================
# fig f4i
#plot the result of GO/Rea/KEGG
indirect=read.delim("indirect_all_summary.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
direct=read.delim("direct_all_summary.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
indirect$group="indirect"
direct$group="direct"
indirect=indirect[,-c(14)] #remove ABC_cutoff column
both=rbind(indirect, direct)
both=both[-grep("KD",both$test),]
both$mlogq=-log10(both$qvalue)
sboth=spread(both[,c(2,15,20,21)], key=3, value=4)
sboth$direct[is.na(sboth$direct)]=0
mboth=reshape2::melt(sboth, id=c(1,2))
mboth$test=factor(mboth$test, levels=c("TEAD4_ChIP","ONECUT2_ChIP"))

f4i=ggplot() + 
  labs(x="-log10(q-value)", y = "Pathways",title= "Enriched pathways mediated by TEAD4 & ONECUT2", fill=NULL) +
  scale_fill_npg()+
  facet_wrap(vars(test),nrow=2, scales="free_y", space="free_y")+
  geom_bar(data=mboth, mapping=aes(y=reorder(Description,value), x=value, fill=variable), linewidth=0.25, position_dodge(), color = "black", stat="identity") + 
  theme1+theme(legend.position = c(0.8,0.5))
pdf(paste0(path_fig4,"f4i.GO_enrichment.pdf"), width = 3.2, height = 1.7)
print(f4i)
dev.off() 

#=================
# ATAC log2CPM before and after KD for significant CRE, divided into chromatin state
path5=paste0(KDATAC_folder,"edgeR/")
#ATACcount=read.delim("/analysisdata/fantom6/Interactome/single_cell_wallace/Figure/Data_and_code/Fig4.KD_ATAC/snATAC_peak_counts.tsv.gz", header=T, stringsAsFactors = F, check.names = F, rownames=1)
TEAD_KD_ATAC=read.delim(paste0(KDATAC_folder,"edgeR/TEAD24KD.norm.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
TEAD_KD_ATAC$siCtrl_CPM=(TEAD_KD_ATAC$iPSC_scr_rep1+TEAD_KD_ATAC$iPSC_scr_rep2)/2
TEAD_KD_ATAC$siTEAD_CPM=(TEAD_KD_ATAC$iPSC_TEAD24_1_rep1+TEAD_KD_ATAC$iPSC_TEAD24_1_ep2+TEAD_KD_ATAC$iPSC_TEAD24_2_rep1+TEAD_KD_ATAC$iPSC_TEAD24_2_rep2)/4

aCREe=read.delim(paste0(path_fig4_data,"peaks.merged.all.markalone_k16.subclass_TEAD_ONECUT.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
aCREe1=aCREe[which(aCREe$TEAD_motif == "Yes"),]
aCREe1=aCREe1[which(aCREe1$TEAD_KD=="Yes"),]
aCREe1=left_join(aCREe1,TEAD_KD_ATAC[,c(7:9)], by=c("V4"="peakID"),copy=F)
aCREe1a=aCREe1%>%group_by(k16ips_all)%>%dplyr::summarise(ctrl=median(siCtrl_CPM),tead=median(siTEAD_CPM))

#==================






