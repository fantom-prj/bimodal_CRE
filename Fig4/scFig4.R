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
pdf("f4b.FE_KD_ChIP.pdf", width = 1.9, height = 1.6)
print(f4b)
dev.off() 

#===============
# fig4c
maCREe=read.delim(paste0(path_fig4_data,"proportion_chromatin_state.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)

maCREe1=maCREe[which(maCREe$variable!="KD"),]%>%group_by(group, variable, chromatinState)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
maCREe1$chromatinState=gsub("_","\n",maCREe1$chromatinState)
maCREe1$chromatinState=factor(maCREe1$chromatinState, levels=c("Promoter","Bivalent\npromoter", "Active\nenhancer", "Primed\nenhancer", "Repressed\nenhancer", "Others") )
maCREe1$variable=factor(maCREe1$variable, levels=c("motif","FootPrint","ChIP","Motif_ChIP","FootPrint_ChIP","FootPrint_KD","FootPrint_ChIP_KD"))
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
pdf("f4c.proportion_chromatin_state.pdf", width = 3.2, height = 1.7)
print(f4c)
dev.off() 

#================
# fig 4d
# plot data from related CRE as subset motif activity
# Prime /active /promoter
both=read.delim(paste0(path_fig4_data,"TEAD_ONECUT.neuronal.trajectory.chromatin_state.aCRE.vs.tCRE.activity.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
both$TF=factor(both$TF, levels=c("TEAD","ONECUT"))
both$chromatin_state=factor(both$chromatin_state, levels=c("Active_enhancer","Primed_enhancer","All_promoter"))

f4d= ggplot(both[which(both$trajectory == "Trajectory 1" & both$chromatin_state %in% c("Active_enhancer","Primed_enhancer","All_promoter") & both$bimodal == "ATAC"),]) +
  stat_smooth(mapping=aes(x =pseudotime, y = value, color=group), method = "lm", formula = y ~ poly(x, 21), se = FALSE, alpha=1, linewidth=0.2) +
  scale_color_npg(labels=c("ChIP_KD"="ChIP & KD"))+
  facet_grid(row=vars(TF), col=vars(chromatin_state), scales = "free")+
  labs(color=NULL, y="Average activity", x = "Ranked by pseudotime", title = "Chromatin activity of interacting CREs") +
  theme1+theme(axis.text.x = element_blank())
pdf(paste0(path_fig4,"f4d.both.neuronal.trajectory.chromatin_state.aCRE.vs.tCRE.activity.pdf"), width = 3.3, height = 1.7)
f4d
dev.off()

#===========
# fig 4e
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

f4e1= ggplot() +
  geom_bar(data=mboth2, aes(x=variable, y=percentage, fill=value, alpha=ChIP), linewidth=0.25, color ="black", stat="identity") + 
  scale_y_continuous(labels = scales::percent)+
  facet_grid(cols=vars(TF))+
  labs(title="TF binds to closed chromatin", x=NULL, y="% of ChIP-seq peaks", fill="Opened\nchromatin")+
  scale_fill_manual(values=c("Yes"="#4DBBD5FF","No"="white"))+
  scale_alpha_manual(values=c("Yes"=1,"No"=0.5),guide=NULL)+
  geom_text(data=mboth2, aes(x=variable, label = label, y = 0.95), color = "black", size = 1.8)+
  theme1+theme(axis.text.x = element_text(angle=25, hjust=1, vjust=1), legend.position = c(0.6,0.3))
pdf(paste0(path_fig4,"f4e.ChIP_ATAC.pdf"), width = 1.7, height = 1.3)
print(f4e1)
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

f4e2= ggplot() +
  geom_line(data=overall1, aes(x=cell, y=Rel.Exp, group=1), linewidth=0.25, color ="black") + 
  facet_grid(cols=vars(motifName))+
  scale_y_continuous(breaks=c(0,0.5,1))+
  labs(title="TF expression", x=NULL)+
  theme1+theme(axis.text.x = element_text(angle=25, hjust=1, vjust=1), legend.position = c(0.6,0.6))
pdf(paste0(path_fig4,"f4e.ChIP_ATAC2.pdf"), width = 1.6, height = 0.85)
print(f4e2)
dev.off()

#=========



