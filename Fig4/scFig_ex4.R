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
FP_output=paste0(primary_folder,"Data_and_code/Foot_Printing/output/")
primary_FP_path="/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/subset_BAM/BINDetect_compare_output/"

ChIP_path=paste0(primary_folder,"Data_and_code/ChIP_seq/")

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
# ex4a
# TEAD4 heatmap
TEAD4final0=read.delim(paste0(path_fig4_data,"heatmap.aCRE.TEAD4.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
colnames(TEAD4final0)[2:4]=c("iPSC","NSC","Neuron")
rownames(TEAD4final0)=TEAD4final0$peakID
TEAD4final0=TEAD4final0[,-1]
TEAD4final0 <- TEAD4final0 %>% mutate(across(c(iPSC, NSC, Neuron),
    ~ case_when(.x == "Promoter" ~ 3, .x == "Bivalent_promoter" ~ 3, .x == "Active_enhancer" ~ 2, 
                .x == "Repressed_enhancer" ~ 1, .x == "Primed_enhancer" ~ 1, TRUE ~ 0)))

TEAD4final0[, 1:3] = TEAD4final0[, 1:3] * 10
row_dist0 <- dist(TEAD4final0)
row_clust0 <- hclust(row_dist0)
row_order0 <- row_clust0$labels[row_clust0$order]
dend_row0 <- as.dendrogram(row_clust0)
p_rowtree0 <- ggdendrogram(dend_row0, rotate = T)+ 
  scale_y_reverse()+ 
  theme_void()+
  theme(plot.margin = unit(c(0.18, 0, 0, 0), "cm"))
p_rowtree0$layers[[2]]$aes_params$size <- 0.1

TEAD4final0$peakID=rownames(TEAD4final0)
TEAD4final0m=reshape2::melt(TEAD4final0, id=760)
TEAD4final0m$group=sapply(strsplit(as.character(TEAD4final0m$variable),"_kk"),"[",2)
TEAD4final0m$SEACell=sapply(strsplit(as.character(TEAD4final0m$variable),"_kk"),"[",1)
TEAD4final0m=left_join(TEAD4final0m, pseudotime1[,c(1,5)],by="SEACell",copy=F)
TEAD4final0m$peakID=factor(TEAD4final0m$peakID, levels=row_order0)

df0 = TEAD4final0m[is.na(TEAD4final0m$group), ]
df0$value = factor(df0$value, levels = c("10","20","30","0"))
df0$group="Chromatin state"

out00=ggplot(df0, aes(x = variable, y = peakID, fill = value)) +
  labs(x=NULL, y=NULL, fill="Chromatin state")+
  scale_fill_manual(labels=c("10"="Primed / repressed\nenhancer","20"="Active enhancer","30"="Promoter","0"="NA"), values=c("#00A087FF","#E64B35FF","#B09C85FF","grey"))+
  facet_grid(cols=vars(group))+
  geom_tile_rast() + 
  new_scale_fill() +
  geom_tile_rast(data = unique(df0[,c(2,4)]), aes(x = variable, y = length(levels(TEAD4final0m$peakID)) + 100, fill = variable), height = 100,inherit.aes = FALSE)+
  scale_fill_manual(values = c("iPSC"="lightgrey" , "NSC"="darkgrey","Neuron"="black"), name="Cell-type") +
  theme1 +
  theme(panel.grid.major = element_blank(),  axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), legend.position="right", strip.background =element_blank())

df0b=TEAD4final0m[which(TEAD4final0m$group == "ATAC"), ]
df0b$SEACell=factor(df0b$SEACell, levels=pseudotime1$SEACell)
anno_df0b=unique(df0b[,c("SEACell","label")])
anno_df0b$label=factor(anno_df0b$label, levels=cluster_info$label)
out00a=ggplot(df0b, aes(x = SEACell, y = peakID, fill = value)) +
  labs(x=NULL, y=NULL, fill="ATAC")+
  scale_fill_gradientn(colours = c("white", "firebrick3", "firebrick3"),values = scales::rescale(c(0, 1.4, 1.72)), limits = c(0, 1.72),breaks = c(0, 0.5, 1, 1.5))+
  #scale_fill_gradient2(low = "white",high = "firebrick3", limits = c(0,0.6), oob = scales::squish) +
  facet_grid(cols=vars(group))+
  geom_tile_rast() + 
  new_scale_fill() +
  geom_tile_rast(data = anno_df0b, aes(x = SEACell, y = length(levels(TEAD4final0m$peakID)) + 100, fill = label), height = 100,inherit.aes = FALSE)+
  scale_fill_manual(values = c("iPSC_Sphase"="#E64B35FF","iPSC_G2M"="#4DBBD5FF", "NSC_stem"="#00A087FF", "NSC_1"="#3C5488FF", "NSC_2"="#F39B7FFF","NPC_like"="#8491B4FF", "Immature_neuron"="#91D1C2FF", "Mature_neuron_1"="#DC0000FF", "Mature_neuron_2"="#7E6148FF"), name="Cell-type", guide=NULL) +
  theme1 +
  theme(panel.grid.major = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), legend.position="right", strip.background =element_blank())

df0c=TEAD4final0m[which(TEAD4final0m$group == "TSS"), ]
df0c$SEACell=factor(df0c$SEACell, levels=pseudotime1$SEACell)
anno_df0c=unique(df0c[,c("SEACell","label")])
anno_df0c$label=factor(anno_df0c$label, levels=cluster_info$label)
out00b=ggplot(df0c, aes(x = SEACell, y = peakID, fill = value)) +
  labs(x=NULL, y=NULL, fill="TSS")+
  scale_fill_gradientn(colours = c("white", "navy", "navy"),values = scales::rescale(c(0, 2, 4.6)), limits = c(0, 4.6),breaks = c(0,  1, 2, 3,4))+
  #scale_fill_gradient2(low = "white",high = "navy", limits = c(0,1)) +
  facet_grid(cols=vars(group))+
  geom_tile_rast() + 
  new_scale_fill() +
  geom_tile_rast(data = anno_df0c, aes(x = SEACell, y = length(levels(TEAD4final0m$peakID)) + 100, fill = label), height = 100,inherit.aes = FALSE)+
  scale_fill_manual(values = c("iPSC_Sphase"="#E64B35FF","iPSC_G2M"="#4DBBD5FF", "NSC_stem"="#00A087FF", "NSC_1"="#3C5488FF", "NSC_2"="#F39B7FFF","NPC_like"="#8491B4FF", "Immature_neuron"="#91D1C2FF", "Mature_neuron_1"="#DC0000FF", "Mature_neuron_2"="#7E6148FF"), name="Cell-type") +
  theme1 +
  theme(panel.grid.major = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), legend.position="right", strip.background =element_blank())

pdf(paste0(path_fig4,"ex4a.heatmap.aCRE.TEAD4.pdf"), width=4.6, height=2)
print(plot_grid(p_rowtree0, out00, out00a, out00b, ncol = 4, rel_widths = c(1, 9.5,8,10.6)))
dev.off()

#==================
# ex4a
# ONECUT both with and without transcription
ONECUTfinal0=read.delim(paste0(path_fig4_data,"heatmap.aCRE.ONECUT.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
colnames(ONECUTfinal0)[2:4]=c("iPSC","NSC","Neuron")
rownames(ONECUTfinal0)=ONECUTfinal0$peakID
ONECUTfinal0=ONECUTfinal0[,-1]
ONECUTfinal0 <- ONECUTfinal0 %>% mutate(across(c(iPSC, NSC, Neuron),
                                             ~ case_when(.x == "Promoter" ~ 3, .x == "Bivalent_promoter" ~ 3, .x == "Active_enhancer" ~ 2, 
                                                         .x == "Repressed_enhancer" ~ 1, .x == "Primed_enhancer" ~ 1, TRUE ~ 0)))

ONECUTfinal0[, 1:3] = ONECUTfinal0[, 1:3] * 10
row_dist4 <- dist(ONECUTfinal0)
row_clust4 <- hclust(row_dist4)
row_order4 <- row_clust4$labels[row_clust4$order]
dend_row4 <- as.dendrogram(row_clust4)
p_rowtree4 <- ggdendrogram(dend_row4, rotate = T)+ 
  scale_y_reverse()+ 
  theme_void()+
  theme(plot.margin = unit(c(0.2, 0, 0, 0), "cm"))
p_rowtree4$layers[[2]]$aes_params$size <- 0.1

ONECUTfinal0$peakID=rownames(ONECUTfinal0)
ONECUTfinal0m=reshape2::melt(ONECUTfinal0, id=760)
ONECUTfinal0m$group=sapply(strsplit(as.character(ONECUTfinal0m$variable),"_kk"),"[",2)
ONECUTfinal0m$SEACell=sapply(strsplit(as.character(ONECUTfinal0m$variable),"_kk"),"[",1)
ONECUTfinal0m=left_join(ONECUTfinal0m, pseudotime1[,c(1,5)],by="SEACell",copy=F)
ONECUTfinal0m$peakID=factor(ONECUTfinal0m$peakID, levels=row_order4)

df4 = ONECUTfinal0m[is.na(ONECUTfinal0m$group), ]
df4$value = factor(df4$value, levels = c("10","20","30","0"))
df4$group="Chromatin state"

out04=ggplot(df4, aes(x = variable, y = peakID, fill = value)) +
  labs(x=NULL, y=NULL, fill="Chromatin state")+
  scale_fill_manual(labels=c("10"="Primed / repressed\nenhancer","20"="Active enhancer","30"="Promoter","0"="NA"), values=c("#00A087FF","#E64B35FF","#B09C85FF","grey"))+
  facet_grid(cols=vars(group))+
  geom_tile_rast() + 
  new_scale_fill() +
  geom_tile_rast(data = unique(df4[,c(2,4)]), aes(x = variable, y = length(levels(ONECUTfinal0m$peakID)) + 50, fill = variable), height = 50,inherit.aes = FALSE)+
  scale_fill_manual(values = c("iPSC"="lightgrey" , "NSC"="darkgrey","Neuron"="black"), name="Cell-type") +
  theme1 +
  theme(panel.grid.major = element_blank(),  axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), legend.position="right", strip.background =element_blank())

df4b=ONECUTfinal0m[which(ONECUTfinal0m$group == "ATAC"), ]
df4b$SEACell=factor(df4b$SEACell, levels=pseudotime1$SEACell)
anno_df4b=unique(df4b[,c("SEACell","label")])
anno_df4b$label=factor(anno_df4b$label, levels=cluster_info$label)
out04a=ggplot(df4b, aes(x = SEACell, y = peakID, fill = value)) +
  labs(x=NULL, y=NULL, fill="ATAC")+
  scale_fill_gradientn(colours = c("white", "firebrick3", "firebrick3"),values = scales::rescale(c(0, 1.3, 1.42)), limits = c(0, 1.5),breaks = c(0, 0.5, 1, 1.5))+
  facet_grid(cols=vars(group))+
  geom_tile_rast() + 
  new_scale_fill() +
  geom_tile_rast(data = anno_df4b, aes(x = SEACell, y = length(levels(ONECUTfinal0m$peakID)) + 50, fill = label), height = 50,inherit.aes = FALSE)+
  scale_fill_manual(values = c("iPSC_Sphase"="#E64B35FF","iPSC_G2M"="#4DBBD5FF", "NSC_stem"="#00A087FF", "NSC_1"="#3C5488FF", "NSC_2"="#F39B7FFF","NPC_like"="#8491B4FF", "Immature_neuron"="#91D1C2FF", "Mature_neuron_1"="#DC0000FF", "Mature_neuron_2"="#7E6148FF"), name="Cell-type", guide=NULL) +
  theme1 +
  theme(panel.grid.major = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), legend.position="right", strip.background =element_blank())

df4c=ONECUTfinal0m[which(ONECUTfinal0m$group == "TSS"), ]
df4c$SEACell=factor(df4c$SEACell, levels=pseudotime1$SEACell)
anno_df4c=unique(df4c[,c("SEACell","label")])
anno_df4c$label=factor(anno_df4c$label, levels=cluster_info$label)
out04b=ggplot(df4c, aes(x = SEACell, y = peakID, fill = value)) +
  labs(x=NULL, y=NULL, fill="TSS")+
  scale_fill_gradientn(colours = c("white", "navy", "navy"),values = scales::rescale(c(0, 2, 4.7)), limits = c(0, 4.7),breaks = c(0,  1, 2, 3,4))+
  #scale_fill_gradient2(low = "white",mid = "navy",high = "navy", midpoint = 2.5) +
  facet_grid(cols=vars(group))+
  geom_tile_rast() + 
  new_scale_fill() +
  geom_tile_rast(data = anno_df4c, aes(x = SEACell, y = length(levels(ONECUTfinal0m$peakID)) + 50, fill = label), height = 50,inherit.aes = FALSE)+
  scale_fill_manual(values = c("iPSC_Sphase"="#E64B35FF","iPSC_G2M"="#4DBBD5FF", "NSC_stem"="#00A087FF", "NSC_1"="#3C5488FF", "NSC_2"="#F39B7FFF","NPC_like"="#8491B4FF", "Immature_neuron"="#91D1C2FF", "Mature_neuron_1"="#DC0000FF", "Mature_neuron_2"="#7E6148FF"), name="Cell-type") +
  theme1 +
  theme(panel.grid.major = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), legend.position="right", strip.background =element_blank())

pdf(paste0(path_fig4,"ex4a.heatmap.aCRE.ONECUT2.pdf"), width=4.6, height=2)
print(plot_grid(p_rowtree4, out04, out04a, out04b, ncol = 4, rel_widths = c(1, 9.5,8,10.6)))
dev.off()

#===============================================================================
# ex4b
# directly extracted from footprint results

#===============================================================================
# ex4c
library(pROC)

tead_FP1=read.delim(paste0(path_fig4_data,"tobias_chip_tead4_cluster1.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
roc_data1 <- roc(tead_FP1$ChIP, tead_FP1$V13, levels=c("Yes","No"))
auc1=auc(roc_data1) #0.937
roc_df1 <- data.frame(
  specificity = roc_data1$specificities,
  sensitivity = roc_data1$sensitivities,
  threshold = roc_data1$thresholds)
TEADt=ggplot(roc_df1, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(color = "blue", size=0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth=0.2)+
  annotate("text", x = 0.4, y = 0.75, label = paste("AUC =", signif(auc1, 3)), size = 2, color = "black")+
  labs(title = "ROC Curve against TEAD4 ChIP-seq", x = "False Positive Rate (1 - Specificity)", y = "True Positive Rate (Sensitivity)") +
  theme1
pdf(paste0(path_fig4,"ex4c.Tobias_ChIP_ROC.cluster1.pdf"), width = 1.7, height = 1.7)
print(TEADt)
dev.off() 

#=======
#ONECUT2
onecut_FP1=read.delim(paste0(path_fig4_data,"tobias_chip_onecut2_cluster7.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
roc_data3 <- roc(onecut_FP1$ChIP, onecut_FP1$V13, levels=c("Yes","No"))
auc3=auc(roc_data3) #0.9302

roc_df3 <- data.frame(
  specificity = roc_data3$specificities,
  sensitivity = roc_data3$sensitivities,
  threshold = roc_data3$thresholds)

OCt=ggplot(roc_df3, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(color = "purple", size=0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size=0.2)+
  annotate("text", x = 0.4, y = 0.75, label = paste("AUC =", signif(auc3, 3)), size = 2, color = "black")+
  labs(title = "ROC Curve against ONECUT2 ChIP-seq", x = "False Positive Rate (1 - Specificity)", y = "True Positive Rate (Sensitivity)") +
  theme1
pdf(paste0(path_fig4,"ex4c.Tobias_ChIP_ROC.cluster7.pdf"), width = 1.7, height = 1.7)
print(OCt)
dev.off() 


#===============================================================================
#ex4d
qpcr=read.delim(paste0(path_fig4_data, "qpcr.txt"), header=T, stringsAsFactors = F, check.names = F)
qpcr$Cell[which(qpcr$Cell == "Cortical neuron")]="Neuron"
qpcr$Primer=gsub("\\.p","\nprimer",qpcr$Primer)

qpcr1=qpcr%>%group_by(Cell, Sample, Primer)%>%dplyr::summarise(mean_FC=mean(FC), SD_FC=sd(FC))
qpcr1$group="Sample"
qpcr1$group[grep("NC1",qpcr1$Sample)]="Control"
qpcr1$group=factor(qpcr1$group, levels=c("Control","Sample"))
qpcr1$Cell=factor(qpcr1$Cell, levels=c("iPSC","Neuron"))
qpcr1$Sample[which(qpcr1$group == "Control")]=gsub("ONECUT2 ","",qpcr1$Sample[which(qpcr1$group == "Control")])
qpcr1$Sample[which(qpcr1$group == "Control")]=gsub("TEAD2/4 ","",qpcr1$Sample[which(qpcr1$group == "Control")])

ex4d= ggplot(qpcr1, aes(x=Sample, y=mean_FC, fill=group)) +
  geom_bar(linewidth=0.25, color =NA, stat="identity", alpha=0.7) + 
  geom_errorbar(aes(ymin=mean_FC-SD_FC, ymax=mean_FC+SD_FC), linewidth=0.15, width=0.25)+
  facet_wrap(Cell ~ Primer, ncol=6, scale="free_x")+
  labs(title="Knockdown efficiency by RT-qPCR", x=NULL, y="Relative fold change compared to NC")+
  scale_fill_manual(values=c("grey","black"),guide=NULL)+
  theme1+theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1))
pdf(paste0(path_fig4,"ex4d.qpcr.pdf"), width = 4.1, height = 2)
print(ex4d)
dev.off()

#===============================================================================
#ex4e
#volcano plot again motif or not
DE1=read.delim(paste0(path_fig4_data,"TEAD24KD.DE.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
DE1$sig="n.s."
DE1$sig[which(DE1$FDR < 0.05 & DE1$logFC>0.5)]="Up"
DE1$sig[which(DE1$FDR < 0.05 & DE1$logFC<(-0.5))]="Down"
DE2=read.delim(paste0(path_fig4_data,"ONECUT2KD.DE.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
DE2$sig="n.s."
DE2$sig[which(DE2$FDR < 0.05 & DE2$logFC>0.5)]="Up"
DE2$sig[which(DE2$FDR < 0.05 & DE2$logFC<(-0.5))]="Down"
DE1%>%group_by(sig)%>%dplyr::summarise(count=n())
DE2%>%group_by(sig)%>%dplyr::summarise(count=n())

v1=ggplot()+
  geom_point_rast(data=DE1, aes(shape=Motif, x=logFC, y=-log10(FDR), color=sig), stroke=0.1, fill="white", size=0.5, alpha=0.4)+
  scale_color_manual(values = c("Down"="blue", "n.s."="black","Up"="red"), labels=c("Down"="Down\nn=7300","Up"="Up\nn=613"), name="Sig.DE")+
  scale_shape_manual(values=c("yes"=19, "no"=22), name="contain\nmotif")+
  coord_cartesian(ylim=c(0,25))+
  labs(color = NULL, title="Differential accessibility\nfrom TEAD2/4 KD")+
  theme1+theme(legend.position = c(0.65,0.8), legend.box = "horizontal")
v2=ggplot()+
  geom_point_rast(data=DE2, aes(shape=Motif, x=logFC, y=-log10(FDR), color=sig), stroke=0.1, fill="white", size=0.5, alpha=0.4)+
  scale_color_manual(values = c("Down"="blue", "n.s."="black",'Up'="red"), labels=c("Down"="Down\nn=253","Up"="Up\nn=37"), name="Sig.DE")+
  scale_shape_manual(values=c("yes"=19, "no"=22), name="contain\nmotif")+
  coord_cartesian(ylim=c(0,5))+
  labs(color = NULL, title="Differential accessibility\nfrom ONECUT2 KD")+
  theme1+theme(legend.position = c(0.65,0.8), legend.box = "horizontal")
pdf(paste0(path_fig4,"ex4e.volcano_KD_both.pdf"), width = 3.2, height = 2)
grid.arrange(arrangeGrob(v1,v2, ncol=2, nrow = 1, widths = c(2,2)))
dev.off()

#==================
# fig ex4f

DE1=read.delim(paste0(path_fig4_data,"TEAD_CAGE.DE.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
DE1$sig="n.s."
DE1$sig[which(DE1$FDR < 0.05 & DE1$logFC>0.5)]="Up"
DE1$sig[which(DE1$FDR < 0.05 & DE1$logFC<(-0.5))]="Down"
DE1%>%group_by(sig)%>%dplyr::summarise(count=n())

ex4f=ggplot()+
  geom_point_rast(data=DE1, aes(shape=TEAD_motif, x=logFC, y=-log10(FDR), color=sig), stroke=0.1, fill="white", size=0.8, alpha=0.6)+
  scale_color_manual(values = c("Down"="blue", "n.s."="black","Up"="red"), labels=c("Down"="Down\nn=247","Up"="Up\nn=56"), name="Sig.DE")+
  scale_shape_manual(values=c("Yes"=19, "No"=22), name="contain\nmotif")+
  coord_cartesian(ylim=c(0,3.5))+
  labs(color = NULL, title="Differential expression\nfrom TEAD2/4 KD")+
  theme1+theme(legend.position = c(0.75,0.7), legend.box = "horizontal")

pdf(paste0(path_fig4,"ex4f.volcano_KD_CAGE_TEAD24.pdf"), width = 1.7, height = 1.7)
print(ex4f)
dev.off()

ex4f2=ggplot()+
  geom_point_rast(data=DE1[which(DE1$sig %in% c("Up","Down") & DE1$TEAD_motif == "Yes" & DE1$k16ips %in% c("enhancer","promoter")),], aes(shape=TEAD_KD, x=logFC, y=-log10(FDR), color=k16ips), stroke=0.2, fill="white", size=1.8, alpha=0.9)+
  scale_color_manual(values = c("enhancer"="red", "promoter"="#B09C85FF"), labels=c("enhancer"="Enhancer","promoter"="Promoter"), name="Promoter type")+
  scale_shape_manual(values=c("Yes"=19, "No"=22), name="KD ATAC")+
  coord_cartesian(ylim=c(0,3.5))+
  labs(color = NULL, title="Differential expression\nfrom TEAD2/4 KD")+
  theme1+theme(legend.position = c(0.75,0.7), legend.box = "horizontal")
pdf(paste0(path_fig4,"ex4f2.volcano_KD_CAGE_TEAD24.pdf"), width = 1.7, height = 1.7)
print(ex4f2)
dev.off()

#==================
# fig ex4g
final3=read.delim(paste0(path_fig4_data,"ABC_summary_state.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

final3$state=factor(final3$state, levels=c("Active","Repressed","Primed","Others"))
final3$sample=factor(final3$sample, levels=c("iPSC","NSC","Neuron"))
final4=final3%>%group_by(sample, group)%>%dplyr::summarise(count=sum(count))

ex4f=ggplot() + 
  labs(x=NULL, y = "% of enhancer",title= "Chromatin state of ABC-linked enhancers", fill=NULL) +
  scale_fill_manual(values=c("Active"="#E64B35FF", "Primed"="#4DBBD5FF", "Repressed"="#00A087FF", "Others"="grey"))+
  scale_alpha_manual(values=c("All_enhancer"=0.5, "Predicted_enhancer"=0.9), guide=NULL)+
  facet_grid(cols=vars(sample))+
  scale_y_continuous(labels = scales::percent, limits=c(0, 1.2), breaks=c(0, 0.25,0.5,0.75,1))+
  geom_bar(data=final3, mapping=aes(y=percent, x=group, fill=state, alpha=group), linewidth=0.25, color = "black", stat="identity") + 
  geom_text(data=final4, mapping=aes(y=1.02, x=group, label=paste0("n = ",count)), size=1.8, angle=25, hjust=0, vjust=0)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig4,"ex4g.ABC_state.pdf"), width = 2, height = 1.7)
print(ex4g)
dev.off() 

#==========
# fig ex4h
final3b=read.delim(paste0(path_fig4_data,"ABC_summary_transcription.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

final3b$TSS=factor(final3b$TSS, levels=c("Transcribed","Non-transcribed"))
final3b$sample=factor(final3b$sample, levels=c("iPSC","NSC","Neuron"))
final4b=final3b%>%group_by(sample, group)%>%dplyr::summarise(count=sum(count))

ex4h=ggplot() + 
  labs(x=NULL, y = "% of enhancer",title= "Transcription of ABC-linked enhancers", fill=NULL) +
  scale_fill_manual(values=c("Transcribed"="#F39B7FFF", "Non-transcribed"="#8491B4FF"))+
  scale_alpha_manual(values=c("All_enhancer"=0.5, "Predicted_enhancer"=0.9), guide=NULL)+
  facet_grid(cols=vars(sample))+
  scale_y_continuous(labels = scales::percent, limits=c(0, 1.2), breaks=c(0, 0.25,0.5,0.75,1))+
  geom_bar(data=final3b, mapping=aes(y=percent, x=group, fill=TSS, alpha=group), linewidth=0.25, color = "black", stat="identity") + 
  geom_text(data=final4b, mapping=aes(y=1.02, x=group, label=paste0("n = ",count)), size=1.8, angle=25, hjust=0, vjust=0)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig4,"ex4h.ABC_transcription.pdf"), width = 2.2, height = 1.7)
print(ex4h)
dev.off() 



