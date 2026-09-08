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
chromvar_folder=paste0(primary_folder,"Data_and_code/ChromVar/")
path_fig3_data=paste0(primary_folder,"Fig3/data/")
path_fig3=paste0(primary_folder,"Fig3/out/")

setwd(path_fig3)

#=====================
theme1=theme(  panel.background = element_rect(fill = "white", color = "white", linewidth = 0.25, linetype = "solid"),
               panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid', color = "grey90"), 
               axis.text.x = element_text(color ="black"),
               axis.text.y = element_text(color ="black"),
               plot.title = element_text(hjust = 0.5,  color ="black", margin = margin(0.1,0.1,0.1,0.1, "cm")),
               plot.subtitle = element_text(color="black", hjust=0.5),
               text = element_text(size=6),
               strip.text = element_text(size=6),
               legend.background = element_rect(fill="white", color="black", size=0.25),
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
# fig 3a
# brief venn, variability and correlation
together2=read.delim(paste0(path_fig3_data,"auto_correlation_to_p_en.z.summarise_noFP.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
together2=together2%>%group_by(group1,group2,motifID)%>%slice_max(abs(acf_exVall_opt))
#this file pass expression level only, some of them (21) pass neither variability nor correlation
#247 of them pass both at the same condition (enhancer ATAC, enhancer RNA, promoter ATAC, promoter RNA)
#26 of them pass both across different conditions

together2a=together2
together2a$ID=paste0(together2a$motifID,"_",together2a$motifName)

Corr=unique(together2a$ID[which(abs(together2a$acf_exVall_opt) >= 0.5)])  #324
Var=unique(together2a$ID[which(together2a$chromVar_Var == "yes")]) #338
nboth=unique(together2a$ID[which(abs(together2a$acf_exVall_opt) >= 0.5 & together2a$chromVar_Var == "yes")]) #247
nboth_any=intersect(Corr, Var) #273

library(VennDiagram)
grid.newpage()
f3a <- draw.pairwise.venn(area1 = length(Corr), area2 = length(Var), cross.area = length(nboth_any), euler.d = TRUE, scaled = T,
                          category=c("Correlation","ChromVar"), alpha=c(0.5,0.5),
                          fill = c("#4DBBD5FF","#E64B35FF"), lty = "blank", cex = 0.7, cat.cex = 0.75, cat.col = "black")
pdf(paste0(path_fig3,"f3a.Corr_Var.venn.pdf"), width = 1, height = 1)
grid.draw(f3a)
dev.off()

#===============================================================================
# fig 3b
# plot variability heatmap with sign(correlation)

together2=read.delim(paste0(path_fig3_data,"auto_correlation_to_p_en.z.summarise_noFP.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
together2=together2%>%group_by(group1,group2,motifID)%>%slice_max(abs(acf_exVall_opt))
together2a=together2[which(abs(together2$acf_exVall_opt) >= 0.5 & together2$chromVar_Var == "yes"),]

together2a$ID=paste0(together2a$motifID,"_",together2a$motifName)

#together2a=together2a%>%group_by(group1, group2)%>%dplyr::mutate(scale_var_sign=chromVar_variability*sign(acf_exVall_opt))
together2a=together2a%>%group_by(group1, group2)%>%dplyr::mutate(scale_var_sign=zscaled_var*sign(acf_exVall_opt))

together2a$group1[which(together2a$group1 == "ATAC-signal")]="ATAC"
together2a$group1[which(together2a$group1 == "TSS-signal")]="TSS"
together2a$group=paste0(together2a$group1,"_",together2a$group2)
together2c=data.frame(spread(together2a[,c(20,22,21)], key=2, value=3))

rownames(together2c)=together2c$ID
together2c=together2c[,c(2:5)]
together2c[is.na(together2c)]=0

library(dendextend)
library(ggdendro)

row_dist <- dist(together2c)
row_clust <- hclust(row_dist)
row_order <- row_clust$labels[row_clust$order]

dend_row <- as.dendrogram(row_clust)
p_rowtree <- ggdendrogram(dend_row, rotate = T)+ 
  scale_y_reverse()+ 
  scale_x_reverse()+ 
  theme_void()+
  theme(plot.margin = unit(c(0.66, 0, 0.66, 0), "cm"))
p_rowtree$layers[[2]]$aes_params$size <- 0.2

together2c$ID=rownames(together2c)
together2b=reshape2::melt(together2c, id=5)

together2b$ID <- factor(together2b$ID, levels = rev(row_order))
together2b$variable=factor(together2b$variable,levels=c("ATAC_Enhancer","TSS_Enhancer","ATAC_Promoter","TSS_Promoter"))

rows_to_label=c("Pou5f1::Sox2","TEAD4","TCF12","NEUROD1","ONECUT2","ZEB1","MEIS2","Sox3","SOX13")
together2b$motifName=sapply(strsplit(as.character(together2b$ID),"_"),"[",2)
together2b$motifName[-which(together2b$motifName %in% rows_to_label)]=""
together2b$motifName[which(together2b$variable != "TSS_Promoter")]=""

library(ggrepel)
library(ggnewscale)
library(cowplot)
label_data <- together2b %>%dplyr::filter(motifName != "")

f3b=ggplot() +
  labs(x=NULL, y=NULL, fill="Scaled_variability\nx sign(Corr)", title="247 TF motif variability from 4 contexts")+
  geom_tile(data=together2b, mapping=aes(x = variable, y = ID, fill = value)) +
  geom_text_repel(data=together2b, aes(x=4.67,  y = ID, label = motifName), hjust = 0, nudge_x = 2.4, nudge_y = 20, direction = "y", force_pull = 0, segment.size = 0.2 , size=1.6) +
  geom_segment(data=label_data, aes(x=4.5, xend=4.7, y=ID, yend=ID), linewidth=0.2) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick3", midpoint = 0) +
  theme1 +
  theme(panel.grid.major = element_blank(), axis.text.y = element_blank(),axis.text.x = element_text(angle=25, hjust=1, vjust=1), axis.ticks = element_blank(), axis.line = element_blank(), legend.position="right")

pdf(paste0(path_fig3,"f3b.pheatmap.var.ChromVar.score.p_en_var_247TFmotif.pdf"), width=2.8, height=2.6)
print(plot_grid(p_rowtree, f3b, ncol = 2, rel_widths = c(1, 8)))
dev.off()

#===============================================================================
# fig 3c
# dot plot #z-score version
overall=read.delim(paste0(path_fig3_data,"plot.table.all.noFP.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
overall2=overall%>%group_by(group1, group2, variable,motifName,T2)%>%dplyr::summarise(value2=mean(value), TF1=mean(expression_TF1), TF2=mean(expression_TF2)) 
overall2=overall2 %>%group_by(variable)  %>%dplyr::mutate (rel.TF1=TF1*0.5/max(TF1), rel.TF2=TF2*0.5/max(TF2))

selection=c("Pou5f1::Sox2","TEAD4","FOS::JUND","NEUROD1","ONECUT2","ZEB1","MEIS3","Sox3")
overall1=overall2[which(overall2$motifName %in% selection),]
overall1$motifName=factor(overall1$motifName, levels=c("Pou5f1::Sox2","TEAD4","FOS::JUND","NEUROD1","ONECUT2","INSM1","ZEB1","MEIS3","Sox3"))
f3c= ggplot(overall1) +
  geom_point(mapping=aes(x =T2, y = value2, color=group2), size=0.2, shape=19) + 
  stat_smooth(mapping=aes(x =T2, y = rel.TF1*20), method = "lm", formula = y ~ poly(x, 21), se = FALSE, color="black", alpha=1, linewidth=0.2) +
  stat_smooth(mapping=aes(x =T2, y = rel.TF2*20), method = "lm", formula = y ~ poly(x, 21), se = FALSE, color="grey", alpha=0.6, linewidth=0.2) +
  scale_color_manual(values=c("#7E6148FF","#B09C85FF"))+
  facet_grid(row=vars(group1), col=vars(motifName), scales = "free_x")+
  labs(color=NULL, y="Relative TF motif activity & relative TF expression", x = "Ranked by pseudotime for trajectory 1", title = "TF motif activity by chromatin accessibility and transcription") +
  theme1+theme(axis.text.x = element_blank())
pdf(paste0(path_fig3,"f3c.Chromvar.Pseudotime.z.neuronal.trajectory.aCRE.vs.tCRE.expression.pdf"), width = 6, height = 2.5)
f3c
dev.off()

#trajectory2
boverall=read.delim(paste0(path_fig3_data,"plot.table.all.trajectory2.noFP.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
boverall2=boverall%>%group_by(group1, group2, variable,motifName,T2)%>%dplyr::summarise(value2=mean(value), TF1=mean(expression_TF1), TF2=mean(expression_TF2)) 
boverall2=boverall2 %>%group_by(variable)  %>%dplyr::mutate (rel.TF1=TF1*0.5/max(TF1), rel.TF2=TF2*0.5/max(TF2))

selection=c("TCF7L2","Hand1","FOS::JUND*","NEUROD1*","BCL11A","CREB3L1*","REST*","GATA6","OSR1*")
boverall1=boverall2[which(boverall2$motifName %in% selection),]
boverall1$motifName=factor(boverall1$motifName, levels=c("TCF7L2","Hand1","BCL11A","GATA6"))

xi= ggplot(boverall1) +
  #scale_y_continuous(limits=c(-0.4, 1.05), breaks=c(0,0.5,1))+
  geom_point(mapping=aes(x =T2, y = value2, color=group2), size=0.2, shape=19) + 
  stat_smooth(mapping=aes(x =T2, y = rel.TF1*20), method = "lm", formula = y ~ poly(x, 11), se = FALSE, color="black", alpha=1, linewidth=0.2) +
  #stat_smooth(mapping=aes(x =T2, y = rel.TF2*20), method = "lm", formula = y ~ poly(x, 11), se = FALSE, color="grey", alpha=0.6, size=0.2) +
  scale_color_manual(values=c("#7E6148FF","#B09C85FF"))+
  scale_x_continuous(limits=c(43,84))+
  facet_grid(row=vars(group1), col=vars(motifName), scales = "free_x")+
  labs(color=NULL, y="Mean TF motif activity & TF relative expression", x = "Ranked by pseudotime for trajectory 2 (excluding iPSC)", title = "TF motif activity by chromatin accessibility and transcription") +
  theme1+theme(axis.text.x = element_blank())
pdf(paste0(path_fig3,"f3c2.Chromvar.Pseudotime.z.OPC.trajectory.aCRE.vs.tCRE.expression.pdf"), width = 3.5, height = 2.5)
xi
dev.off()

#===============================================================================
# fig 3d
# main fig with selected TF
ddata1=read.delim(paste0(path_fig3_data,"heatmap_z_12motif.plot.tsv.gz") , header=T, stringsAsFactors = F, check.names = F)
ddata1$V2=sapply(strsplit(as.character(ddata1$V2),"_"),"[",2)
need=c("SOX13","ONECUT2","NEUROG2","NEUROD1","JUND","Jun","CREB3L1","NFYC","TFDP1","TEAD4","ZEB1","Pou5f1::Sox2")
cluster_info=read.delim(paste0(primary_folder,"Fig1/cluster_info.tsv"),header=T, stringsAsFactors = F, check.names = F)

ddata1$V2 <- factor(ddata1$V2, levels = need)
ddata1$set=factor(ddata1$set,levels=c("en_aCRE","en_tCRE","p_aCRE","p_tCRE"))
ddata1$label=factor(ddata1$label, levels=cluster_info$label)

anno_df <- ddata %>%
  distinct(set, label) %>%
  mutate(anno1 = ifelse(grepl("tCRE", set), "TSS", "ATAC"), anno2 = ifelse(grepl("en", set), "Enhancer", "Promoter"))
anno_df$anno1a=substr(anno_df$anno1,1,1)
anno_df$anno2a=substr(anno_df$anno2,1,1)

library(ggnewscale)
top_y <- length(levels(ddata1$V2)) + 0.5

out01=ggplot(ddata1, aes(x = set, y = V2, fill = scaledz)) +
  labs(x=NULL, y=NULL, fill="Scaled z-score", title="TF activity across 10 cell types")+
  facet_grid(cols=vars(label))+
  geom_tile_rast() +
  scale_y_discrete(position = "right") +
  scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick3", midpoint = 0) +
  new_scale_fill() +
  geom_tile(data = anno_df, aes(x = set, y = -0.1, fill = anno1), height = 0.5,inherit.aes = FALSE)+
  geom_text(data = anno_df, aes(x = set, y = -0.1, label=anno1a), size=1.4,inherit.aes = FALSE)+
  scale_fill_manual(values = c("ATAC" = "green", "TSS" = "deepskyblue"), name="Signal") +
  new_scale_fill() +
  geom_tile(data = anno_df, aes(x = set, y = -0.6, fill = anno2), height = 0.5,inherit.aes = FALSE)+
  geom_text(data = anno_df, aes(x = set, y = -0.6, label=anno2a), size=1.4,inherit.aes = FALSE)+
  scale_fill_manual(values = c("Enhancer"= "#7E6148FF", "Promoter"= "#B09C85FF"), name="CRE type")+
  new_scale_fill() +
  geom_tile(data = anno_df, aes(x = set, y = top_y, fill = label), height = 0.5,inherit.aes = FALSE)+
  scale_fill_manual(values = c("iPSC_Sphase"="#E64B35FF" , "iPSC_G2M"="#4DBBD5FF", "NSC_stem"="#00A087FF", "NSC_1"="#3C5488FF", 
                               "NSC_2"="#F39B7FFF", "NPC_like"="#8491B4FF","Immature_neuron"="#91D1C2FF", "Mature_neuron_1"="#DC0000FF","Mature_neuron_2"="#7E6148FF", "OPC_like"="#B09C85FF"), name="Cell-type") +
  theme1 +
  theme(panel.spacing = unit(0.12, "lines"), panel.grid.major = element_blank(), axis.text.x = element_blank(),  axis.ticks = element_blank(), axis.line = element_blank(), legend.position="right", strip.background =element_blank(), strip.text.x.top = element_blank())
pdf(paste0(path_fig3,"f3d.heatmap_z_12motif.pdf"), width=3.3, height=1.7)
print(out01)
dev.off()

#===============================================================================
# fig 3e
overall=read.delim(paste0(path_fig3_data,"plot.table.all.noFP.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
overall1=overall[which(overall$group1 == "ATAC-signal"),]
overall2=overall[which(overall$group1 == "TSS-signal"),]
overall.all=left_join(overall1[,c("metacell","variable","group2","value")],overall2[,c("metacell","variable","group2","value")],by=c("metacell","variable","group2"),suffix=c("_TSS","_ATAC"),copy=F)
overall.alla=overall.all%>%group_by(variable,group2)%>%dplyr::summarise(pearson_r=cor.test(value_TSS,value_ATAC,alternative="greater",method="pearson")$estimate,
                                                                        pearson_p=cor.test(value_TSS,value_ATAC,alternative="greater",method="pearson")$p.value)
n247=read.delim(paste0(path_fig3_data,"heatmap_z_247motif.plot.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
overall.allb=overall.alla[which(overall.alla$variable %in% unique(n247$variable)),]
overall.allb$group2=gsub("Enhancer","Transcribed\nenhancer",overall.allb$group2)
overall.allb$group2=gsub("Promoter","Transcribed\npromoter",overall.allb$group2)

f3e=ggplot() + 
  labs(x=NULL, y = "Pearson r (247 TF motifs)",title= "Correlation between motif activities\nderived from ATAC & TSS signal") +
  scale_fill_npg(guide=NULL)+
  ggdist::stat_halfeye(data=overall.allb, mapping=aes(y=pearson_r, x=group2, fill=group2),adjust = .5, width = .4, .width = 0, justification = -.4, point_colour = NA, alpha=0.5)+
  geom_boxplot(data=overall.allb, mapping=aes(y=pearson_r, x=group2, fill=group2), size=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.2, outlier.shape = NA) + 
  geom_hline(yintercept=0.3, linetype="dashed", linewidth=0.2)+
  theme1
pdf(paste0(path_fig3,"f3e.corr_TFmotif_ATAC_TSS_zscore.pdf"), width = 1.6, height = 1.7)
print(f3e)
dev.off()

#===============================================================================
# fig 3f
# upsetter
outputfile1=read.delim(paste0(path_fig3_data,"TF_grouping_by_chromvar_correlation_and_literature.n247.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

outputfile3=outputfile1[,c("label","motifID","ATAC+","ATAC-","TSS+","TSS-")]
for (i in 1:nrow(outputfile3)){
  number=which(outputfile3[i,c(3:6)] == "Yes")
  outputfile3$result[i]=list(colnames(outputfile3)[number+2])}
options(scipen=999)
outputfile3$coupled="No"
outputfile3$coupled[which(outputfile3$`ATAC+`=="Yes" & outputfile3$`TSS+` == "Yes")]="Yes"
outputfile3$coupled[which(outputfile3$`ATAC-`=="Yes" & outputfile3$`TSS-` == "Yes")]="Yes"

library(ggupset)
library(ComplexUpset)
f3f=ggplot(outputfile3, aes(x = result)) +
  geom_bar(aes(fill=coupled), linewidth=0.25, color="black", alpha=1,  width=0.8,  stat="count") + 
  scale_fill_manual(values=c("Yes"="black", "No"="white"))+
  scale_x_upset()+
  labs(x="", y = "Number of TF motifs", title="TF regulatory modes", fill="Coupled") +
  theme1+theme(legend.position = c(0.8,0.7))+
  theme_combmatrix(combmatrix.panel.point.color.fill = "black",
                   combmatrix.panel.point.size = 1.5,
                   combmatrix.panel.line.size = 0.5,
                   combmatrix.label.text = element_text(color ="black", size=6),
                   combmatrix.label.extra_spacing = 0.3,
                   combmatrix.label.make_space = FALSE)
pdf(paste0(path_fig3,"f3f.upseter.TF.mode.by.expression.correlation.pdf"), width = 1.9, height = 1.7)
print(f3f)
dev.off()

#upsetter side -> venn diagram is better
motif_ATAC=together2b1$motifID[which(together2b1$`ATAC+` == "Yes" | together2b1$`ATAC-` == "Yes")]
motif_RNA=together2b1$motifID[which(together2b1$`TSS+` == "Yes" | together2b1$`TSS-` == "Yes")]
both=intersect(motif_ATAC,motif_RNA)
both_real=together2b1$motifID[which(together2b1$coupled == "Yes")]

library(VennDiagram)
grid.newpage()
f3f2 <- draw.pairwise.venn(area1 = length(motif_ATAC), area2 = length(motif_RNA), cross.area = length(both), euler.d = TRUE, scaled = T,
                           category=c("ATAC","TSS"), alpha=c(0.5,0.5),
                           fill = c("green","deepskyblue"), lty = "blank", cex = 0.7, cat.cex = 0.75, cat.col = "black")
pdf(paste0(path_fig3,"f3f2.ATAC_RNA.venn.pdf"), width = 1, height = 1)
grid.draw(f3f2)
dev.off()

#===============================================================================
# fig 3g
data6=read.delim(paste0(path_fig3_data,"ATAC_RNA.venn.with.cuttoff.tsv.gz"), header=T, stringsAsFactors=F, check.names = F)
data6$uncoupled=data6$union-data6$coupled
data6$var_cutoff=factor(data6$var_cutoff, levels=c("Var > P25","Var > Median","Var > P75"))
mdata6=reshape2::melt(data6, id=c(1,2,3,5))
mdata6$modal_specific=paste0(signif(mdata6$modal_specific*100,2),"%")
mdata6$modal_specific[which(mdata6$variable == "uncoupled")]=NA

f3g=ggplot(mdata6, aes(x = as.character(corr_cutoff), y=value, fill=variable)) +
  geom_bar(linewidth=0.25, alpha=1,  width=0.8,  color="black", stat="identity") + 
  facet_grid(cols=vars(var_cutoff))+
  scale_fill_manual(values=c("coupled"="black", "uncoupled"="white"))+
  geom_text(data=mdata6, aes(label=modal_specific, y=5),angle=90, hjust=0, size=1.8,color="black")+
  #scale_y_continuous(labels = scales::percent)+
  labs(x="Cutoff of correlation", y = "Number of TF motifs", title="Modality-specific TF motif identification", fill=NULL) +
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1), legend.position = c(0.8,0.8))
pdf(paste0(path_fig3,"f3g.ATAC_RNA.venn.with.cuttoff.pdf"), width = 2.8, height = 1.7)
print(f3g)
dev.off()

#===============================================================================
# fig 3h
# FE with literature only considering those with results from both literature and correlation  
FEsum3=read.delim(paste0(path_fig3_data,"TF_grouping_by_chromvar_correlation_and_literature.FE.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

FEsum3a=FEsum3[which(FEsum3$variable.y != "TF_closer" & FEsum3$var_cutoff == "Var > Median" & FEsum3$corr_cutoff == 0.5),]
FEsum3a$variable.y = gsub("TF_pioneer","Pioneer",FEsum3a$variable.y)
FEsum3a$variable.y = gsub("TF_activator","Activator",FEsum3a$variable.y)
FEsum3a$variable.y = gsub("TF_repressor","Repressor",FEsum3a$variable.y)
FEsum3a$variable.y = factor(FEsum3a$variable.y, levels=c("Pioneer","Activator","Repressor"))
FEsum3a$variable.x = factor(FEsum3a$variable.x, levels=c("ATAC+","ATAC-","TSS+","TSS-"))
f3h=ggplot() + 
  labs(x="From literature", y ="From co-expression", title="TF action enrichment", fill="logOR")+
  scale_y_discrete(limits=rev)+
  geom_tile(data=FEsum3a, mapping=aes(x=variable.y, fill=logOR, y=variable.x), size=0.6)+
  scale_fill_gradient2(low="#3C5488FF",mid="white",high="#F0AB00FF", midpoint=0)+
  theme1+theme(panel.grid.major = element_blank(),legend.position = "right",
               axis.ticks = element_blank(), axis.line = element_blank(),
               axis.text.x = element_text(hjust=1, vjust=1, angle=25))
pdf(paste0(path_fig3,"f3h.TF_action_enrichment167.pdf"), width = 1.5, height = 1.7)
print (f3h)
dev.off() 

#===============================================================================
# get number of TF from 247 motif
outputfile=read.delim(paste0(path_fig3_data,"variance_corr_p_en_400motif.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
outputfile1=read.delim(paste0(path_fig3_data,"TF_grouping_by_chromvar_correlation_and_literature.n247.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
outputfilea=outputfile[which(outputfile$motif %in% outputfile1$motif),]
length(unique(outputfilea$TF_symbol))#214

#===============================================================================


