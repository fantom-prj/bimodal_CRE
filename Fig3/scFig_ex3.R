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
# fig ex3a -> schematics

#===============================================================================
# fig ex3b
data2=read.delim(paste0(path_fig3_data,"variability_ChromVar_noFP.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data2=data2[which(data2$set1 %in% c("Enhancer","Promoter")),]
data2$set1=factor(data2$set1, levels=c("Enhancer","Promoter"))
data2$V2[which(data2$rank >5)]=NA

ex3b=ggplot(data2[which(data2$set1 %in% c("Enhancer","Promoter")),], aes(x = rank, y = variability)) +
  facet_grid(cols=vars(set1), rows=vars(set2), scales = "free") +
  coord_cartesian(xlim=c(1,879)) +
  geom_errorbar(aes(ymin = bootstrap_lower_bound, ymax = bootstrap_upper_bound), linewidth = 0.2, color = "grey25") +
  geom_point(shape=20, size = 0.05, color = "black") +
  geom_text_repel(aes(label=V2), direction = "y", segment.size = 0.2 , nudge_x=250, hjust = 0, size=1.6, force=0.1) +
  labs(x = "Sorted 879 TF motifs", y = "Variability", title="Variability of TF motif activity") +
  theme1+ theme(axis.text.x = element_blank(), axis.ticks = element_blank())
pdf(paste0(path_fig3,"ex3b.variability_ChromVar_noFP.pdf"), width = 2.6, height = 2.3)
print (ex3b)
dev.off() 
#-> no cutoff, show all variability

#=========================
# fig ex3c
# plot x-y variability and correlation

together2=read.delim(paste0(path_fig3_data,"auto_correlation_to_p_en.z.summarise_noFP.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
outputfile1=read.delim(paste0(path_fig3_data,"TF_grouping_by_chromvar_correlation_and_literature.n247.tsv.gz"),  header=T, stringsAsFactors = F, check.names = F)

# together2a=unique(together2[,c("group1", "group2","median_zscaled_var")])
# choose the strongest correlation if there is two TFs
together2=together2%>%group_by(group1,group2,motifID)%>%slice_max(abs(acf_exVall_opt))

ex3c=ggplot() + 
  labs(x="ChromVar Variability (Scaled)", y ="Abs(Correlation)", alpha="Pass chromVar\n& Correlation", title="Evaluable TF mmotifs")+
  facet_grid(rows=vars(group1), cols=vars(group2))+
  geom_point(data=together2, mapping=aes(x=zscaled_var, color=interaction(group1, group2), y=abs(acf_exVall_opt),  alpha = pass_both), size=0.2, shape=19)+
  scale_color_npg(guide=NULL)+
  scale_alpha_manual(values=c("Yes"=1, "No"=0.3))+
  geom_hline(yintercept=0.5, linetype="dashed", linewidth=0.2)+
  geom_vline(data=together2a,aes(xintercept=median_zscaled_var), linetype="dashed", linewidth=0.2)+
  theme1+theme(legend.position = "right")
pdf(paste0(path_fig3,"ex3c.scale_variability_n_correlation_xy_noFP.pdf"), width = 2.6, height = 1.8)
print (ex3c)
dev.off() 

#======================
# fig ex3d
# plotting zscore heatmap showing all cell types (n=10)
ddata=read.delim(paste0(path_fig3_data,"heatmap_z_247motif.plot.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
cluster_info=read.delim(paste0(primary_folder,"Fig1/cluster_info.tsv"),header=T, stringsAsFactors = F, check.names = F)

mddata=data.frame(spread(ddata[,c(4,7,9)], key=2 , value=3))

rownames(mddata)=mddata$V2
together2b=mddata[,c(2:41)]

cn=unique(ddata[,c(5,6,7)])
cn$set=factor(cn$set,levels=c("en_aCRE","en_tCRE","p_aCRE","p_tCRE"))
cn$label=factor(cn$label, levels=cluster_info$label)
cn=cn[order(cn$label, cn$set),]

together2b=together2b[,cn$label2]
data2= read.delim(paste0(path_fig3_data,"variability_ChromVar_noFP.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data2$V2[which(data2$rank >5)]=NA

#varfile=read.delim(paste0(path_fig3,"variability_for_all_analyses.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
#varfile$V2[which(varfile$rank >5)]=NA
need=unique(varfile$label[!is.na(data2$V2)])
need=need[!need %in% c("MA1962.1_POU2F1::SOX2","MA0754.3_CUX1","MA1642.2_NEUROG2","MA0755.2_CUX2")]
need=c(need,"MA0809.3_TEAD4","MA0138.3_REST","MA1109.2_NEUROD1","MA0514.3_Sox3")
rows_to_label=unique(need)

library(dendextend)
library(ggdendro)

row_dist <- dist(together2b)
row_clust <- hclust(row_dist)
row_order <- row_clust$labels[row_clust$order]

dend_row <- as.dendrogram(row_clust)
p_rowtree <- ggdendrogram(dend_row, rotate = T)+ 
  scale_y_reverse()+ 
  theme_void()+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
p_rowtree$layers[[2]]$aes_params$size <- 0.2

ddata$V2 <- factor(ddata$V2, levels = row_order)
ddata$set=factor(ddata$set,levels=c("en_aCRE","en_tCRE","p_aCRE","p_tCRE"))
ddata$label=factor(ddata$label, levels=cluster_info$label)

ddata$motifName=sapply(strsplit(as.character(ddata$V2),"_"),"[",2)
ddata$motifName[-which(ddata$V2 %in% rows_to_label)]=""
ddata$motifName[which(ddata$cluster != "OPC_like")]=""
ddata$motifName[which(ddata$set != "p_tCRE")]=""

anno_df <- ddata %>%
  distinct(set, label) %>%
  mutate(anno1 = ifelse(grepl("tCRE", set), "TSS", "ATAC"), anno2 = ifelse(grepl("en", set), "Enhancer", "Promoter"))
anno_df$anno1a=substr(anno_df$anno1,1,1)
anno_df$anno2a=substr(anno_df$anno2,1,1)

library(ggrepel)
library(ggnewscale)
top_y <- length(levels(ddata$V2)) + 2.5

out01=ggplot(ddata, aes(x = set, y = V2, fill = scaledz)) +
  labs(x=NULL, y=NULL, fill="Z-score")+
  facet_grid(cols=vars(label))+
  geom_tile_rast() +
  #geom_text_repel(aes(label = motifName), nudge_x = 5, direction = "y", hjust = 0, segment.size = 0.2 , size=1.6) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick3", midpoint = 0) +
  new_scale_fill() +
  geom_tile(data = anno_df, aes(x = set, y = -5.5, fill = anno1), height = 5,inherit.aes = FALSE)+
  geom_text(data = anno_df, aes(x = set, y = -5.5, label=anno1a), size=1.4,inherit.aes = FALSE)+
  scale_fill_manual(values = c("ATAC" = "green", "TSS" = "deepskyblue"), name="Signal") +
  new_scale_fill() +
  geom_tile(data = anno_df, aes(x = set, y = -10.5, fill = anno2), height = 5,inherit.aes = FALSE)+
  geom_text(data = anno_df, aes(x = set, y = -10.5, label=anno2a), size=1.4,inherit.aes = FALSE)+
  scale_fill_manual(values = c("Enhancer"= "#7E6148FF", "Promoter"= "#B09C85FF"), name="CRE type")+
  new_scale_fill() +
  geom_tile(data = anno_df, aes(x = set, y = top_y, fill = label), height = 5,inherit.aes = FALSE)+
  scale_fill_manual(values = c("iPSC_Sphase"="#E64B35FF" , "iPSC_G2M"="#4DBBD5FF", "NSC_stem"="#00A087FF", "NSC_1"="#3C5488FF", 
                               "NSC_2"="#F39B7FFF", "NPC_like"="#8491B4FF","Immature_neuron"="#91D1C2FF", "Mature_neuron_1"="#DC0000FF","Mature_neuron_2"="#7E6148FF", "OPC_like"="#B09C85FF"), name="Cell-type") +
  theme1 +
  theme(panel.grid.major = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), legend.position="right", strip.background =element_blank(), strip.text.x.top = element_blank())

pdf(paste0(path_fig3,"ex3d.heatmap_z_247motif.pdf"), width=4, height=3.6)
print(plot_grid(p_rowtree, out01, ncol = 2, rel_widths = c(1, 8)))
dev.off()

label_data <- ddata %>%filter(motifName != "")
out02=ggplot(ddata, aes(x = set, y = V2, fill = scaledz)) +
  labs(x=NULL, y=NULL, fill="Z-score")+
  geom_tile_rast() +
  geom_text_repel(aes(x=4.78, label = motifName), nudge_x = 3, direction = "y", hjust = 1, segment.size = 0.2 , size=1.6) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick3", midpoint = 0) +
  geom_segment(data=label_data, aes(x=4.5, xend=4.8, y=V2, yend=V2), linewidth=0.2) +
  #facet_grid(cols=vars(label))+
  new_scale_fill() +
  geom_tile(data = anno_df, aes(x = set, y = -5.5, fill = anno1), height = 5,inherit.aes = FALSE)+
  geom_text(data = anno_df, aes(x = set, y = -5.5, label=anno1a), size=1.4,inherit.aes = FALSE)+
  scale_fill_manual(values = c("ATAC" = "green", "TSS" = "deepskyblue"), name="Signal") +
  new_scale_fill() +
  geom_tile(data = anno_df, aes(x = set, y = -10.5, fill = anno2), height = 5,inherit.aes = FALSE)+
  geom_text(data = anno_df, aes(x = set, y = -10.5, label=anno2a), size=1.4,inherit.aes = FALSE)+
  scale_fill_manual(values = c("Enhancer"= "#7E6148FF", "Promoter"= "#B09C85FF"), name="CRE type")+
  new_scale_fill() +
  geom_tile(data = anno_df, aes(x = set, y = top_y, fill = label), height = 5,inherit.aes = FALSE)+
  scale_fill_manual(values = c("iPSC_Sphase"="#E64B35FF" , "iPSC_G2M"="#4DBBD5FF", "NSC_stem"="#00A087FF", "NSC_1"="#3C5488FF", 
                               "NSC_2"="#F39B7FFF", "NPC_like"="#8491B4FF","Immature_neuron"="#91D1C2FF", "Mature_neuron_1"="#DC0000FF","Mature_neuron_2"="#7E6148FF", "OPC_like"="#B09C85FF"), name="Cell-type") +
  theme1 +
  theme(panel.grid.major = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), legend.position="right", strip.background =element_blank(), strip.text.x.top = element_blank())
pdf(paste0(path_fig3,"ex3d2.heatmap_z_247motif.label.pdf"), width=2.3, height=3.6)
print(out02)
dev.off()

#====================
# fig ex3e
# FE with with all cutoffs
FEsum3=read.delim(paste0(path_fig3_data,"TF_grouping_by_chromvar_correlation_and_literature.FE.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

FEsum3a=FEsum3[which(FEsum3$variable.y != "TF_closer" ),]
FEsum3a$variable.y = gsub("TF_pioneer","Pioneer",FEsum3a$variable.y)
FEsum3a$variable.y = gsub("TF_activator","Activator",FEsum3a$variable.y)
FEsum3a$variable.y = gsub("TF_repressor","Repressor",FEsum3a$variable.y)
FEsum3a$variable.y = factor(FEsum3a$variable.y, levels=c("Pioneer","Activator","Repressor"))
FEsum3a$variable.x = factor(FEsum3a$variable.x, levels=c("ATAC+","ATAC-","TSS+","TSS-"))
FEsum3a$var_cutoff=factor(FEsum3a$var_cutoff, levels=c("Var > P25","Var > Median","Var > P75"))

mx <- max(FEsum3a$logOR[is.finite(FEsum3a$logOR)])
mn <- min(FEsum3a$logOR[is.finite(FEsum3a$logOR)])
FEsum3a$logOR[FEsum3a$logOR == Inf]  <- mx + 0.1
FEsum3a$logOR[FEsum3a$logOR == -Inf] <- mn - 0.1

ex3e=ggplot() + 
  labs(x="From literature", y ="From co-expression", title="TF action enrichment", fill="logOR")+
  facet_grid(rows=vars(var_cutoff),cols=vars(corr_cutoff))+
  scale_y_discrete(limits=rev)+
  geom_tile(data=FEsum3a, mapping=aes(x=variable.y, fill=logOR, y=variable.x), size=0.6)+
  scale_fill_gradient2(low="#3C5488FF",mid="white",high="#F0AB00FF", midpoint=0)+
  theme1+theme(panel.grid.major = element_blank(),legend.position = "top", legend.justification = "right",
               axis.ticks = element_blank(), axis.line = element_blank(),
               axis.text.x = element_text(hjust=1, vjust=1, angle=25))
pdf(paste0(path_fig3,"ex3e.TF_action_enrichment_allcutoffs.pdf"), width = 4.2, height = 2.6)
print (ex3e)
dev.off() 

#======================
# fig ex3f
# upsetter
outputfile1=read.delim(paste0(path_fig3_data,"TF_grouping_by_joint_pen_chromvar_correlation_and_literature.n181.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

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
ex3f=ggplot(outputfile3, aes(x = result)) +
  geom_bar(aes(fill=coupled), linewidth=0.25, color="black", alpha=1,  width=0.8,  stat="count") + 
  scale_fill_manual(values=c("Yes"="black", "No"="white"))+
  scale_x_upset()+
  labs(x="", y = "Number of TF motifs", title="TF regulatory modes from\npromoter enhancer joint analysis", fill="Coupled") +
  theme1+theme(legend.position = c(0.8,0.7))+
  theme_combmatrix(combmatrix.panel.point.color.fill = "black",
                   combmatrix.panel.point.size = 1.5,
                   combmatrix.panel.line.size = 0.5,
                   combmatrix.label.text = element_text(color ="black", size=6),
                   combmatrix.label.extra_spacing = 0.3,
                   combmatrix.label.make_space = FALSE)
pdf(paste0(path_fig3,"ex3f.joint_pen_upseter.TF.mode.by.expression.correlation.pdf"), width = 2.1, height = 1.7)
print(ex3f)
dev.off()

#=====================
# fig ex3g
data6=read.delim(paste0(path_fig3_data,"ATAC_RNA.venn.with.cuttoff_joint_pen.tsv.gz"), header=T, stringsAsFactors=F, check.names = F)

data6$uncoupled=data6$union-data6$coupled
data6$var_cutoff=factor(data6$var_cutoff, levels=c("Var > P25","Var > Median","Var > P75"))
mdata6=reshape2::melt(data6, id=c(1,2,3,5))
mdata6$modal_specific=paste0(signif(mdata6$modal_specific*100,2),"%")
mdata6$modal_specific[which(mdata6$variable == "uncoupled")]=NA

ex3g=ggplot(mdata6, aes(x = as.character(corr_cutoff), y=value, fill=variable)) +
  geom_bar(linewidth=0.25, alpha=1,  width=0.8,  color="black", stat="identity") + 
  facet_grid(cols=vars(var_cutoff))+
  scale_fill_manual(values=c("coupled"="black", "uncoupled"="white"))+
  geom_text(data=mdata6, aes(label=modal_specific, y=5),angle=90, hjust=0, size=1.8,color="black")+
  labs(x="Cutoff of correlation", y = "Number of TF motifs", title="Modality-specific TF motifs identified\nfrom promoter enhancer joint motif analysis", fill=NULL) +
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1), legend.position = c(0.8,0.8))
pdf(paste0(path_fig3,"ex3g.joint_pen_ATAC_RNA.venn.with.cuttoff.pdf"), width = 3.2, height = 1.7)
print(ex3g)
dev.off()

#=====================
# fig ex3h
# FE with literature only considering those with results from both literature and correlation  
FEsum3=read.delim(paste0(path_fig3_data,"TF_grouping_by_joint_pen_chromvar_correlation_and_literature.FE.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

FEsum3a=FEsum3[which(FEsum3$variable.y != "TF_closer" & FEsum3$var_cutoff == "Var > Median" & FEsum3$corr_cutoff == 0.5),]
FEsum3a$variable.y = gsub("TF_pioneer","Pioneer",FEsum3a$variable.y)
FEsum3a$variable.y = gsub("TF_activator","Activator",FEsum3a$variable.y)
FEsum3a$variable.y = gsub("TF_repressor","Repressor",FEsum3a$variable.y)
FEsum3a$variable.y = factor(FEsum3a$variable.y, levels=c("Pioneer","Activator","Repressor"))
FEsum3a$variable.x = factor(FEsum3a$variable.x, levels=c("ATAC+","ATAC-","TSS+","TSS-"))
ex3h=ggplot() + 
  labs(x="From literature", y ="From co-expression", title="TF action enrichment", fill="logOR")+
  scale_y_discrete(limits=rev)+
  geom_tile(data=FEsum3a, mapping=aes(x=variable.y, fill=logOR, y=variable.x), size=0.6)+
  scale_fill_gradient2(low="#3C5488FF",mid="white",high="#F0AB00FF", midpoint=0)+
  theme1+theme(panel.grid.major = element_blank(),legend.position = "right",
               axis.ticks = element_blank(), axis.line = element_blank(),
               axis.text.x = element_text(hjust=1, vjust=1, angle=25))
pdf(paste0(path_fig3,"ex3h.joint_pen_TF_action_enrichment122.pdf"), width = 1.5, height = 1.7)
print (ex3h)
dev.off() 

#=================

