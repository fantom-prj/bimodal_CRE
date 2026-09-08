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
library(cowplot)

###color#####
library("ggsci")
mypal = pal_npg("nrc", alpha = 1)(9)
mypal
library("scales")
show_col(mypal)

#####################
primary_folder=[primary_folder]
#primary_folder="/analysisdata/fantom6/Interactome/single_cell_wallace/Figure/"
path_fig2_data=paste0(primary_folder,"Fig2/data/")
path_fig2=paste0(primary_folder,"Fig2/out/")

setwd(path_fig2)

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

#=====================================
#Fig.ex2a
#tCRE vs aCRE
atCRE=read.delim(paste0(path_fig2_data,"separate_identification_aCRE_tCRE.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)
atCRE$label=gsub(" with","\nwith",atCRE$label)
atCRE1=atCRE%>%group_by(label, support)%>%dplyr::summarise(count=n())%>%mutate(percent=count/sum(count))

ex2a=ggplot() + 
  labs(x=NULL, y ="% of CRE", title="Detection of aCRE and tCRE\noverlapped by each other", fill="Overlap")+
  geom_bar(data=atCRE1, mapping=aes(x=label, y=percent, fill=support), stat="identity", linewidth=0.25, color="black", alpha=0.8)+
  scale_fill_manual(values=c("white", "grey"))+
  geom_text(data=atCRE1, mapping=aes(x=label, percent, label=count, group=support), position = position_stack(vjust = 0.5), size =1.8 )+
  scale_y_continuous(labels = scales::percent, limits=c(0,1), breaks=c(0,0.25,0.5,0.75,1))+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig2,"ex2a_taCRE.pdf"), width = 1.7, height = 1.7)
print(ex2a)
dev.off() 

#=====================================
#Fig.ex2b
#CRE features
aCRE=read.delim(paste0(path_fig2_data,"peaks.merged.all.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
aCRE=aCRE[which(nchar(aCRE$V1)<=5),] #only primary chromosome
aCRE$CRE_length=aCRE$V3-aCRE$V2
summary(aCRE$CRE_length)

aCRE%>%group_by(tCRE)%>%dplyr::summarise(median=median(CRE_length),count=n())
aCRE$CRE_length_bin=round(aCRE$CRE_length/100)*100
aCRE1=aCRE%>%group_by(tCRE,CRE_length_bin)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
aCRE1$label=paste0(signif(aCRE1$percent*100,2),"%")
aCRE1$label[which(aCRE1$percent<0.02)]=""
ex2b=ggplot() + 
  labs(x="Region length (nt)", y ="count", color="Transcription", title="Length of CREs (n = 463,866)")+
  geom_line(data=aCRE1, mapping=aes(x=CRE_length_bin, y=percent, color=tCRE), linewidth=0.25)+
  geom_vline(xintercept = 902, color = "black", linetype = "dotted", size=0.3)+
  scale_y_continuous(labels = percent)+
  scale_x_continuous(breaks=c(0,1000,2000,3000,4000), labels=c(0,1000,2000,3000,4000))+
  scale_color_npg()+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1), legend.position = c(0.7, 0.7))
pdf(paste0(path_fig2,"ex2b.length_CRE_ATAC.pdf"), width = 1.6, height = 1.7)
print(ex2b)
dev.off() 


#=====================================
#Fig.ex2d
#bimodal in metacell
data1=read.delim(paste0(path_fig2_data,"metacell_UMAP.tsv.gz"), header = T, stringsAsFactors = F, check.names = F)

data1$label=factor(data1$label, levels=cluster_info$label)
data1$variable=factor(data1$variable, levels=c("Gene","ATAC-signal","TSS-signal"))

ex2d=ggplot() + 
  labs(x="UMAP1", y ="UMAP2", title="UMAP of metacells using different signals", color=NULL)+
  facet_wrap(vars(variable), ncol=3, scales="free")+
  scale_color_npg()+
  geom_point(data=data1, mapping=aes(x=UMAP1, color=label, y=UMAP2), shape=16, alpha=0.75, size=0.45)+
  #geom_text_repel(data=umap_info, mapping=aes(x=UMAP1, y=UMAP2, label=label), size=1.8, color="black", bg.color = "white", bg.r=0.1)+
  theme1+theme(legend.position = "bottom",  panel.grid.major = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank())
pdf(paste0(path_fig2,"ex2d.RNA.aCRE.tCRE.umap.pdf"), width = 3.2, height = 1.7) #ext fig 2c
print (ex2d)
dev.off() 

#===================
# fig ex2f
both_table=read.delim(paste0(path_fig2_data,"CRE_463866_sample_expression_ATAC_TSS.plot.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
mboth_table=reshape2::melt(both_table, id=c(1:3,6:13))
mboth_table=mboth_table[which(mboth_table$sample_promoter_type %in% c("promoter-like","enhancer-like")),]
mboth_table$group3="Transcribed"
mboth_table$group3[which(mboth_table$group2=="Without transcription")]="Non-transcribed"
mboth_table$group3=factor(mboth_table$group3, levels=c("Transcribed","Non-transcribed"))
mboth_table$state=factor(mboth_table$state, levels=c("Active","Repressed","Primed","Bivalent"))

mboth_table2=mboth_table[which(mboth_table$sample_promoter_type == "enhancer-like" & mboth_table$state %in% c("Repressed","Primed")),]%>%group_by(variable,group3)%>%dplyr::summarise(p=wilcox.test(value ~state, alternative = "two.sided")$p.value)
mboth_table2$label="***"
mboth_table2$label[which(mboth_table2$p > 0.05)]="n.s."
mboth_table3=mboth_table[which(mboth_table$sample_promoter_type == "enhancer-like" & mboth_table$state %in% c("Repressed","Primed")),]%>%group_by(variable,group3,sample)%>%dplyr::summarise(p=wilcox.test(value ~state, alternative = "two.sided")$p.value)

a19a=ggplot() + 
  labs(x=NULL, y = "ATAC signal",title= NULL, fill=NULL) +
  coord_cartesian(ylim=c(0,25))+
  scale_fill_manual(values=c("Active"="#E64B35FF", "Primed"="#4DBBD5FF", "Repressed"="#00A087FF"))+
  facet_grid(cols=vars(group3), rows=vars(variable), scales="free", space = "free")+
  geom_boxplot(data=mboth_table[which(mboth_table$variable == "ATAC" & mboth_table$sample_promoter_type =="enhancer-like"),], mapping=aes(y=value, x=state, fill=state), size=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.8, outlier.shape = NA) + 
  geom_text(data=mboth_table2[which(mboth_table2$variable == "ATAC"),],mapping=aes(x=2.5, y=11, label=label), size=2.2)+
  annotate("segment", x=2,xend=3,y=10,yend=10, size=0.25)+
  theme1+theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(), legend.key.size = unit(0.3, 'cm'), legend.position = c(0.7,0.7))

a19b=ggplot() + 
  labs(x=NULL, y = "TSS signal",title= NULL, fill=NULL) +
  coord_cartesian(ylim=c(0,0.75))+
  scale_fill_manual(values=c("Active"="#E64B35FF", "Primed"="#4DBBD5FF", "Repressed"="#00A087FF"))+
  facet_grid(cols=vars(group3), rows=vars(variable), scales="free", space = "free")+
  geom_boxplot(data=mboth_table[which(mboth_table$variable == "TSS" & mboth_table$sample_promoter_type =="enhancer-like"),], mapping=aes(y=value, x=state, fill=state), size=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.8, outlier.shape = NA) + 
  geom_text(data=mboth_table2[which(mboth_table2$variable == "TSS"),],mapping=aes(x=2.5, y=0.45, label=label), size=2.2)+
  annotate("segment", x=2,xend=3,y=0.4,yend=0.4, size=0.25)+
  theme1+theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(), legend.key.size = unit(0.3, 'cm'), legend.position = c(0.7,0.7))

a19c=ggplot() + 
  labs(x=NULL, y = NULL,title= NULL, fill=NULL) +
  coord_cartesian(ylim=c(0,100))+
  scale_fill_manual(values=c("Active"="#E64B35FF", "Bivalent"="#00A087FF"))+
  facet_grid(cols=vars(group3), rows=vars(variable), scales="free", space = "free")+
  geom_boxplot(data=mboth_table[which(mboth_table$variable == "ATAC" & mboth_table$sample_promoter_type =="promoter-like"),], mapping=aes(y=value, x=state, fill=state), size=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.8, outlier.shape = NA) + 
  theme1+theme(axis.text.x = element_blank(),axis.ticks.x=element_blank() , legend.key.size = unit(0.3, 'cm'), legend.position = c(0.7,0.7))

a19d=ggplot() + 
  labs(x=NULL, y = NULL,title= NULL, fill=NULL) +
  coord_cartesian(ylim=c(0,125))+
  scale_fill_manual(values=c("Active"="#E64B35FF", "Bivalent"="#00A087FF"))+
  facet_grid(cols=vars(group3), rows=vars(variable), scales="free", space = "free")+
  geom_boxplot(data=mboth_table[which(mboth_table$variable == "TSS" & mboth_table$sample_promoter_type =="promoter-like"),], mapping=aes(y=value, x=state, fill=state), size=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.8, outlier.shape = NA) + 
  theme1+theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(), legend.key.size = unit(0.3, 'cm'), legend.position = c(0.7,0.7))

b19=grid.arrange(arrangeGrob(a19b,a19a, nrow = 2, heights=c(1,1)))
b20=grid.arrange(arrangeGrob(a19d,a19c, nrow = 2, heights=c(1,1)))

pdf(paste0(path_fig2,"ex2f.TSS_ATAC_signal_chromatin_state.pdf"), width = 3, height = 1.7)
print(plot_grid(b19, b20, ncol = 2, rel_widths = c(3,2.6)))
dev.off()

#===================
# fig ex2g
# correlation between aCRE and tCRE singal
CRE_rho=read.delim(paste0(path_fig2_data,"CRE_43697_rho_p_e.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

CRE_rho_summary=CRE_rho%>%group_by(promoter_type2)%>%dplyr::summarise(median_rho=median(rho),median_mlp=median(-log10(pv)),count=n())
k=wilcox.test(CRE_rho$rho[which(CRE_rho$promoter_type2 %in% c ("promoter-like","enhancer-like"))] ~CRE_rho$promoter_type2[which(CRE_rho$promoter_type2 %in% c ("promoter-like","enhancer-like"))], alternative = "two.sided")$p.value
CRE_rho$promoter_type2=factor(CRE_rho$promoter_type2,levels=c("promoter-like","enhancer-like","CTCF-alone","unclassed"))

ex2g=ggplot() + 
  labs(x=NULL, y = "Rho from transcribed CREs",title= "Correlation between ATAC\n& TSS signals") +
  coord_cartesian(ylim=c(0,0.8))+
  scale_fill_npg(guide=NULL)+
  ggdist::stat_halfeye(data=CRE_rho, mapping=aes(y=rho, x=promoter_type2, fill=promoter_type2),adjust = 0.5, width = 0.4, .width = 0, justification = -0.45, point_colour = NA, alpha=0.5)+
  geom_boxplot(data=CRE_rho, mapping=aes(y=rho, x=promoter_type2, fill=promoter_type2), linewidth=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.25, outlier.shape = NA) + 
  geom_text(data=CRE_rho_summary, mapping=aes(y=median_rho, x=promoter_type2, label=signif(median_rho,3)), angle=90, size=1.8, vjust=-1.3, hjust=0.5, color="black")+
  annotate("text",x=1.5, y=0.76, label="***", size=2.2)+
  annotate("segment", x=1,xend=2,y=0.75,yend=0.75, linewidth=0.25)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig2,"ex2g.CRE_43697_rho_p_e.pdf"), width = 1.4, height = 1.8) 
print(ex2g)
dev.off()

#==================
# fig ex2h
library(ComplexHeatmap)
library(magick)
library("RColorBrewer")

aCREe=read.delim(paste0(path_fig2_data,"peaks.merged.all.markalone_k16.subclass.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
aCREe$group="Others"
aCREe$group[which(aCREe$promoter_type2=="enhancer-like" & aCREe$tCRE == "Yes")]="Transcribed enhancer"
aCREe$group[which(aCREe$promoter_type2=="enhancer-like" & aCREe$tCRE == "No")]="Untranscribed enhancer"
aCREe$group[which(aCREe$promoter_type2=="promoter-like" & aCREe$tCRE == "Yes")]="Transcribed promoter"
aCREe$group[which(aCREe$promoter_type2=="promoter-like" & aCREe$tCRE == "No")]="Untranscribed promoter"

aCREe2=aCREe[which(aCREe$group != "Others"),c(3,5,7,10)]
colnames(aCREe2)[1:3]=c("iPSC","NSC","Neuron")
aCREe2[aCREe2=="NA"]=0
aCREe2[is.na(aCREe2)]=0
aCREe2[aCREe2=="CTCF"]=0
aCREe2$iPSC[intersect(grep("enhancer",aCREe2$iPSC),grep("promoter",aCREe2$group))]=1
aCREe2$NSC[intersect(grep("enhancer",aCREe2$NSC),grep("promoter",aCREe2$group))]=1
aCREe2$Neuron[intersect(grep("enhancer",aCREe2$Neuron),grep("promoter",aCREe2$group))]=1 #beige
aCREe2[aCREe2=="Primed_enhancer"]=2 #blue
aCREe2[aCREe2=="Repressed_enhancer"]=3 #green
aCREe2[aCREe2=="Active_enhancer"]=4 #red

aCREe2[aCREe2=="Promoter"]=4 #red
aCREe2[aCREe2=="Flanking_promoter"]=4 #red
aCREe2[aCREe2=="Bivalent_promoter"]=3 #green

aCREe2[,c(1:3)] <- aCREe2[,c(1:3)] %>% mutate_if(is.character, as.numeric)
aCREe3=aCREe2%>%group_by(group,iPSC,NSC,Neuron)%>%dplyr::summarise(count=n())
kk=aCREe2%>%group_by(group)%>%dplyr::summarise(count=n(),active_count=sum(rowSums(across(1:3) == 4) > 0, na.rm = TRUE))
kk$percent=kk$active_count/kk$count*100
write.table(kk,"heatmap.markalone_all_count.tsv", col.names=T, row.names=F, sep="\t", quote=F)
aCREe3$count=round(aCREe3$count/100)
aCREe4 <- aCREe3[rep(row.names(aCREe3), aCREe3$count), c(1:4)]

out3a=Heatmap(as.matrix(aCREe4[which(aCREe4$group=="Transcribed enhancer"),c(2:4)]),
              column_title = "Transcribed\nenhancer\n20,951",column_title_gp = gpar(fontsize = 7),
              column_names_gp = gpar(fontsize = 5), column_names_rot = 25, 
              col = c("0" = "white", "1" = "beige", "2" = "#4DBBD5FF", "3" = "#00A087FF", "4" = "#E64B35FF"),
              use_raster = T,
              cluster_columns = F,
              show_heatmap_legend = F,
              row_dend_reorder = FALSE,
              row_dend_width = unit(4, "mm"),
              row_dend_gp = gpar(lwd = 0.5))
out3b=Heatmap(as.matrix(aCREe4[which(aCREe4$group=="Untranscribed enhancer"),c(2:4)]),
              column_title = "Non-transcribed\nenhancer\n339,753",column_title_gp = gpar(fontsize = 7),
              column_names_gp = gpar(fontsize = 5), column_names_rot = 25, 
              col = c("0" = "white", "1" = "beige", "2" = "#4DBBD5FF", "3" = "#00A087FF", "4" = "#E64B35FF"),
              use_raster = T,
              cluster_columns = F,
              show_heatmap_legend = F,
              row_dend_reorder = FALSE,
              row_dend_width = unit(4, "mm"),
              row_dend_gp = gpar(lwd = 0.5))

out3c=Heatmap(as.matrix(aCREe4[which(aCREe4$group=="Transcribed promoter"),c(2:4)]),
              column_title = "Transcribed\npromoter\n22,341",column_title_gp = gpar(fontsize = 7),
              column_names_gp = gpar(fontsize = 5), column_names_rot = 25, 
              col = c("0" = "white", "1" = "beige", "2" = "#4DBBD5FF", "3" = "#00A087FF", "4" = "#E64B35FF"),
              use_raster = T,
              cluster_columns = F,
              show_heatmap_legend = F,
              row_dend_reorder = FALSE,
              row_dend_width = unit(4, "mm"),
              row_dend_gp = gpar(lwd = 0.5))

out3d=Heatmap(as.matrix(aCREe4[which(aCREe4$group=="Untranscribed promoter"),c(2:4)]),
              column_title = "Non-transcribed\npromoter\n12,094",column_title_gp = gpar(fontsize = 7),
              column_names_gp = gpar(fontsize = 5), column_names_rot = 25, 
              col = c("0" = "white", "1" = "beige", "2" = "#4DBBD5FF", "3" = "#00A087FF", "4" = "#E64B35FF"),
              use_raster = T,
              cluster_columns = F,
              show_heatmap_legend = F,
              row_dend_reorder = FALSE,
              row_dend_width = unit(4, "mm"),
              row_dend_gp = gpar(lwd = 0.5))

legend_labels = c("Active","Repressive/Bivalent","Primed","Context-Specific Enhancers")
legend_colors = c("4" = "#E64B35FF", "3" = "#00A087FF", "2" = "#4DBBD5FF","1" = "beige")

global_legend = Legend(
  labels = legend_labels,
  legend_gp = gpar(fill = legend_colors),
  direction = "vertical",
  ncol = 4,
  labels_gp = gpar(fontsize = 5))

g1 = grid.grabExpr(draw(out3a))
g2 = grid.grabExpr(draw(out3b))
g3 = grid.grabExpr(draw(out3c))
g4 = grid.grabExpr(draw(out3d))
gl = grid.grabExpr(draw(global_legend))
plot4 = plot_grid(g1, g2, g3, g4, ncol = 4, rel_widths = c(1, 1,1,1))
ex2h = plot_grid(plot4, gl, nrow = 2, rel_heights = c(1.5, 0.1))

# Display
grid.newpage()
grid.draw(ex2h)

pdf(paste0(path_fig2,"ex2h.heatmap.chromatin_state_all.pdf"), width=4.5, height=1.95)
grid.draw(ex2h)
dev.off()

#===========
#ex2i
#genomic location of CRE
scACREv3=read.delim(paste0(path_fig2_data,"all_aCRE.PE.final.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
scACREv3$group2="Transcribed"
scACREv3$group2[which(scACREv3$tCRE == "No")]="Non-\ntranscribed"
scACREv3$group2=factor(scACREv3$group2, levels=c("Transcribed","Non-\ntranscribed"))
data3=scACREv3[which(scACREv3$analysis == "included"),]%>%group_by(group2,promoter_type_CT,genomic_region)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
data3$genomic_region=factor(data3$genomic_region, levels=c("5'end","5'UTR","3'UTR","Exon","Intron","Intergenic"))
data4=data3%>%group_by(group2,promoter_type_CT)%>%dplyr::summarise(count=paste0("n=",sum(count)))

ex2i=ggplot() + 
  facet_wrap(vars(group2), ncol=2, scales="free_y")+
  labs(x=NULL, y ="% of CRE", title="Location of CREs", fill="Annotated\nfeature")+
  geom_bar(data=data3, mapping=aes(x=promoter_type_CT, y=percent, fill=genomic_region), stat="identity", linewidth=0.25, color="black", alpha=0.8)+
  scale_fill_manual(values=c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF" ,"white"))+
  geom_text(data=data4, mapping=aes(x=promoter_type_CT, y=1.02, label=count), angle=90, hjust=(0), vjust=(0.5), size =1.8 )+
  scale_y_continuous(labels = scales::percent, limits=c(0,1.4), breaks=c(0, 0.25, 0.5, 0.75,1))+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig2,"ex2i.Annotated_feature_CRE.pdf"), width = 2.2, height = 2)
print(ex2i)
dev.off() 

#=======================
#fig ex2j TATA-INR
dd=read.delim(paste0(path_fig2_data,"TATA_INR_distance.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
dd95 = dd[dd$prob >= 0.95 * max(dd$prob), ]
range(dd95$distance) #25-33

ex2j=ggplot(dd, aes(distance, prob)) +
  geom_line(linewidth=0.2)+
  geom_vline(xintercept=c(25,33), linetype="dashed", linewidth=0.2)+
  labs(x = "TATA-Inr distance (nt)",y = "Predicted probability", title="Predicted distance\nbetween TATA and Inr") +
  theme1
pdf(paste0(path_fig2,"ex2j.distance_TATA_INR.pdf"), width = 1.4, height = 1.6)
print(ex2j)
dev.off() 

#=======================
# fig ex2k
# fig2g promoter-like CRE that linked to a ENSG
aCRE_cell_both2=read.delim(paste0(path_fig2_data,"aCRE_cell_Reg_State_FE.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
aCRE_cell_both2=aCRE_cell_both2[which(aCRE_cell_both2$promoter_type_CT == "promoter-like"),]
aCRE_cell_both2$cell=factor(aCRE_cell_both2$cell, levels=c("iPSC","NSC","Neuron"))

aCRE_cell_both2$label="n.s."
aCRE_cell_both2$label[which(aCRE_cell_both2$p.val<0.05)]="*"
aCRE_cell_both2$label[which(aCRE_cell_both2$p.val<0.01)]="**"
aCRE_cell_both2$label[which(aCRE_cell_both2$p.val<0.001)]="***"
aCRE_cell_both2$label=factor(aCRE_cell_both2$label, levels=c("n.s.","*","**","***"))
aCRE_cell_both2$logOR=log(aCRE_cell_both2$OR)
aCRE_cell_both2$group2=gsub("Non-transcribed","Non-\ntranscribed",aCRE_cell_both2$group2)
aCRE_cell_both2$group2=factor(aCRE_cell_both2$group2, levels=c("Transcribed","Non-\ntranscribed"))

ex2k=ggplot() + 
  labs(x=NULL, y =NULL ,title= "Enrichment between promoter\nstates & regulatory elements", size="p_value") +
  facet_grid(cols=vars(group2),rows=vars(cell))+
  scale_y_discrete(limits=rev)+
  geom_point(data=aCRE_cell_both2, mapping=aes(y=State, x=Reg, fill=logOR, size=label),shape=21)+
  scale_fill_gradient2(low = "#3C5488FF", mid = "white", high = "#DC0000FF", midpoint = 0)+ 
  scale_size_manual(values=c(0.4,1.2,2.4))+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig2,"ex2k.regE_promoter_chromatin_state.FE.pdf"), width = 1.8, height = 2)
print(ex2k)
dev.off() 

#================
#fig ex2l
kk=read.delim(paste0(path_fig2_data,"heatmap.markalone_ZBTB6_KLF13.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
kk=kk[which(kk$value == "TRUE"),]
kk=kk[which(kk$sample != "Neuron"),]
kk=kk[which(kk$state %in% c("Primed","All","Repressed")),]
kk$state=factor(kk$state,levels=c("Primed","All","Repressed"))
kk$group2=gsub(" enhancer","",kk$group2)
kk$group2=gsub("Non-transcribed","Non-\ntranscribed",kk$group2)
kk$group2=factor(kk$group2, levels=c("Transcribed","Non-\ntranscribed"))
ex2l=ggplot() + 
  facet_grid(cols=vars(group2), rows=vars(sample))+
  labs(x=NULL, y ="% of enhancer", title="Transition to active", fill="Scope")+
  geom_bar(data=kk, mapping=aes(x=state, y=percent, fill=state), stat="identity", linewidth=0.25, color="black", alpha=0.8)+
  geom_text(data=kk, mapping=aes(x=state, y=1, label=paste0(signif(percent,3),"%")), angle=90, hjust=(0), vjust=(0.5), size =1.8 )+
  #geom_text(data=kk, mapping=aes(x=state, y=percent+3, label=paste0("n=",active_count)), angle=90, hjust=(0), vjust=(0.5), size =1.8 )+
  scale_y_continuous(limits=c(0,100))+
  scale_fill_manual(values=c("All"="grey", "Primed"="#4DBBD5FF", "Repressed"="#00A087FF"), guide=NULL)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig2,"ex2l.iPSNSC_tCRE_active_percentage.pdf"), width = 1.5, height = 2)
print(ex2l)
dev.off() 




