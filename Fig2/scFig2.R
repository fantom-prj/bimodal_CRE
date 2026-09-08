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

#===============================================================================
#fig 2a
state=c("Active enhancer","CTCF","Promoter","Promoter flank","Promoter","Primed enhancer","Bivalent promoter","Repressed region",
"Genomic region","Repressed enhancer","Genomic region","Genomic region","Primed enhancer","Active enhancer","Active enhancer","Genomic region")
state_number=16:1
state_df=data.frame(cbind(state_number, state))
state_df$group="Others"
state_df$group[grep("nhancer",state_df$state)]="Enhancers"
state_df$group[grep("romoter",state_df$state)]="fPromoters"
state_df=state_df[order(state_df$group,state_df$state),]
state_df$group[grep("romoter",state_df$state)]="Promoters"

# plot overall figure #all cell types
emission1=read.delim(paste0(path_fig2_data,"chromatin_state_emission.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
emission1$V4=factor(emission1$V4, levels=c("ATAC","CTCF","K27ac","K27m3","K4m1","K4m3"))
emission1$V2=factor(emission1$V2, levels=rev(state_df$state_number))
z21=ggplot(emission1, aes(y=V2, x=V4, fill=V6)) + 
  labs(y="Chromatin state cluster", x =NULL, title="Emission parameter", color=NULL)+
  geom_tile(alpha=0.8)+
  geom_text(mapping=aes(label=label), size=1.8)+
  scale_fill_gradient(low = "white", high = "#3C5488FF")+
  theme1+theme(panel.grid.major = element_blank(),  axis.text.x = element_text(color ="black", angle = 90, hjust=1, vjust=0.5),
               legend.position = "none", axis.text.y =element_blank(),
               plot.margin = unit(c(0.1, 0.1, 1.21, 0.25), "cm"))

mregion2=read.delim(paste0(path_fig2_data,"chromatin_state_TSS.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
mregion2$cell=factor(mregion2$cell, levels=c("iPS","NSC","NRN"))
mregion2$variable=factor(mregion2$variable, levels=rev(state_df$state_number))
z4=ggplot(mregion2, aes(y=variable, x=coord, fill=value)) + 
  labs(y=NULL, x =NULL, title="RefSeq TSS", fill=NULL)+
  geom_tile()+
  facet_wrap(vars(cell), ncol=3, strip.position="right")+
  scale_fill_gradient(low = "white", high = "#E64B35FF")+
  theme1+ theme(panel.grid.major = element_blank(),   axis.text.x = element_text(color ="black", angle = 90, hjust=1, vjust=0.5), axis.text.y =element_blank(),
                legend.position = "none", strip.background = element_blank(), strip.text.y = element_blank(),
                plot.margin = unit(c(0.1, 0.1, 1.27, 0.25), "cm"))

content2=read.delim(paste0(path_fig2_data,"chromatin_state_genome.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
content2$cell=factor(content2$cell, levels=c("iPS","NSC","NRN"))
content2$V4=factor(content2$V4, levels=rev(state_df$state_number))
z22=ggplot(content2, aes(y=as.factor(V4), x=cell, fill=cover)) + 
  labs(y=NULL, x =NULL, title="Genome", color=NULL)+
  geom_tile()+
  geom_text(mapping=aes(label=label), size=1.8)+
  scale_fill_gradient(low = "white", high = "#00A087FF", limits = c(0, 0.022))+
  theme1+theme(panel.grid.major = element_blank(),  axis.text.x = element_text(color ="black", angle = 90, hjust=1, vjust=0.5),
               legend.position = "none", axis.text.y =element_blank(),
               plot.margin = unit(c(0.1, 0.1, 1.393, 0.25), "cm"))

feature2=read.delim(paste0(path_fig2_data,"chromatin_state_feature.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
feature2$cell=factor(feature2$cell, levels=c("iPS","NSC","NRN"))  
feature2$V4=factor(feature2$V4, levels=rev(state_df$state_number))
z23=ggplot(feature2, aes(y=V4, x=variable, fill=value)) + 
  labs(y=NULL, x =NULL, title="Feature enrichment", color=NULL)+
  geom_tile()+
  #geom_text(mapping=aes(label=label), size=2.2)+
  facet_wrap(vars(cell), ncol=3, strip.position="right")+
  scale_fill_gradient(low = "white", high = "#7E6148FF")+
  theme1+theme(panel.grid.major = element_blank(),  axis.text.x = element_text(color ="black", angle = 90, hjust=1, vjust=0.5), axis.text.y =element_blank(),
               strip.background = element_blank(), strip.text.y = element_blank(), legend.position = "none",
               plot.margin = unit(c(0.1, 0.1, 0.2, 0.25), "cm"))

pdf(paste0(path_fig2,"f2a.chromatin_state.pdf"), width = 6.8, height = 2.5)
grid.arrange(z21,z4,z22,z23,ncol=4, widths=c(2,1.8,1.6,4.3))
dev.off() 

#================
# fig 2b
library(ComplexHeatmap)
library(magick)
library("RColorBrewer")

scACREv3=read.delim(paste0(path_fig2_data,"all_aCRE.PE.final.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
scACREv4=scACREv3[which(scACREv3$analysis == "included"),c("iPSC_CT","NSC_CT","Neuron_CT","tCRE")]
colnames(scACREv4)=c("iPSC","NSC","Neuron","TSS")
scACREv4[is.na(scACREv4)]=1
scACREv4[scACREv4=="CTCF"]=1
scACREv4$iPSC[grep("nhancer",scACREv4$iPSC)]=2
scACREv4$NSC[grep("nhancer",scACREv4$NSC)]=2
scACREv4$Neuron[grep("nhancer",scACREv4$Neuron)]=2
scACREv4$iPSC[grep("romoter",scACREv4$iPSC)]=3
scACREv4$NSC[grep("romoter",scACREv4$NSC)]=3
scACREv4$Neuron[grep("romoter",scACREv4$Neuron)]=3
scACREv4[,c(1:3)] <- scACREv4[,c(1:3)] %>% mutate_if(is.character, as.numeric)

taCRE=scACREv4[which(scACREv4$TSS == "Yes"),c(1:3)]
taCRE1=taCRE%>%group_by(iPSC,NSC,Neuron)%>%dplyr::summarise(count=n())
taCRE1$count=round(taCRE1$count/20) #take 1/20 for clustering
taCRE2 <- taCRE1[rep(row.names(taCRE1), taCRE1$count), c(1:3)]
taCRE3=taCRE2
taCRE3$anno <- apply(taCRE3, 1, max)

group_col <- c("1" = "white", "2" = "#7E6148FF", "3" = "#B09C85FF")
ha_right <- rowAnnotation(
  Group = anno_simple(taCRE3$anno, col = group_col),
  width = unit(1, "mm"),
  annotation_name_gp = gpar(fontsize = 0) )

f2b1=Heatmap(as.matrix(taCRE2),
        column_title = "43,697 CREs\nwith transcription", column_title_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 5), column_names_rot = 25, 
        col = colorRampPalette(c("white", "#7E6148FF", "#B09C85FF"))(3),
        use_raster = T,
        cluster_columns = F,
        right_annotation = ha_right,
        show_heatmap_legend = F,
        row_dend_reorder = FALSE,
        row_dend_width = unit(5, "mm") )

pdf(paste0(path_fig2,"f2b.heatmap.43697taCRE.pdf"), width=1, height=2.2, useDingbats = FALSE)
draw(f2b1)
dev.off()

#
aaCRE=scACREv4[which(scACREv4$TSS == "No"),c(1:3)]
aaCRE1=aaCRE%>%group_by(iPSC,NSC,Neuron)%>%dplyr::summarise(count=n())
aaCRE1$count=round(aaCRE1$count/200) #take 1/200 for clustering
aaCRE2 <- aaCRE1[rep(row.names(aaCRE1), aaCRE1$count), c(1:3)]
aaCRE3=aaCRE2
aaCRE3$anno <- apply(aaCRE3, 1, max)

ha2_right <- rowAnnotation(
  Group = anno_simple(aaCRE3$anno, col = group_col),
  width = unit(1, "mm"),
  annotation_name_gp = gpar(fontsize = 0) )

f2b2=Heatmap(as.matrix(aaCRE2),
             column_title = "420,169 CREs\nwithout transcription", column_title_gp = gpar(fontsize = 7),
             column_names_gp = gpar(fontsize = 5), column_names_rot = 25, 
             col = colorRampPalette(c("white", "#7E6148FF", "#B09C85FF"))(3),
             use_raster = T,
             cluster_columns = F,
             right_annotation = ha2_right,
             show_heatmap_legend = F,
             row_dend_reorder = FALSE,
             row_dend_width = unit(5, "mm") )
pdf(paste0(path_fig2,"f2b.heatmap.420169aaCRE.pdf"), width=1, height=2.2, useDingbats = FALSE)
draw(f2b2)
dev.off()

#=====================
# fig 2c
library(ComplexHeatmap)
library(magick)
library("RColorBrewer")

aCREe=read.delim(paste0(path_fig2_data,"peaks.merged.all.markalone_k16.subclass.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
aCREe1=aCREe[which(aCREe$promoter_type2=="enhancer-like"),]
aCREe2=aCREe1[,c(3,5,7)]
colnames(aCREe2)=c("iPSC","NSC","Neuron")
aCREe2[aCREe2=="NA"]=0
aCREe2[is.na(aCREe2)]=0
aCREe2[aCREe2=="CTCF"]=0
aCREe2[aCREe2=="Repressed_enhancer"]=1 #green
aCREe2[aCREe2=="Primed_enhancer"]=2 #blue
aCREe2[aCREe2=="Active_enhancer"]=3 #red
aCREe2 <- aCREe2 %>% mutate_if(is.character, as.numeric)
aCREe3=aCREe2%>%group_by(iPSC,NSC,Neuron)%>%dplyr::summarise(count=n())
aCREe3$count=round(aCREe3$count/100) #take 1/100 for clustering
aCREe4 <- aCREe3[rep(row.names(aCREe3), aCREe3$count), c(1:3)]

aCREe2%>%group_by(iPSC)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
aCREe2%>%group_by(NSC)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
aCREe2%>%group_by(Neuron)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))

f2c=Heatmap(as.matrix(aCREe4),
             column_title = "360,704 enhancer-\nlike CREs", column_title_gp = gpar(fontsize = 7),
             column_names_gp = gpar(fontsize = 5), column_names_rot = 25, 
             col = colorRampPalette(c("white","#00A087FF","#4DBBD5FF", "#E64B35FF"))(4),
             use_raster = T,
             cluster_columns = F,
             show_heatmap_legend = T,
             row_dend_reorder = FALSE,
             row_dend_width = unit(5, "mm"), 
             heatmap_legend_param = list(
               title = NULL,
               labels_gp = gpar(fontsize = 5),
               at = c("1", "2", "3"),
               labels = c("Repressive", "Primed", "Active"),
               grid_height = unit(2, "mm"),
               grid_width = unit(2, "mm"),
               direction = "vertical"))

pdf(paste0(path_fig2,"f2c.heatmap_360704eaCRE.pdf"), width=1.1, height=2.2)
draw(f2c,heatmap_legend_side = "bottom")
dev.off()

#===================
# fig 2d
# FE use the ChromHMM result directly
summaryCRE=read.delim(paste0(path_fig2_data,"all_aCRE.TSS.access.markeralone_cuttag.FEresult.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
summaryCRE$variable=factor(summaryCRE$variable, levels=c("SE","Bivalent","Repressed","Primed","Active"))
summaryCRE$label="n.s."
summaryCRE$label[which(summaryCRE$p.val<0.05)]="*"
summaryCRE$label[which(summaryCRE$p.val<0.01)]="**"
summaryCRE$label[which(summaryCRE$p.val<0.001)]="***"
summaryCRE$label=factor(summaryCRE$label, levels=c("n.s.","*","**","***"))
summaryCRE$logOR=log(summaryCRE$OR)
summaryCRE$cell=factor(summaryCRE$cell, levels=c("iPSC","NSC","Neuron"))

f2d=ggplot() + 
  labs(x=NULL, y =NULL ,title= "Enrichment with transcription ", size="p_value") +
  facet_grid(cols=vars(promoter_type_CT))+
  geom_point(data=summaryCRE, mapping=aes(y=variable, x=cell, fill=logOR, size=label),shape=21)+
  scale_fill_gradient2(low = "#3C5488FF", mid = "white", high = "#DC0000FF", midpoint = 0)+ 
  scale_size_manual(values=c(1.2,2.4))+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig2,"f2d.transcription_chromatin_state_SE.FE.pdf"), width = 2, height = 1.7)
print(f2d)
dev.off() 

#==========
# fig 2e
GC=read.delim(paste0(path_fig2_data,"both_GCcontent_result.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
GC1=GC[which(GC$promoter_type %in% c("promoter-like","enhancer-like") & GC$bimodal %in% c("transcrib_midpoint","untranscrib_midpoint")),]
GC1$bimodal=gsub("untranscrib_midpoint","Non-transcribed",GC1$bimodal)
GC1$bimodal=gsub("transcrib_midpoint","Transcribed",GC1$bimodal)
GC1$bimodal=factor(GC1$bimodal, levels=c("Transcribed","Non-transcribed"))

f2e = ggplot()+
  scale_color_npg()+
  geom_line(data = GC1, aes(x=as.numeric(position), y=GCprecent/100, color=group, group=group),alpha=0.75, linewidth=0.25)+
  facet_grid( rows=vars(promoter_type), cols=vars(bimodal))+
  labs(color=NULL, x=NULL, y=NULL, title="GC content")+
  scale_y_continuous(labels = scales::percent)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig2,"f2e.GCcontent.aCRE.pdf"), width = 2.1, height = 1.7)
print(f2e)
dev.off()

#=====================================
#fig 2f

scACREv3=read.delim(paste0(path_fig2_data,"all_aCRE.PE.final.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
scACREv3$Reg_Ele="Null"
scACREv3$Reg_Ele[which(scACREv3$CGI == "Yes")]="CGI"
scACREv3$Reg_Ele[which(scACREv3$TATA_box == "Yes")]="TATA"
scACREv3$Reg_Ele[which(scACREv3$CGI == "Yes" & scACREv3$TATA_box == "Yes")]="Both"

scACREv3$group2="Transcribed"
scACREv3$group2[which(scACREv3$tCRE == "No")]="Non-transcribed"
data1=scACREv3[which(scACREv3$analysis == "included"),]%>%group_by(promoter_type_CT,group2, Reg_Ele)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
data1$Reg_Ele=factor(data1$Reg_Ele, levels=c("CGI","Both","TATA","Null"))
data1$promoter_type_CT=factor(data1$promoter_type_CT, levels=c("promoter-like","enhancer-like","CTCF-alone","unclassed"))
data1$group2=factor(data1$group2, levels=c("Transcribed","Non-transcribed"))
data2=data1%>%group_by(promoter_type_CT,group2)%>%dplyr::summarise(count=paste0("n=",sum(count)))

f2f=ggplot() + 
  facet_wrap(vars(group2), ncol=2, scales="free_y")+
  labs(x=NULL, y ="% of CRE", title="Regulatory elements from CREs")+
  geom_bar(data=data1, mapping=aes(x=promoter_type_CT, y=percent, fill=Reg_Ele), stat="identity", linewidth=0.25, color="black", alpha=0.8)+
  geom_text(data=data2, mapping=aes(x=promoter_type_CT, y=1.02, label=count), angle=90, hjust=(0), vjust=(0.5), size =1.8 )+
  scale_y_continuous(labels = scales::percent, limits=c(0,1.3), breaks=c(0,0.25,0.5,0.75,1))+
  scale_fill_manual(values=c("#4DBBD5FF","darkorchid3","#E64B35FF","white"))+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig2,"f2f.regulatory_element_class_tCRE.pdf"), width = 2.2, height = 2)
print(f2f)
dev.off() 

#======================================
# fig 2g
aCRE_cell_both2=read.delim(paste0(path_fig2_data,"aCRE_cell_Reg_State_FE.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
aCRE_cell_both2=aCRE_cell_both2[which(aCRE_cell_both2$promoter_type_CT == "enhancer-like"),]
aCRE_cell_both2$cell=factor(aCRE_cell_both2$cell, levels=c("iPSC","NSC","Neuron"))

aCRE_cell_both2$label="n.s."
aCRE_cell_both2$label[which(aCRE_cell_both2$p.val<0.05)]="*"
aCRE_cell_both2$label[which(aCRE_cell_both2$p.val<0.01)]="**"
aCRE_cell_both2$label[which(aCRE_cell_both2$p.val<0.001)]="***"
aCRE_cell_both2$label=factor(aCRE_cell_both2$label, levels=c("n.s.","*","**","***"))
aCRE_cell_both2$logOR=log(aCRE_cell_both2$OR)
aCRE_cell_both2$group2=factor(aCRE_cell_both2$group2, levels=c("Transcribed","Non-transcribed"))

f2g=ggplot() + 
  labs(x=NULL, y =NULL ,title= "Enrichment between enhancer\nstates & regulatory elements", size="p_value") +
  facet_grid(cols=vars(group2),rows=vars(cell))+
  scale_y_discrete(limits=rev)+
  geom_point(data=aCRE_cell_both2[which(aCRE_cell_both2$promoter_type_CT == "enhancer-like"),], mapping=aes(y=State, x=Reg, fill=logOR, size=label),shape=21)+
  scale_fill_gradient2(low = "#3C5488FF", mid = "white", high = "#DC0000FF", midpoint = 0)+ 
  scale_size_manual(values=c(0.4,1.2,2.4))+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig2,"f2g.regE_enhancer_chromatin_state.FE.pdf"), width = 2, height = 2)
print(f2g)
dev.off() 

#====================
# f2h
enricha=read.delim(paste0(path_fig2_data,"F_enrichemnt_for_xT.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
# from Fig2.motif_enrichment/motif_enrich.R

enrichbscale=data.frame(spread(enricha[,c(2,10,12)], key=2, value=3))
rownames(enrichbscale)=enrichbscale$ALT_ID
enrichbscale=enrichbscale[,-1]
row_order <- rownames(enrichbscale)[hclust(dist(enrichbscale))$order]

enricha$ALT_ID=factor(enricha$ALT_ID,levels=row_order)
enricha$Set5=paste0(enricha$Set3,"_",enricha$Set2)
enricha$Set5=gsub("_en","\nenhancer",enricha$Set5)
enricha$Set5=gsub("_p","\npromoter",enricha$Set5)
enricha$Set5=factor(enricha$Set5, levels=c("active\npromoter","bivalent\npromoter","active\nenhancer","repressed\nenhancer","primed\nenhancer"))
enricha$Set1=factor(enricha$Set1, levels=c("iPSC","NSC","NRN"))
f2h=ggplot() + 
  facet_wrap(~Set5, ncol=5)+
  labs(x=NULL, y =NULL, title="Motif enrichment of non-transcribed CREs", fill="Scaled\nenrichment")+
  scale_y_discrete(position = "right", limits=rev)+
  geom_tile(data=enricha, mapping=aes(x=Set1, fill=scaled_value, y=ALT_ID), size=0.6)+
  scale_fill_gradient2(low="white",high="firebrick3")+
  theme1+theme(panel.grid.major = element_blank(),
               axis.ticks = element_blank(), axis.line = element_blank(),
               axis.text.x = element_text(color ="black",angle=90,hjust=1, vjust=0.5),
               axis.text.y = element_text(color ="black", size=5))
pdf(paste0(path_fig2,"f2h.heatmap_TF_enrichemnt_for_xT.pdf"), width = 2.4, height = 2)
print (f2h)
dev.off() 

#===================
# fig 2i
kk=read.delim(paste0(path_fig2_data,"heatmap.markalone_ZBTB6_KLF13.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
kk=kk[which(kk$value == "TRUE"),]
kk=kk[which(kk$sample == "Neuron"),]
kk=kk[which(kk$state %in% c("Primed","All","Repressed","Repressed_ZBTB6")),]
kk$state=factor(kk$state,levels=c("Primed","All","Repressed","Repressed_ZBTB6"))
kk$group2=gsub(" enhancer","",kk$group2)
kk$group2=gsub("Non-transcribed","Non-\ntranscribed",kk$group2)
kk$group2=factor(kk$group2, levels=c("Transcribed","Non-\ntranscribed"))
f2i=ggplot() + 
  facet_grid(cols=vars(group2))+
  labs(x=NULL, y ="% of enhancer in neuron", title="Transition to active", fill="Scope")+
  geom_bar(data=kk, mapping=aes(x=state, y=percent, fill=state), stat="identity", linewidth=0.25, color="black", alpha=0.8)+
  geom_text(data=kk, mapping=aes(x=state, y=1, label=paste0(signif(percent,3),"%")), angle=90, hjust=(0), vjust=(0.5), size =1.8 )+
  geom_text(data=kk, mapping=aes(x=state, y=percent+3, label=paste0("n=",active_count)), angle=90, hjust=(0), vjust=(0.5), size =1.8 )+
  scale_y_continuous(limits=c(0,100))+
  scale_fill_manual(values=c("All"="grey", "Primed"="#4DBBD5FF", "Repressed"="#00A087FF", "Repressed_ZBTB6"="#3C5488FF"), guide=NULL)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig2,"f2i.neuron_tCRE_active_percentage.pdf"), width = 1.4, height = 2)
print(f2i)
dev.off() 

#==============















# figure not included
#===============================================================================
#===============================================================================
# fig 2c-promoter
aCREe=read.delim(paste0(path_fig2_data,"peaks.merged.all.markalone_k16.subclass.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

aCREp1=aCREe[which(aCREe$promoter_type2=="promoter-like"),]
aCREp2=aCREp1[,c(3,5,7)]
colnames(aCREp2)=c("iPSC","NSC","Neuron")
aCREp2[is.na(aCREp2)]=0
aCREp2[aCREp2=="CTCF"]=0
aCREp2$iPSC[grep("Promoter",aCREp2$iPSC)]=3
aCREp2$NSC[grep("Promoter",aCREp2$NSC)]=3
aCREp2$Neuron[grep("Promoter",aCREp2$Neuron)]=3
aCREp2$iPSC[grep("Flanking_promoter",aCREp2$iPSC)]=3
aCREp2$NSC[grep("Flanking_promoter",aCREp2$NSC)]=3
aCREp2$Neuron[grep("Flanking_promoter",aCREp2$Neuron)]=3 #red
aCREp2$iPSC[grep("Bivalent_promoter",aCREp2$iPSC)]=2
aCREp2$NSC[grep("Bivalent_promoter",aCREp2$NSC)]=2
aCREp2$Neuron[grep("Bivalent_promoter",aCREp2$Neuron)]=2 #green
aCREp2$iPSC[grep("enhancer",aCREp2$iPSC)]=1
aCREp2$NSC[grep("enhancer",aCREp2$NSC)]=1
aCREp2$Neuron[grep("enhancer",aCREp2$Neuron)]=1 #brown

aCREp2 <- aCREp2 %>% mutate_if(is.character, as.numeric)
aCREp3=aCREp2%>%group_by(iPSC,NSC,Neuron)%>%dplyr::summarise(count=n())
aCREp3$count=round(aCREp3$count/20)
aCREp4 <- aCREp3[rep(row.names(aCREp3), aCREp3$count), c(1:3)]

out4=Heatmap(as.matrix(aCREp4),
             column_title = "34,435 promoter-\nlike CREs", column_title_gp = gpar(fontsize = 7),
             column_names_gp = gpar(fontsize = 5), column_names_rot = 25, 
             col = colorRampPalette(c("white","#00A087FF","#4DBBD5FF", "#E64B35FF"))(4),
             use_raster = T,
             cluster_columns = F,
             show_heatmap_legend = T,
             row_dend_reorder = FALSE,
             row_dend_width = unit(5, "mm"), 
             heatmap_legend_param = list(
               title = NULL,
               labels_gp = gpar(fontsize = 5),
               at = c("1", "2", "3"),
               labels = c("Repressive", "Primed", "Active"),
               grid_height = unit(2, "mm"),
               grid_width = unit(2, "mm"),
               direction = "vertical"))

pdf(paste0(path_fig2,"f2c.heatmap_34435paCRE.pdf"), width=1.1, height=2.2)
draw(out4,heatmap_legend_side = "bottom")
dev.off()

#======================
# repressed enhancer heatmap

mboth_table2=read.delim(paste0(motif_folder,"repressedEnhancer_motif/all_repressed_enhancer_state.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

mboth_table3=mboth_table2%>%group_by(group2,iPSC,NSC,Neuron)%>%dplyr::summarise(count=n())
mboth_table3$count=round(mboth_table3$count/10)
mboth_table4 <- mboth_table3[rep(row.names(mboth_table3), mboth_table3$count), c(1:4)]

out4a=Heatmap(as.matrix(mboth_table4[which(mboth_table4$group2=="With transcription"),c(2:4)]),
              column_title = "Transcribed\nenhancer 1,519",column_title_gp = gpar(fontsize = 7),
              column_names_gp = gpar(fontsize = 5), column_names_rot = 25, 
              col = c("0" = "white", "1" = "beige", "2" = "#4DBBD5FF", "3" = "#00A087FF", "4" = "#E64B35FF"),
              use_raster = T,
              cluster_columns = F,
              show_heatmap_legend = F,
              row_dend_reorder = FALSE,
              row_dend_width = unit(4, "mm"),
              row_dend_gp = gpar(lwd = 0.5))
out4b=Heatmap(as.matrix(mboth_table4[which(mboth_table4$group2=="Without transcription"),c(2:4)]),
              column_title = "Non-transcribed\nenhancer 18,324",column_title_gp = gpar(fontsize = 7),
              column_names_gp = gpar(fontsize = 5), column_names_rot = 25, 
              col = c("0" = "white", "1" = "beige", "2" = "#4DBBD5FF", "3" = "#00A087FF", "4" = "#E64B35FF"),
              use_raster = T,
              cluster_columns = F,
              show_heatmap_legend = F,
              row_dend_reorder = FALSE,
              row_dend_width = unit(4, "mm"),
              row_dend_gp = gpar(lwd = 0.5))

legend_labels = c("Primed","Repressive","Active")
legend_colors = c("2" = "#4DBBD5FF", "3" = "#00A087FF", "4" = "#E64B35FF")
global_legend = Legend(
  labels = legend_labels,
  legend_gp = gpar(fill = legend_colors),
  direction = "vertical",
  ncol = 1,
  labels_gp = gpar(fontsize = 5),
  title_gp = gpar(fontsize = 6))

g1 = grid.grabExpr(draw(out4a))
g2 = grid.grabExpr(draw(out4b))
gl = grid.grabExpr(draw(global_legend))

final <- plot_grid(g1, g2, gl, ncol = 3, rel_widths = c(1, 1,0.5))

grid.draw(final)



