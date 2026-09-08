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
library(Seurat)

###color#####
library("ggsci")
mypal = pal_npg("nrc", alpha = 1)(10)
mypal2=pal_bmj()(9)
mypal
library("scales")
show_col(mypal)
show_col(mypal2)

#####################
primary_folder=[primary_folder]
#primary_folder="/analysisdata/fantom6/Interactome/single_cell_wallace/Figure/"
path_fig1_data=paste0(primary_folder,"Fig1/data/")
path_fig1=paste0(primary_folder,"Fig1/out/")

setwd(path_fig1)
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

#============================================
cluster=c("iPS_0_undifferentiaed","iPS_1_neuroepithelial","NSC_0_iPS_like", "NSC_1_differentiating", "NSC_2_astrocyte_like", "neuron_progenitor_n_schwann", "neuron_immature", "neuron_mature","neuron_mature_2", "OPC_like")
cluster1=c("iPSC_0","iPSC_1","NSC_0", "NSC_1", "NSC_2", "Neuron_0", "Neuron_1", "Neuron_2","Neuron_3", "Neuron_4")
label=c("iPSC_Sphase","iPSC_G2M","NSC_stem","NSC_1","NSC_2","NPC_like","Immature_neuron","Mature_neuron_1","Mature_neuron_2","OPC_like")
cluster_info=data.frame(cbind(cluster, cluster1, label))

#=======
# fig 1b
umap3=read.delim(paste0(path_fig1_data,"harmony.scRNA.umap.with.metacell.tsv.gz"), header=T, stringsAsFactors =F, check.names = F)
umap3=left_join(umap3,cluster_info, by="cluster", copy=F)
umap3$label=factor(umap3$label, levels=label)
umap3$tech="UMAP of sc5nRNA-seq"
colnames(umap3)[c(4,5)]=c("UMAP1", "UMAP2")
Aumap3=read.delim(paste0(path_fig1_data,"harmony.snATAC.umap.with.metacell.tsv.gz"), header=T, stringsAsFactors =F, check.names = F)
Aumap3=left_join(Aumap3,cluster_info, by="cluster", copy=F)
Aumap3$label=factor(Aumap3$label, levels=label)
Aumap3$tech="UMAP of snATAC-seq"
colnames(Aumap3)[c(4,5)]=c("UMAP1", "UMAP2")

umap4=umap3[which(umap3$group == "single cell"),]
Aumap4=Aumap3[which(Aumap3$group == "single cell"),]
umap5=rbind(umap4,Aumap4)
umap_info=umap5%>%group_by(tech,label)%>%dplyr::summarise(UMAP1=mean(UMAP1),UMAP2=mean(UMAP2))

f1b=ggplot() + 
  labs(x="UMAP1", y ="UMAP2", title=NULL, color="Cluster")+
  facet_wrap(vars(tech), ncol=2, scales="free")+
  scale_color_npg()+
  geom_point_rast(data=umap5, mapping=aes(x=UMAP1, color=label, y=UMAP2), shape=16, alpha=0.6, size=0.15)+
  geom_text_repel(data=umap_info, mapping=aes(x=UMAP1, y=UMAP2, label=label), segment.size=0.2, size=1.8, color="black", bg.color = "white", bg.r=0.1)+
  theme1+theme(legend.position = "none",  panel.grid.major = element_blank())
pdf(paste0(path_fig1,"f1b.harmony.scRNA.snATAC.umap.pdf"), width = 3.4, height = 1.7)
print (f1b)
dev.off() 


#=====================================
# fig 1c
cluster=c("iPS_0_undifferentiaed","iPS_1_neuroepithelial","NSC_0_iPS_like", "NSC_1_differentiating", "NSC_2_astrocyte_like", "neuron_progenitor_n_schwann", "neuron_immature", "neuron_mature","neuron_mature_2", "OPC_like")
cluster1=c("iPSC_0","iPSC_1","NSC_0", "NSC_1", "NSC_2", "Neuron_0", "Neuron_1", "Neuron_2","Neuron_3", "Neuron_4")
label=c("iPSC_Sphase","iPSC_G2M","NSC_stem","NSC_1","NSC_2","NPC_like","Immature_neuron","Mature_neuron_1","Mature_neuron_2","OPC_like")
cluster_info=data.frame(cbind(cluster, cluster1, label))

meta=read.delim(paste0(path_fig1_data,"metadata.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
meta$celltype="iPSC"
meta$celltype[grep("NSC",meta$V1)]="NSC"
meta$celltype[grep("NRN",meta$V1)]="Neuron"
meta$celltype[grep("_5",meta$V1)]="NSC"
meta$celltype[grep("_6",meta$V1)]="NSC"
meta$celltype[grep("_3",meta$V1)]="Neuron"
meta$celltype[grep("_4",meta$V1)]="Neuron"

meta$dataType[which(meta$dataType == "ATAC")]="snATAC-seq"
meta$dataType[which(meta$dataType == "RNA")]="sc5nRNA-seq"
meta=left_join(meta, cluster_info, by="cluster", copy=F)
meta1=meta%>%group_by(dataType,celltype,label)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
meta1$label=factor(meta1$label, levels=label)
meta1$celltype=factor(meta1$celltype, levels=c("iPSC","NSC","Neuron"))

f1c=ggplot() + 
  labs(x=NULL, y = "% of cells",title= "Cell clusters from samples", fill=NULL) +
  scale_fill_npg()+
  scale_y_continuous(labels = scales::percent, limits=c(0,1), breaks=c(0,0.25,0.5,0.75,1))+
  facet_grid(cols=vars(dataType))+
  geom_bar(data=meta1, mapping=aes(y=percent, x=celltype, fill=label), linewidth=0.25, color = "black", stat="identity", alpha=1) + 
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1), 
               legend.box.margin = margin(t=-10, r=0, b=0, l=-10),legend.position = "right", legend.direction = "vertical")
pdf(paste0(path_fig1,"f1c.cluster_per_sample.pdf"), width = 2.3, height = 1.7)
print(f1c)
dev.off()

#============
# fig 1d
umap3=read.delim(paste0(path_fig1_data,"harmony.scRNA.umap.with.metacell.tsv.gz"), header=T, stringsAsFactors =F, check.names = F)
umap3=left_join(umap3,cluster_info, by="cluster", copy=F)
umap3$label=factor(umap3$label, levels=label)
umap3$tech="UMAP of sc5nRNA-seq"
colnames(umap3)[c(4,5)]=c("UMAP1", "UMAP2")
umap4=umap3[which(umap3$group == "single cell"),]

edges=read.delim(paste0(path_fig1_data,"harmony.scRNA.umap.with.monocle.pseudotime.edge.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
pseudotime=read.delim(paste0(path_fig1_data,"pseudotime_table_scRNA.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
umap4=left_join(umap4,pseudotime[,c(37,38)], by="cellID", copy=F)

f1d=ggplot() + 
  labs(x="UMAP1", y ="UMAP2", title="Pseudotime from sc5nRNA-seq", color=NULL)+
  geom_point_rast(data=umap4, mapping=aes(x=UMAP1, color=pseudotime, y=UMAP2), shape=16, alpha=0.6, size=0.15)+
  scale_color_viridis_c(option = "plasma", direction = (-1))+
  geom_segment(data = edges, aes(x = x_start, y = y_start, xend = x_end, yend = y_end),lineend="round", color = "grey20", linewidth = 0.5) +
  #geom_text_repel(data=umap_info[which(umap_info$tech=="UMAP of sc5nRNA-seq"),], mapping=aes(x=UMAP1, y=UMAP2, label=label), size=1.8, color="black", bg.color = "white", bg.r=0.1)+
  theme1+theme(legend.position = c(0.2,0.2), panel.grid.major = element_blank())
pdf(paste0(path_fig1,"f1d.harmony.scRNA.umap.with.monocle.pseudotime.pdf"), width = 1.7, height = 1.7)
print (f1d)
dev.off()

#========
# f1e heatmap for sc5 and ATAC imput RNA
clusters_of_interest=c("iPSC","iPSC","iPSC","iPSC","NSC","NSC","NSC","NPC_like","Immature_neuron","Immature_neuron","Mature_neuron_1","Mature_neuron_1","Mature_neuron_2","Mature_neuron_2","OPC_like")
genes_of_interest=c("POU5F1","EPCAM","HIST1H1E","UBE2C","CDKN1A","WLS","TAGLN","NEUROD1","NHLH1","NEFM","STMN2","ONECUT2","RTN1","HOXB9","DCN")
interest=data.frame(cbind(clusters_of_interest,genes_of_interest))
interest$header=paste0(interest$clusters_of_interest,":\n",interest$genes_of_interest)

Aumap3=read.delim(paste0(path_fig1_data,"snATAC_n15048.with_UMAP_meta.tsv.gz"), header=T, stringsAsFactors =F, check.names = F)
Bumap3=read.delim(paste0(path_fig1_data,"snATAC_n15048.with_UMAP_meta_imputed_RNA.tsv.gz"), header=T, stringsAsFactors =F, check.names = F)
umap3=read.delim(paste0(path_fig1_data,"scRNA_n24434.with_UMAP_meta.tsv.gz"), header=T, stringsAsFactors =F, check.names = F)
# generated from Fig1.metacell/meta.cell.object.R

Aumap3a=reshape2::melt(Aumap3[,c(1,4,5,7,9:23)], id=c(1:4))
Bumap3a=reshape2::melt(Bumap3[,c(1,4,5,7,9:23)], id=c(1:4))
umap3a=reshape2::melt(umap3[,c(1,4,5,7,9:23)], id=c(1:4))

Aumap3a$group="snATAC_gene_activity"
Bumap3a$group="snATAC_imputed_RNA"
umap3a$group="sc5nRNA_expression"
umap3a=rbind(Aumap3a,Bumap3a,umap3a)
umap3b=umap3a%>%group_by(group,label,variable)%>%dplyr::summarise(value=mean(value))
umap3b=umap3b%>%group_by(group,variable)%>%dplyr::mutate(scaled_value=value/max(value))
umap3b$variable=factor(umap3b$variable, levels=genes_of_interest)
umap3b$label=factor(umap3b$label, levels=label)

f1e=ggplot() + 
  labs(x=NULL, y =NULL, title=NULL, fill="Scaled\nexpression")+
  facet_grid(cols=vars(group))+
  scale_y_discrete(position = "right", limits=rev)+
  geom_tile(data=umap3b[which(umap3b$group %in% c("sc5nRNA_expression","snATAC_imputed_RNA")),], mapping=aes(x=label, fill=scaled_value, y=variable), linewidth=0.01)+
  scale_fill_gradient2(low="#3C5488FF",mid="white",high="#F0AB00FF", midpoint=0.6)+
  theme1+theme(panel.grid.major = element_blank(),legend.position = "right",axis.text.x = element_text(color ="black", angle=90, hjust=1, vjust=0.5))
pdf(paste0(path_fig1,"f1e.heatmap_cluster.marker.gene.pdf"), width = 2.8, height = 2.2)
print (f1e)
dev.off() 

#===============
# fig 1f
pca1=read.delim(paste0(path_fig1_data,"scDART_SEACell_for_all_clusters.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
# generated from Fig1.coembed/conembed.R

nrow(unique(pca1[,c(13,14)])) #395
pca1=left_join(pca1, cluster_info, by="cluster",copy=F)
pca1$SEACell=sapply(strsplit(pca1$SEACell,"_"), tail, 1)
pca1$SEACell2=paste0(pca1$label,"_",pca1$SEACell)
pca2=pca1%>%group_by(SEACell2,label,dataType)%>%dplyr::summarise(count=n())
pca3=spread(pca2,key=3, value=4)
pca3[is.na(pca3)]=0
exclude=unique(pca3$SEACell2[which(pca3$ATAC<3 | pca3$RNA<3 )])
pca1$label=factor(pca1$label, levels=label)
pca1$dataType=paste0(pca1$dataType,"-seq")

pca3$include="Yes"
pca3$include[which(pca3$SEACell2 %in% exclude)]="No"
pca3$RNA=pca3$RNA*(-1)
pca3$both=pca3$ATAC-pca3$RNA
pca4=melt(pca3, id=c(1,2,5,6))
pca4$label=factor(pca4$label, levels=label)

f1f=ggplot() + 
  labs(x="Number of single cell (24,342 |15,046)", y ="Metacell (n = 391)", title="Content of metacells", fill=NULL)+
  geom_bar(data=pca4[which(pca4$include == "Yes"),], mapping=aes(y=reorder(SEACell2,-both), x=value, fill=as.factor(variable)), colour=NA,size=0, width = 1, stat="identity")+
  scale_fill_manual(labels=c("scATAC","sc5nRNA"), values=c("#00A087FF","#2A6EBBFF"))+
  facet_grid(rows=vars(label), scale="free_y", space="free")+
  scale_x_continuous(breaks=c(-500, -250, 0, 250), labels=c(500,250,0,250))+
  theme1 + theme(axis.text.y = element_blank(), strip.text.y = element_text(size=5, angle=0, margin = margin(0,0,0,0, "cm")), 
                 legend.position = c(0.1,0.9), axis.line.y = element_blank(), axis.ticks.y = element_blank(),panel.grid.major = element_blank(),
                 panel.spacing = unit(0.1, "lines"))
pdf(paste0(path_fig1,"f1f.metacell.content.pdf"), width = 2.2, height = 2.2)
print (f1f)
dev.off()

#=======
# fig 1g
umap3=read.delim(paste0(path_fig1_data,"harmony.scRNA.umap.with.metacell.tsv.gz"), header=T, stringsAsFactors =F, check.names = F)
umap3=left_join(umap3,cluster_info, by="cluster", copy=F)
umap3$label=factor(umap3$label, levels=label)
umap3$tech="UMAP of sc5nRNA-seq"
colnames(umap3)[c(4,5)]=c("UMAP1", "UMAP2")

f1g=ggplot() + 
  labs(x="UMAP1", y ="UMAP2", title="Metacells overlay sc5nRNA-seq UMAP", color="Cluster")+
  geom_point_rast(data=umap3, mapping=aes(x=UMAP1, alpha=as.factor(group), size=as.factor(group), color=label, y=UMAP2),shape=16)+
  scale_color_npg()+
  scale_alpha_manual(values=c("meta cell"=0.95, "single cell"=0.05), guide="none")+
  scale_size_manual(values=c("meta cell"=0.25, "single cell"=0.1), guide="none")+
  geom_text_repel(data=umap_info[which(umap_info$tech=="UMAP of sc5nRNA-seq"),], mapping=aes(x=UMAP1, y=UMAP2, label=label), segment.size=0.2, size=1.8, color="black", bg.color = "white", bg.r=0.1)+
  theme1+theme(legend.position = "none", panel.grid.major = element_blank())
pdf(paste0(path_fig1,"f1g.harmony.scRNA.umap.with.metacell.pdf"), width = 2.1, height = 2.2)
print (f1g)
dev.off()

#=======================


#figure not included
#===============================================================================
Aumap3=read.delim(paste0(path_fig1_data,"harmony.snATAC.umap.with.metacell.tsv.gz"), header=T, stringsAsFactors =F, check.names = F)
Aumap3=left_join(Aumap3,cluster_info, by="cluster", copy=F)
Aumap3$label=factor(Aumap3$label, levels=label)
Aumap3$tech="UMAP of snATAC-seq"
colnames(Aumap3)[c(4,5)]=c("UMAP1", "UMAP2")

c3=ggplot() + 
  labs(x="UMAP1", y ="UMAP2", title="Metacells overlay snATAC-seq UMAP", color="Cluster")+
  geom_point_rast(data=Aumap3, mapping=aes(x=UMAP1, alpha=as.factor(group), size=as.factor(group), color=label, y=UMAP2),shape=16)+
  scale_color_npg()+
  scale_alpha_manual(values=c("meta cell"=0.95, "single cell"=0.05), guide="none")+
  scale_size_manual(values=c("meta cell"=0.25, "single cell"=0.1), guide="none")+
  geom_text_repel(data=umap_info[which(umap_info$tech=="UMAP of snATAC-seq"),], mapping=aes(x=UMAP1, y=UMAP2, label=label), segment.size=0.2, size=1.8, color="black", bg.color = "white", bg.r=0.1)+
  theme1+theme(legend.position = "none", panel.grid.major = element_blank())
pdf("harmony.snATAC.umap.with.metacell.pdf", width = 2.1, height = 2.2)
print(c3)
dev.off()

