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

#======
# ex1a
clusters_of_interest=c("iPSC","iPSC","iPSC","iPSC","NSC","NSC","NSC","NPC_like","Immature_neuron","Immature_neuron","Mature_neuron_1","Mature_neuron_1","Mature_neuron_2","Mature_neuron_2","OPC_like")
genes_of_interest=c("POU5F1","EPCAM","HIST1H1E","UBE2C","CDKN1A","WLS","TAGLN","NEUROD1","NHLH1","NEFM","STMN2","ONECUT2","RTN1","HOXB9","DCN")
interest=data.frame(cbind(clusters_of_interest,genes_of_interest))
interest$header=paste0(interest$clusters_of_interest,":\n",interest$genes_of_interest)
interest1=interest[c(1,6,8,9,12,14,15),] # select representative for UMAP

Aumap3=read.delim(paste0(path_fig1_data,"snATAC_n15048.with_UMAP_meta.tsv.gz"), header=T, stringsAsFactors =F, check.names = F)
Bumap3=read.delim(paste0(path_fig1_data,"snATAC_n15048.with_UMAP_meta_imputed_RNA.tsv.gz"), header=T, stringsAsFactors =F, check.names = F)
umap3=read.delim(paste0(path_fig1_data,"scRNA_n24434.with_UMAP_meta.tsv.gz"), header=T, stringsAsFactors =F, check.names = F)

Aumap3a=reshape2::melt(Aumap3[,c(1,4,5,7,9:23)], id=c(1:4))
Bumap3a=reshape2::melt(Bumap3[,c(1,4,5,7,9:23)], id=c(1:4))
umap3a=reshape2::melt(umap3[,c(1,4,5,7,9:23)], id=c(1:4))

Aumap3a$group="snATAC_gene_activity"
Bumap3a$group="snATAC_imputed_RNA"
umap3a$group="sc5nRNA_expression"
umap3a=rbind(Aumap3a,Bumap3a,umap3a)

umap3c=umap3a[which(umap3a$variable %in% interest1$genes_of_interest),]
umap3c=left_join(umap3c,interest1, by=c("variable"="genes_of_interest"),copy=F)
umap3c$header=factor(umap3c$header,levels=interest1$header)
umap3c=umap3c%>%group_by(group,header)%>%dplyr::mutate(scaled_value=value/max(value))
umap3c=umap3c[sample(nrow(umap3c),nrow(umap3c)),]
ex1a=ggplot() + 
  labs(x="UMAP1", y ="UMAP2", title=NULL, color="Scaled\nexpression")+
  facet_grid(cols=vars(header), rows=vars(group))+
  geom_point_rast(data=umap3c[which(umap3c$group %in% c("sc5nRNA_expression")),], mapping=aes(x=UMAP1, color=scaled_value, y=UMAP2), shape=16, alpha=1, size=0.1)+
  scale_color_gradient(low="grey90",high="#3C5488FF")+
  theme1+theme(panel.grid.major = element_blank(),legend.position = "right")
pdf(paste0(path_fig1,"ex1a.RNA.marker.gene.express.umap.pdf"), width = 7, height = 1.2)
print (ex1a)
dev.off() 
#
ex1a2=ggplot() + 
  labs(x="UMAP1", y ="UMAP2", title=NULL, color="Scaled\nexpression")+
  facet_grid(cols=vars(header), rows=vars(group))+
  geom_point_rast(data=umap3c[which(umap3c$group %in% c("snATAC_imputed_RNA")),], mapping=aes(x=UMAP1, color=scaled_value, y=UMAP2), shape=16, alpha=1, size=0.1)+
  scale_color_gradient(low="grey90",high="darkgreen")+
  theme1+theme(panel.grid.major = element_blank(),legend.position = "right")
pdf(paste0(path_fig1,"ex1a.ATAC.marker.gene.impute.umap.pdf"), width = 7, height = 1.2)
print (ex1a2)
dev.off() 

#======
# ex1b-c
# co-embed with scDART
pca1=read.delim(paste0(path_fig1_data,"scDART_SEACell_for_all_clusters.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
# generated from Fig1.coembed/conembed.R

nrow(unique(pca1[,c(13,14)])) #391 #removed metacell with <3 ATAC or <3 RNA single-cells
pca1=left_join(pca1, cluster_info, by="cluster",copy=F)
pca1$SEACell=sapply(strsplit(pca1$SEACell,"_"), tail, 1)
pca1$SEACell2=paste0(pca1$label,"_",pca1$SEACell)
pca2=pca1%>%group_by(SEACell2,label,dataType)%>%dplyr::summarise(count=n())
pca3=spread(pca2,key=3, value=4)
pca3[is.na(pca3)]=0
pca1$label=factor(pca1$label, levels=label)
pca1$dataType=paste0(pca1$dataType,"-seq")

ex1b=ggplot() + 
  labs(x="scDART_UMAP1", y ="scDART_UMAP2", title="scDART co-embedding (colored by data type)", color="Data Type")+
  geom_point_rast(data=pca1[sample(nrow(pca1),nrow(pca1)),], mapping=aes(x=UMAP1, color=as.factor(dataType), y=UMAP2), shape=16, size=0.15, alpha=0.6)+
  scale_color_npg(labels=c("snATAC-seq","sc5nRNA-seq"))+
  facet_wrap(facets=vars(label), ncol=4,nrow=3, scale="free")+
  theme1 + theme(panel.grid.major = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), legend.position = c(0.7,0.15))
pdf(paste0(path_fig1,"ex1b.coembed.scDART.colored.by.tech.UMAP.pdf"), width = 3.6, height = 3)
print (ex1b)
dev.off() 

ex1c=ggplot() + 
  labs(x="scDART_UMAP1", y ="scDART_UMAP2", title="scDART co-embedding (colored by SEACell)")+
  geom_point_rast(data=pca1[sample(nrow(pca1),nrow(pca1)),], mapping=aes(x=UMAP1, color=as.factor(SEACell), y=UMAP2), shape=16, size=0.15, alpha=0.6)+
  #geom_mark_ellipse(data=pca1, aes(x=X0, color=as.factor(SEACell), y=X1), expand = unit(1, "mm") )+
  facet_wrap(facets=vars(label), ncol=4,nrow=3, scale="free")+
  theme1 + theme(panel.grid.major = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), legend.position = "none")
pdf(paste0(path_fig1,"ex1c.coembed.scDART.colored.by.SEACell.UMAP.pdf"), width = 3.6, height = 3)
print (ex1c)
dev.off() 

#==================


