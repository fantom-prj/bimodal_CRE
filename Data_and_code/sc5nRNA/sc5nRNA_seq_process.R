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
library(ggnewscale)
library(Seurat)
set.seed(2022)
library(Signac)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(GenomicRanges)
library(future)


###color#####
library("ggsci")
mypal = pal_npg("nrc", alpha = 1)(9)
mypal
library("scales")
show_col(mypal)
###############

#primary_folder=[primary_folder]
primary_folder="/analysisdata/fantom6/Interactome/single_cell_wallace/Figure/"
RNA=paste0(primary_folder,"Data_and_code/sc5nRNA/")
setwd(RNA)

#the single-cell 5'end RNAseq from 10x was processed by 4 R script in the order of:
#scRNA_pre_processing.Rmd
#scRNA_QC_filtering.Rmd
#scRNA_integration.Rmd
#scRNA_monocle.Rmd

#final scRNAseq object is saved as paste0(primary_folder,"Data_and_code/DATA/20220331_All.rna.flt.rds") 


#===============================================================================
## Trajectory inference using Monocle3
All.rna.flt <- readRDS(paste0(primary_folder,"Data_and_code/DATA/20220331_All.rna.flt.rds"))

# update the cluster name for metacell calling
library(forcats)
All.rna.flt$celltype_manual <- fct_collapse(All.rna.flt$celltype_manual, neuron_progenitor_n_schwann = c("neuron_progenitor", "neuron_schwann_like"))
saveRDS(All.rna.flt, paste0(primary_folder,"Data_and_code/DATA/20220331_All.rna.flt.rds"))

All.rna.cp <- All.rna.flt

All.rna.cp@reductions[["umap"]] <- All.rna.cp@reductions[["umap_harmony"]]
All.rna.cp@reductions[["pca"]] <- All.rna.cp@reductions[["harmony"]]

# Convert to cds
cds <- as.cell_data_set(All.rna.cp)
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)

# learn trajectory
# don't use partition for a fully-connected graph
cds <- learn_graph(cds, use_partition = FALSE)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

# select root cells
ips_undiff_ids <- colnames(subset(All.rna.cp, celltype_manual == "iPS_0_undifferentiaed"))

# find the root principal node
closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
closest_vertex <- closest_vertex[ips_undiff_ids, ]
closest_vertex <- as.numeric(names(which.max(table(closest_vertex))))
mst <- principal_graph(cds)$UMAP
root_pr_nodes <- igraph::V(mst)$name[closest_vertex]

# Calculate pseudotime
cds <- order_cells(cds, root_pr_nodes = root_pr_nodes)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

library(igraph)
nodes <- data.frame(t(cds@principal_graph_aux$UMAP$dp_mst))
nodes$name=rownames(nodes)
colnames(nodes)[c(1:2)]=c("x_coord","y_coord")
edges <- igraph::as_data_frame(cds@principal_graph$UMAP, what = "edges")
edges=left_join(edges, nodes, by=c("from" = "name"), copy=F)
edges=left_join(edges, nodes, by=c("to" = "name"), copy=F)
colnames(edges)[c(4:7)]=c("x_start","y_start","x_end","y_end")
write.table(edges,gzfile(paste0(path_fig1_data,"harmony.scRNA.umap.with.monocle.pseudotime.edge.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
# -> for fig 1d

# DE along pseudotime for Neuron trajectory [trajectory 1]
  
cds1 <- cds[, !(colData(cds)$celltype_manual %in% c("OPC_like","NSC_2_astrocyte_like"))]
plot_cells(cds1, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
DEG1 <- graph_test(cds1, neighbor_graph="principal_graph", cores = 1, alternative = "greater")

DEG1$geneSymbol <- row.names(DEG1)
DEG1 <- arrange(DEG1, desc(morans_I))
setwd(paste0(RNA,"output/"))
write.table(DEG1, file = gzfile("20220406_DEG1_graph_test_neuron.csv.gz"), quote = FALSE, sep = ",", row.names = F, col.names = T)

# DE along pseudotime for OPC trajectory [trajectory 2]

cds2 <- cds[, !(colData(cds)$celltype_manual %in% c("neuron_mature_2","neuron_mature","neuron_immature","neuron_progenitor","neuron_schwann_like"))]
plot_cells(cds2, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
DEG2 <- graph_test(cds2, neighbor_graph="principal_graph", cores = 1, alternative = "greater")

DEG2$geneSymbol <- row.names(DEG2)
DEG2 <- arrange(DEG2, desc(morans_I))
setwd(paste0(RNA,"output/"))
write.table(DEG2, file = gzfile("20220406_DEG2_graph_test_opc.csv.gz"), quote = FALSE, sep = ",", row.names = F, col.names = T)

T_all=data.frame(principal_graph_aux(cds)$UMAP$pseudotime)
T_all2=data.frame(colData(cds))
colnames(T_all)="pseudotime"
T_all$bin=paste0("T",floor(T_all$pseudotime/2)+1)
T_all$cellID=rownames(T_all)
T_all2$cellID=rownames(T_all2)
T_all=full_join(T_all2, T_all, by="cellID", copy=F)
T_all$major_trajectory=1
T_all$major_trajectory[which(T_all$celltype_manual %in% c("OPC_like","NSC_2_astrocyte_like"))]=0
T_all$minor_trajectory=1
T_all$minor_trajectory[which(T_all$celltype_manual %in% c("neuron_mature_2","neuron_mature","neuron_immature","neuron_progenitor","neuron_schwann_like"))]=0

setwd(paste0(RNA,"output/"))
write.table(T_all, gzfile("pseudotime_table_scRNA.tsv.gz"), sep="\t", quote=F, col.names=T, row.names=F)
write.table(T_all,gzfile(paste0(path_fig1_data,"pseudotime_table_scRNA.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
# -> for fig 1d

#===============================================================================
#obtain the cell barcodes for 10 cell clusters, for visualization
setwd(paste0(primary_folder,"Data_and_code/zenbu/subset_TSS"))
all_clusters <- unique(All.rna.flt@meta.data$celltype_manual)
path1=paste0(primary_folder,"Data_and_code/SCAFE/sc5end/solo/")
cluster_info=read.delim(paste0(primary_folder,"Fig1/cluster_info.tsv"),header=T)
cluster_info=cluster_info[c(1:5,10,6:9),]

need=c("iPS1","iPS2","Neuron1","Neuron2","NSC1","NSC2")
TSS1=data.frame()
for (i in 1:length(need)){
  TSS=fread(paste0(path1,need[i],"/bam_to_ctss/",need[i],"/bed/",need[i],".CB.ctss.bed.gz"),head=F)
  TSS$V4=paste0(TSS$V4,"_",i)
  TSS1=rbind(TSS1,TSS)}

for (j in c(1:nrow(cluster_info))) {
  cluster_barcodes <- rownames(All.rna.flt@meta.data)[All.rna.flt@meta.data$celltype_manual == cluster_info$cluster[j]]
  #cluster_barcodes=sapply(strsplit(cluster_barcodes,"_"),"[",1)
  TSS2=TSS1[which(TSS1$V4 %in% cluster_barcodes),]
  TSS3=TSS2%>%group_by(V1,V2,V3,V6)%>%dplyr::summarise(V5=sum(V5))
  TSS3$V4=paste0(TSS3$V1,"_",TSS3$V2,"_",TSS3$V3,"_",TSS3$V6)
  write.table(TSS3[order(TSS3$V1, TSS3$V2),c(1,2,3,6,5,4)],  file = gzfile(paste0("0",j,"_",cluster_info$label[j], ".CTSS.bed.gz")), quote = FALSE, row.names = FALSE, col.names = FALSE)
}

#===============================================================================
#upload to zenbu

#===============================================================================
# gene expression matrix per sample

All.rna.flt <- readRDS(paste0(primary_folder,"Data_and_code/DATA/20220331_All.rna.flt.rds"))
avg_expr <- AggregateExpression(
  All.rna.flt,
  assays = "RNA",
  group.by = "sampleType",
  slot = "data",
  return.seurat = FALSE
)$RNA

avg_expr=data.frame(avg_expr)
write.table(avg_expr, gzfile(paste0(RNA,"output/gene_normalizedExp_perSample.tsv.gz")), col.names=T, row.names=T, sep="\t", quote=F)

#============================
# gene count matrix per sample

All.rna.flt <- readRDS(paste0(primary_folder,"Data_and_code/DATA/20220331_All.rna.flt.rds"))
sum_count <- AggregateExpression(
  All.rna.flt,
  assays = "RNA",
  group.by = "sampleType",
  slot = "count",
  return.seurat = FALSE
)$RNA

sum_count=data.frame(sum_count)
dge <- DGEList(sum_count)
dge <- calcNormFactors(dge, method = "RLE")
CPM <- edgeR::cpm(dge)

write.table(CPM, gzfile(paste0(RNA,"output/gene_RLE_CPM_perSample.tsv.gz")), col.names=T, row.names=T, sep="\t", quote=F)




