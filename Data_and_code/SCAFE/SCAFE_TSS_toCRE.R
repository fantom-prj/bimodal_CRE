library(Seurat)

#####################
#primary_folder=[primary_folder]
primary_folder="/analysisdata/fantom6/Interactome/single_cell_wallace/Figure/"
SCAFE_path=paste0(primary_folder,"Data_and_code/SCAFE/sc5end/")
scRNA_path=paste0(primary_folder,"Data_and_code/sc5nRNA/")

#===============================================================================
#count valid CTSS from SCAFE into scATAC peak
#use count from SCAFE
setwd(paste0(SCAFE_path,"count/"))
system("perl 20230711.batch_run_scafe_interactome,scafe_count.pl")

# get SCAFE count from aCRE_original
seu.list <- list()
for(i in c("iPS1","iPS2","Neuron1","Neuron2","NSC1","NSC2")){
  tmp <- Read10X(data.dir = paste0(SCAFE_path,"count/aCRE_original/",i,"/matrix/"))
  tmp <- CreateSeuratObject(counts = tmp, project = i, min.cells = 0, min.features = 0)
  seu.list[[i]] <- tmp}

tmp<-c()
for(i in 2:length(seu.list)){
  tmp <- c(tmp, seu.list[[i]])}
seu.all <- merge(seu.list[[1]], y = tmp) # "_1" ~ "_6" added to the end of original cellbarcodes
rm(tmp)

seu.all$sampleID <- factor(seu.all$orig.ident, c("iPS1","iPS2","NSC1","NSC2","Neuron1","Neuron2"))

scRNA=readRDS(paste0(scRNA_path,"20220331_All.rna.flt.rds"))
seu.all_filtered <- seu.all[,colnames(seu.all) %in% colnames(scRNA)] #limited to cells used in scRNA-seq
#normalization
hvg <- VariableFeatures(seu.all_filtered)
res <- seq(0.5, 3, 0.5)
ndims <- 30
vars.reg <- c('nCount_RNA')
seu.all_filtered <- NormalizeData(seu.all_filtered)
seu.all_filtered <- FindVariableFeatures(seu.all_filtered, selection.method = 'vst', nfeatures = 2000)
seu.all_filtered %<>%
  ScaleData(vars.to.regress = vars.reg) %>%
  RunPCA(features = hvg) %>%
  FindNeighbors(dims = 1:ndims) %>%
  FindClusters(resolution = res) %>%
  RunUMAP(dims = 1:ndims) %>%
  RunTSNE(dims = 1:ndims)
FeaturePlot(seu.all_filtered, features="nCount_RNA", pt.size = 0.01)
saveRDS(seu.all_filtered, "/analysisdata/fantom6/Interactome/single_cell_wallace/SCAFE/scRNA_aCRE_orig.rds")

#===============================================================================





