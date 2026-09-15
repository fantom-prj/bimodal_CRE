library(dplyr)
library(magrittr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(patchwork)
library(readxl)
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

#####################

#
#primary_folder=[primary_folder]
primary_folder="/analysisdata/fantom6/Interactome/single_cell_wallace/Figure/"
setwd(paste0(primary_folder,"Data_and_code/scATAC/"))

#===============================================================================
plan("multiprocess", workers = 6)
options(future.globals.maxSize = 50 * 1024^3)

IDList <- c("iPS","NSC","NRN")

# =======================================================
# Need a combined peak set. Merge peak.
# =======================================================
## peaks were merged from the 3 libraries

peak.List <- list()
for(i in IDList){
  tmp <- read.table(
    file = paste0("../F6_Interactome_20220128/scATAC/",i,"/outs/peaks.bed"),
    col.names = c("chr", "start", "end")
  )
  peak.List[[i]] <- tmp
}

## convert to genomic ranges
gr.List <- list()
for(i in IDList){
  tmp <- makeGRangesFromDataFrame(peak.List[[i]])
  gr.List[[i]] <- tmp
}

## Create a unified set of peaks to quantify in each dataset
tmp<- GRanges()
for(i in IDList){tmp <- c(tmp, gr.List[[i]])}
combined.peaks <- reduce(x = tmp)

## Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks


# =======================================================
# Create Fragment objects & Count fragments in peaks
# =======================================================

peak_a <- getPeaks(paste0(primary_folder,"Data_and_code/scATAC/merged_peak/peaks.merged.all.bed"), sort_peaks = TRUE)

md.List <- list()
for(i in IDList){
  tmp <- read.table(
    file = paste0("/analysisdata/fantom6/Interactome/scATACcellranger_Kouno/",i,"/outs/singlecell.csv"), 
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  md.List[[i]] <- tmp
}

## perform an initial filtering of low count cells
## passed_filters => number of non-duplicate, usable read-pairs i.e. "fragments"
for(i in IDList){
  md.List[[i]] <- md.List[[i]][md.List[[i]]$passed_filters > 500, ]
}

## create fragment objects
frags.List <- list()
for(i in IDList){
  tmp <- CreateFragmentObject(
    path = paste0("/analysisdata/fantom6/Interactome/scATACcellranger_Kouno/",i,"/outs/fragments.tsv.gz"),
    cells = rownames(md.List[[i]])
  )
  frags.List[[i]] <- tmp
}

## Count fragments in peaks (takes time...) x core ~30 min for 18 samples, 6 cores ~ 1h
counts.List <- list()
for(i in IDList){
  tmp <-  FeatureMatrix(
    fragments = frags.List[[i]],
    features = peak_a,
    cells = rownames(md.List[[i]])
  )
  counts.List[[i]] <- tmp
}

saveRDS(md.List, file = "/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/md.List.3samples.rds")
saveRDS(frags.List, file = "/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/frags.List.3samples.rds")
saveRDS(counts.List, file = "/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/counts.List.3samples.rds")


# =======================================================
# Create and merge object
# =======================================================

chrom_assay.List <- list()
All.atac.List <- list()
for(i in IDList){
  tmp <-  CreateChromatinAssay(counts = counts.List[[i]], 
                               sep = c(":", "-"), 
                               #genome = 'hg38', # Have to disable this in order to work
                               # https://github.com/timoast/signac/issues/687
                               min.cells = 10,  
                               min.features = 200, 
                               fragments = frags.List[[i]])
  chrom_assay.List[[i]] <- tmp
  
  tmp <-  CreateSeuratObject(counts = chrom_assay.List[[i]], assay = "peaks", meta.data = md.List[[i]][,1:18])
  All.atac.List[[i]] <- tmp
  All.atac.List[[i]]$sampleID <- i
}


## merge all datasets, adding a cell ID to make sure cell names are unique
All.atac.raw <- merge(
  x = All.atac.List[[1]],
  y = All.atac.List[2:(length(All.atac.List))],
  add.cell.ids = IDList
)


# =======================================================
# Add annotation
# =======================================================

## extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

## change to UCSC style since the data was mapped to hg38
annotations <- renameSeqlevels(annotations, mapSeqlevels(seqlevels(annotations), "UCSC")) # https://github.com/timoast/signac/issues/826
#seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

## add the gene information to the object
Annotation(All.atac.raw) <- annotations


# =======================================================
# Compute QC metrics
# =======================================================
# compute nucleosome signal score per cell (3 samples using 6 cores takes 2 min...)
All.atac.raw <- NucleosomeSignal(object = All.atac.raw)

# compute TSS enrichment score per cell (3 samples using 6 cores takes 7 min...)
All.atac.raw <- TSSEnrichment(object = All.atac.raw, fast = FALSE) 

# Compute fraction of reads in peaks (FRiP)
All.atac.raw <- FRiP(object = All.atac.raw, assay = 'peaks', total.fragments = "passed_filters")

# add blacklist ratio 
All.atac.raw$blacklist_ratio <- All.atac.raw$blacklist_region_fragments / All.atac.raw$nCount_peaks * 100

saveRDS(All.atac.raw, file = "/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/All.atac.raw.3samples.rds")
rm(frags.List, gr.List, md.List, peak.List, chrom_assay.List, Marrow.atac.List, counts.List)

All.atac.raw <- readRDS("/analysisdata/fantom6/Interactome/single_cell_wallace/ATAC/All.atac.raw.3samples.rds")


# =======================================================
# QC & Filtering
# =======================================================
# Plot TSS enrichment & nucleosome signal

All.atac.raw$high.tss <- ifelse(All.atac.raw$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(All.atac.raw, group.by = 'high.tss') + NoLegend()

All.atac.raw$nucleosome_group <- ifelse(All.atac.raw$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = All.atac.raw, group.by = 'nucleosome_group', region = 'chr1-1-10000000') 

## plot various metrics
VlnPlot(
  object = All.atac.raw,
  features = c('FRiP', 'nCount_peaks', 'TSS.enrichment', 'nucleosome_signal',"blacklist_ratio"),
  pt.size = 0.01,
  ncol = 5,
  group.by = "is__cell_barcode"
)

All.atac.raw@meta.data %>%
  ggplot(., aes(x=nCount_peaks, y=FRiP, color=is__cell_barcode))+
  geom_point()

All.atac.raw@meta.data %>%
  ggplot(., aes(x=log2(nCount_peaks), y=FRiP, color=is__cell_barcode))+
  geom_point()


## Filter cells selected by cellranger and plot again

All.atac.tmp <- subset(
  x = All.atac.raw,
  subset = is__cell_barcode == 1
)

VlnPlot(
  object = All.atac.tmp,
  features = c('FRiP', 'nCount_peaks', 'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0,
  ncol = 4,
  group.by = "sampleID"
)

All.atac.tmp@meta.data %>%
  ggplot(., aes(x=nCount_peaks, y=FRiP, color=sampleID))+
  geom_point()

All.atac.tmp@meta.data %>%
  ggplot(., aes(x=log2(nCount_peaks), y=FRiP, color=sampleID))+
  geom_point()


## Percentile of each metric

print("===== nCount_peaks =====")
All.atac.raw$nCount_peaks %>% quantile(., c(.01, .05, .25, .5, .75, .95, .99))
dplyr::filter(All.atac.raw@meta.data, is__cell_barcode == 1) %>% .$nCount_peaks %>% quantile(., c(.01, .05, .25, .5, .75, .95, .99))

print("===== FRiP =====")
All.atac.raw$FRiP %>% quantile(., c(.01, .05, .25, .5, .75, .95, .99))
dplyr::filter(All.atac.raw@meta.data, is__cell_barcode == 1) %>% .$FRiP %>% quantile(., c(.01, .05, .25, .5, .75, .95, .99))

print("===== TSS.enrichment =====")
All.atac.raw$TSS.enrichment %>% quantile(., c(.01, .05, .25, .5, .75, .95, .99))
dplyr::filter(All.atac.raw@meta.data, is__cell_barcode == 1) %>% .$TSS.enrichment %>% quantile(., c(.01, .05, .25, .5, .75, .95, .99))

print("===== nucleosome_signal =====")
All.atac.raw$nucleosome_signal %>% quantile(., c(.01, .05, .25, .5, .75, .95, .99))
dplyr::filter(All.atac.raw@meta.data, is__cell_barcode == 1) %>% .$nucleosome_signal %>% quantile(., c(.01, .05, .25, .5, .75, .95, .99))


# =======================================================
# Filter for downstream analysis
# =======================================================

All.atac <- subset(
  x = All.atac.raw,
  subset = nCount_peaks > 1000 &
    nCount_peaks < 20000 &
    FRiP > 0.5 &
    is__cell_barcode == 1 & 
    nucleosome_signal < 4 &
    TSS.enrichment > 3
)

## plot metrics for filtered cells

All.atac$sampleID <- factor(All.atac$sampleID, c("iPS","NSC","NRN"))
All.atac@meta.data %>%
  group_by(sampleID) %>%
  mutate(., nCell=n()) %>%
  ggplot(., aes(x=sampleID, y=nCell, fill = sampleID)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_bw()+
  scale_fill_manual(values=sampleType_color)+
  ggtitle(paste0("Total number of cells after filtering: ", nrow(All.atac@meta.data)))+
  xlab("Sample ID")+ylab("N")

VlnPlot(
  object = All.atac,
  features = c('FRiP', 'nCount_peaks', 'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0,
  ncol = 4,
  group.by = "sampleID"
)

All.atac@meta.data %>%
  ggplot(., aes(x=nCount_peaks, y=FRiP, color=sampleID))+
  geom_point()

All.atac@meta.data %>%
  ggplot(., aes(x=log2(nCount_peaks), y=FRiP, color=sampleID))+
  geom_point()


## Normalization, Select variable features, and linear dimensional reduction
All.atac <- RunTFIDF(All.atac, verbose = TRUE)

#===============================================================================
## using top 25% peaks 
All.atac <- FindTopFeatures(All.atac, min.cutoff = 'q75', verbose = TRUE) # use top 25% variable peaks 
All.atac <- RunSVD(All.atac, verbose = TRUE)

# Check: The first LSI component often captures sequencing depth (technical variation) rather than biological variation, so should be removed from downstream analysis
DepthCor(All.atac)

All.atac <- RunUMAP(object = All.atac, reduction = 'lsi', dims = 2:30) # Not using first LSI component
DimPlot(object = All.atac, group.by = "sampleID") + NoLegend()


#===============================================================================
## using all peaks

All.atac <- FindTopFeatures(All.atac, min.cutoff = 'q0', verbose = TRUE) # use all peaks 
All.atac <- RunSVD(All.atac, verbose = TRUE)

# Check: The first LSI component often captures sequencing depth (technical variation) rather than biological variation, so should be removed from downstream analysis
DepthCor(All.atac)

All.atac <- RunUMAP(object = All.atac, reduction = 'lsi', dims = 2:30) # Not using first LSI component
DimPlot(object = All.atac, group.by = "sampleID") 

#clustering

All.atac <- FindNeighbors(object = All.atac, reduction = 'lsi', dims = 2:30)
All.atac <- FindClusters(object = All.atac, algorithm = 3, resolution = 0.8) # resolution default = 0.8
DimPlot(object = All.atac, label = TRUE) + NoLegend()

FeaturePlot(object = All.atac, features = 'nCount_peaks')
FeaturePlot(object = All.atac, features = 'total')
FeaturePlot(object = All.atac, features = 'TSS.enrichment')
FeaturePlot(object = All.atac, features = 'FRiP')

DimPlot(All.atac, group.by = 'sampleID')


# =======================================================
# UMAP with rsli integraiton
# =======================================================
DefaultAssay(All.atac) <- "peaks"
All.atac.split <- SplitObject(All.atac, split.by = "sampleID")

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = All.atac.split,
  anchor.features = All.atac@assays[["peaks"]]@var.features,
  reduction = "rlsi",
  dims = 2:30
)

# integrate LSI embeddings
All.atac.integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = All.atac@reductions[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

# create a new UMAP using the integrated embeddings
All.atac.integrated <- RunUMAP(All.atac.integrated, reduction = "integrated_lsi", dims = 2:30)

All.atac@reductions[["lsi_integrated"]] <- All.atac.integrated@reductions[["integrated_lsi"]]
All.atac@reductions[["umap_integrated"]] <- All.atac.integrated@reductions[["umap"]]

#Compare umap w/wo integration
p1 <- DimPlot(object = All.atac, group.by = 'sampleID', reduction = "umap") + ggtitle("Non-integrated")
p2 <- DimPlot(object = All.atac, group.by = 'sampleID', reduction = "umap_integrated") + ggtitle("Integrated")
p1|p2


#clustering on integrated lsi
All.atac <- FindNeighbors(object = All.atac, reduction = 'lsi_integrated', dims = 2:30, graph.name = c("integrated_nn","integrated_snn"))
All.atac <- FindClusters(object = All.atac, algorithm = 3, resolution = 0.8, graph.name = "integrated_snn") # resolution default = 0.8
DimPlot(object = All.atac, label = TRUE, reduction = 'umap_integrated') + NoLegend()

saveRDS(All.atac, file = "All.atac.rds")

#===============================================================================
#Create gene activity matrix 
# quantify gene activity
gene.activities <- GeneActivity(All.atac, extend.upstream = 2000, extend.downstream = 0, verbose = T) # default: upstream 2kb, downstream 0kb

# add the gene activity matrix to the Seurat object as a new assay
All.atac[['GeneActivity']] <- CreateAssayObject(counts = gene.activities)

# Normalize and scale gene activity matrix

DefaultAssay(All.atac) <- "GeneActivity"
All.atac <- NormalizeData(
  object = All.atac,
  assay = 'GeneActivity',
  normalization.method = 'LogNormalize',
  scale.factor = median(All.atac$nCount_GeneActivity)
)
All.atac <- ScaleData(All.atac, features = rownames(All.atac))


# =======================================================
# Integration with scRNA-seq data
# =======================================================

## Use "harmony_Louvain_res.0.8" clustering. Remove 8, 14, 15 (high pct.mt)
All.rna <- readRDS(paste0(primary_folder,"scrnRNA/20220331_All.rna.flt.rds"))
Idents(All.rna) <- "harmony_Louvain_res.0.8"
All.rna.flt <- subset(All.rna, harmony_Louvain_res.0.8 %in% c("8","14","15"), invert = TRUE)

# Assign cell type manually
All.rna.flt <- RenameIdents(
  object = All.rna.flt,
  '0' = 'NSC_2_astrocyte_like',
  '1' = 'iPS_1_neuroepithelial',
  '2' = 'iPS_0_undifferentiaed',
  '3' = 'neuron_mature',
  '4' = 'iPS_1_neuroepithelial',
  '5' = 'neuron_immature',
  '6' = 'neuron_mature',
  '7' = 'neuron_immature',
  #  '8' = 'pct_mt_high',
  '9' = 'NSC_1_differentiating',
  '10' = 'neuron_mature_2',
  '11' = 'OPC_like',
  '12' = 'NSC_0_iPS_like',
  '13' = 'neuron_progenitor',
  #  '14' = 'pct_mt_high',
  #  '15' = 'pct_mt_high',
  '16' = 'NSC_1_differentiating',
  '17' = 'NSC_0_iPS_like',
  '18' = 'neuron_schwann_like'
)

All.rna.flt$celltype_manual <- Idents(All.rna.flt)


#Visualize

DimPlot(All.rna, reduction = "umap_harmony", label = TRUE, repel = TRUE, group.by = 'harmony_Louvain_res.0.8') 
DimPlot(All.rna.flt, reduction = "umap_harmony", label = TRUE, repel = TRUE, group.by = 'harmony_Louvain_res.0.8') 
DimPlot(All.rna.flt, reduction = "umap_harmony", label = TRUE, repel = TRUE, group.by = 'celltype_manual')

# =======================================================
# Find anchors and transfer label to scATAC
# =======================================================
DefaultAssay(All.rna.flt) <- 'RNA'

transfer.anchors <- FindTransferAnchors(reference = All.rna.flt, 
                                        reference.assay = "RNA",
                                        normalization.method = "LogNormalize",
                                        features = VariableFeatures(object = All.rna.flt), 
                                        query = All.atac, 
                                        query.assay = "GeneActivity", 
                                        reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, 
                                     refdata = All.rna.flt$celltype_manual, 
                                     weight.reduction = All.atac[["lsi_integrated"]],
                                     dims = 2:30)

# =======================================================
# Add predicted label back to ATAC object
# =======================================================
All.atac$celltype_manual_LT <- factor(celltype.predictions$predicted.id)

#Visualize

DimPlot(All.rna.flt, reduction = "umap_harmony", label = TRUE, repel = TRUE, group.by = 'celltype_manual') + ggtitle("scRNA")
DimPlot(All.atac, reduction = "umap_integrated", group.by = 'celltype_manual_LT') + ggtitle("scATAC (LT from scRNA)")
DimPlot(All.atac, reduction = "umap_integrated", label = TRUE, repel = TRUE, group.by = 'integrated_snn_res.0.8') + ggtitle("scATAC clusters")

saveRDS(All.atac, file = "All.atac.rds")

#only the raw & final R object is located in the [primary_folder]/Data_and_code/DATA/, intermediate objects were not included


