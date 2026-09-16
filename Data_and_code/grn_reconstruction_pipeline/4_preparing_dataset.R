#### Assemble the expression/accessibility matrix for lasso modelling ####
#
# Step 4 of the GRN reconstruction pipeline. Restricts the paired scRNA-seq /
# snATAC-seq metacell object to exactly the measurements the base network
# (3_preparing_base_network.R) needs as predictors or targets — gene
# expression for every TF, chromatin accessibility (aCRE) and transcriptional
# activity (tCRE) for every peak involved in a candidate edge — normalises
# each modality, and stacks them into one matrix (rows = one measurement per
# gene/peak/modality, columns = metacells ordered by pseudotime). This matrix
# is what 5_lasso_models.R regresses each target measurement against.
#
# Reads (see ../config.R to set raw_data_folder and primary_data_folder):
#   [raw_data_folder]/metacell.combined.feature.Rds
#     paired scRNA-seq / snATAC-seq Seurat/Signac object, metacells as
#     columns, with RNA, aCRE, and tCRE assays
#   [raw_data_folder]/pseudotime.scRNA.tsv
#     per-metacell pseudotime (column SEACell = metacell id)
#   [primary_data_folder]/3_preparing_base_network/base_network(_no_hic)
#     (_no_footprinting).rds  (from 3_preparing_base_network.R; must match
#     this script's own use_hic / use_footprinting flags)
#
# Writes into output/4_preparing_dataset/:
#   dataset[_no_hic][_no_footprinting].rds
#     matrix, rows named <gene_or_peak>_Expression / _aCRE / _tCRE, columns
#     = metacells ordered by increasing pseudotime

rm(list = ls())
library(Matrix)
library(Seurat)
library(Signac)

# ── CONTROL PANEL ─────────────────────────────────────────────────────────────
source('../config.R')  # defines raw_data_folder, primary_data_folder

use_footprinting <- TRUE  # must match the base_network variant read below
use_hic <- TRUE           # must match the base_network variant read below

base_network_file <- file.path(primary_data_folder, '3_preparing_base_network',
                               paste0('base_network',
                                     ifelse(use_hic, '', '_no_hic'),
                                     ifelse(use_footprinting, '', '_no_footprinting'),
                                     '.rds'))
results_folder <- file.path(primary_data_folder, '4_preparing_dataset')
# ──────────────────────────────────────────────────────────────────────────────

dir.create(results_folder, showWarnings = FALSE, recursive = TRUE)

metacells    <- readRDS(file.path(raw_data_folder, 'metacell.combined.feature.Rds'))
base_network <- readRDS(base_network_file)
pseudo_time  <- read.table(file.path(raw_data_folder, 'pseudotime.scRNA.tsv'), header = TRUE)

pseudo_time <- pseudo_time[pseudo_time$SEACell %in% colnames(metacells), ]

metacells <- NormalizeData(metacells, assay = 'aCRE')
metacells <- NormalizeData(metacells, assay = 'tCRE')

# metacells ordered by pseudotime, each modality restricted to measurements
# with any non-zero signal
rna_data <- as.matrix(metacells@assays$RNA@data[, pseudo_time$SEACell])
rna_data <- rna_data[rowSums(rna_data) > 0, ]

aCRE_data <- as.matrix(metacells@assays$aCRE@data[, pseudo_time$SEACell])
aCRE_data <- aCRE_data[rowSums(aCRE_data) > 0, ]
rownames(aCRE_data) <- gsub('-', '_', rownames(aCRE_data))

tCRE_data <- as.matrix(metacells@assays$tCRE@data[, pseudo_time$SEACell])
tCRE_data <- tCRE_data[rowSums(tCRE_data) > 0, ]
rownames(tCRE_data) <- gsub('-', '_', rownames(tCRE_data))

# restrict to exactly what the base network needs
required_genes <- unique(base_network$source[base_network$source_type == 'TF'])
required_CRE   <- unique(c(base_network$source[base_network$source_type == 'Enhancer'],
                           base_network$target))

rna_data  <- rna_data[required_genes[required_genes %in% rownames(rna_data)], ]
aCRE_data <- aCRE_data[required_CRE[required_CRE %in% rownames(aCRE_data)], ]
tCRE_data <- tCRE_data[required_CRE[required_CRE %in% rownames(tCRE_data)], ]

rownames(rna_data)  <- paste0(rownames(rna_data),  '_Expression')
rownames(aCRE_data) <- paste0(rownames(aCRE_data), '_aCRE')
rownames(tCRE_data) <- paste0(rownames(tCRE_data), '_tCRE')
dataset <- rbind(rna_data, aCRE_data, tCRE_data)

message('Dataset assembled: ', nrow(dataset), ' measurements x ', ncol(dataset), ' metacells',
        ' (', nrow(rna_data), ' TF expression, ', nrow(aCRE_data), ' aCRE, ', nrow(tCRE_data), ' tCRE)')

saveRDS(dataset, file.path(results_folder,
                           paste0('dataset',
                                  ifelse(use_hic, '', '_no_hic'),
                                  ifelse(use_footprinting, '', '_no_footprinting'),
                                  '.rds')))
