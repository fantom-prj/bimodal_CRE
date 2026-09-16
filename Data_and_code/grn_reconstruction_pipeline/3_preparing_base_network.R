#### Assemble the candidate (base) GRN edge set ####
#
# Step 3 of the GRN reconstruction pipeline. Combines the TF-to-peak calls
# from 2_tf_2_peak.R with genomic proximity and Hi-C chromatin-contact
# evidence to build a single candidate edge table spanning three edge
# classes, before any expression/accessibility-based pruning:
#   TF -> Promoter        (a TF-bound promoter peak)
#   TF -> Enhancer        (a TF-bound enhancer peak)
#   Enhancer -> Promoter  (candidate regulatory link between an enhancer and
#                          a promoter within ep_max_distance, either by raw
#                          genomic proximity or by a significant Hi-C contact)
#
# An enhancer is linked to a promoter if either: (a) it is within
# no_hic_filter_region of the promoter (regardless of Hi-C), or (b) it falls
# within ep_max_distance AND has a Hi-C contact with the promoter passing
# hic_qvalue_threshold (only when use_hic = TRUE). This candidate set is
# later pruned by elastic-net regression against expression/accessibility
# data (5_lasso_models.R); it does not yet reflect which edges are actually
# statistically supported.
#
# The `use_hic` / `use_footprinting` control-panel flags produce one of four
# base-network variants (all/neither/either data source included), used
# together to ask how much each data source contributes to the final
# network (see Data_and_code/network_basic_stats/ and
# Data_and_code/chipseq_eqtl_validation/, which compare all four variants
# against all four corresponding lasso-network variants from
# 5_lasso_models.R). Run this script once per combination of the two flags
# to produce every variant those comparisons need.
#
# Reads (see ../config.R to set raw_data_folder and primary_data_folder):
#   [raw_data_folder]/all_aCRE.PE.final.tsv
#     annotated peak table (peak type: enhancer-like / promoter-like /
#     other; gene symbol for promoter peaks)
#   [raw_data_folder]/<iPSC|NSC|NRN>_001_aCRE_HiC_contact5k.rds
#     per-cell-type Hi-C contact tables (aCRE-aCRE pairs with a q-value)
#   [primary_data_folder]/2_tf_2_peak/tf2peak(_no_footprinting).rds
#     (from 2_tf_2_peak.R)
#
# Writes into output/3_preparing_base_network/ (filename suffix reflects
# use_hic / use_footprinting, e.g. base_network_no_hic.rds when
# use_hic = FALSE):
#   base_network[_no_hic][_no_footprinting].rds
#     columns: source, target, source_type, target_type, source_gene,
#     target_gene, link_type, link_score (TOBIAS score for TF_peak edges;
#     genomic distance in bp for Enhancer_promoter edges)
#   base_network[_no_hic][_no_footprinting]_subset.rds / .csv
#     (100-row random sample, for quick review)
#   processed_hic.rds  (cached, distance- and q-value-filtered Hi-C
#                        contacts; shared across variants, reused on rerun)

rm(list = ls())
library(tidyverse)
library(doParallel)
library(foreach)
library(GenomicRanges)
library(data.table)

# ── CONTROL PANEL ─────────────────────────────────────────────────────────────
source('../config.R')  # defines raw_data_folder, primary_data_folder

peak_types_to_keep <- c('enhancer-like', 'promoter-like')
annotated_peaks_file <- file.path(raw_data_folder, 'all_aCRE.PE.final.tsv')
tf2peak_file <- file.path(primary_data_folder, '2_tf_2_peak/tf2peak')
ep_max_distance <- 500000  # maximum enhancer-promoter distance considered (bp)
use_footprinting <- TRUE  # rerun with FALSE for the "no_footprinting" variant
use_hic <- TRUE           # rerun with FALSE for the "no_hic" variant
hic_resolution <- '5k'
no_hic_filter_region <- 5000  # enhancers within this range are linked regardless of Hi-C
hic_qvalue_threshold <- 0.05
hic_folder <- raw_data_folder
hic_files <- c('iPSC', 'NSC', 'NRN')
num_cores <- 20

results_folder <- file.path(primary_data_folder, '3_preparing_base_network')
# ──────────────────────────────────────────────────────────────────────────────

dir.create(results_folder, showWarnings = FALSE, recursive = TRUE)

my.cluster <- parallel::makeCluster(num_cores)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()

#### loading data ####

annotated_peaks <- fread(annotated_peaks_file, data.table = FALSE)

if (use_footprinting) {
  tf2peak <- readRDS(paste0(tf2peak_file, '.rds'))
} else {
  tf2peak <- readRDS(paste0(tf2peak_file, '_no_footprinting.rds'))
}

tf2peak$start  <- tf2peak$start + 1
tf2peak$peakID <- paste0(tf2peak$seqnames, '_', tf2peak$start, '_', tf2peak$end)
tf2peak <- tf2peak[tf2peak$peakID %in% annotated_peaks$peakID, ]

# Hi-C contacts: cached after the first, expensive pass (distance + q-value
# filtering over all chromosome pairs)
if (!file.exists(file.path(results_folder, 'processed_hic.rds'))) {

  hic <- vector('list', length(hic_files))

  for (i in 1:length(hic_files)) {

    print(i)
    hic[[i]] <- readRDS(file.path(hic_folder,
                                  paste0(hic_files[i], '_001_aCRE_HiC_contact', hic_resolution, '.rds')))

    distance <- apply(hic[[i]], 1, function(x) {
      a <- as.numeric(strsplit(x[1], '-')[[1]][3])
      b <- as.numeric(strsplit(x[2], '-')[[1]][2])
      b - a - 1
    })
    hic[[i]]$distance <- distance
    idx <- distance <= ep_max_distance
    hic[[i]] <- hic[[i]][idx, ]

    print(paste0('percentage kept: ', sum(idx) / length(idx)))

    hic[[i]][[1]] <- gsub('-', '_', hic[[i]][[1]])
    hic[[i]][[2]] <- gsub('-', '_', hic[[i]][[2]])
  }

  hic <- do.call(rbind, hic)
  hic <- hic[hic$qvalue <= hic_qvalue_threshold, ]
  hic$qvalue <- NULL
  hic <- unique(hic)

  saveRDS(hic, file = file.path(results_folder, 'processed_hic.rds'))

} else {
  hic <- readRDS(file.path(results_folder, 'processed_hic.rds'))
}

#### filtering peaks to enhancer-like / promoter-like ####

annotated_peaks <- annotated_peaks[annotated_peaks$promoter_type_CT %in% peak_types_to_keep, ]

#### TF -> Promoter edges ####

to_keep <- which(tf2peak$peakID %in%
                   annotated_peaks$peakID[annotated_peaks$promoter_type_CT == 'promoter-like'])

tmp <- merge(tf2peak[to_keep, ], annotated_peaks[annotated_peaks$promoter_type_CT == 'promoter-like', ],
             by.x = 'peakID', by.y = 'peakID', all.x = TRUE)

base_network <- data.frame(source = tmp$tf,
                           target = tmp$peakID,
                           source_type = 'TF',
                           target_type = 'Promoter',
                           source_gene = tmp$tf,
                           target_gene = tmp$geneName,
                           link_type = 'TF_peak',
                           link_score = tmp$max_score)

base_network <- unique(base_network)

#### TF -> Enhancer edges ####

to_keep <- which(tf2peak$peakID %in%
                   annotated_peaks$peakID[annotated_peaks$promoter_type_CT == 'enhancer-like'])

tmp <- merge(tf2peak[to_keep, ], annotated_peaks[annotated_peaks$promoter_type_CT == 'enhancer-like', ],
             by.x = 'peakID', by.y = 'peakID', all.x = TRUE)

to_add <- data.frame(source = tmp$tf,
                     target = tmp$peakID,
                     source_type = 'TF',
                     target_type = 'Enhancer',
                     source_gene = tmp$tf,
                     target_gene = tmp$geneName,
                     link_type = 'TF_peak',
                     link_score = tmp$max_score)

to_add <- unique(to_add)
base_network <- rbind(base_network, to_add)

#### Enhancer -> Promoter edges ####

annotated_peaks$seqnames <- sapply(annotated_peaks$peakID, function(x) strsplit(x, split = '_')[[1]][1])
annotated_peaks$start <- as.numeric(sapply(annotated_peaks$peakID, function(x) strsplit(x, split = '_')[[1]][2]))
annotated_peaks$end   <- as.numeric(sapply(annotated_peaks$peakID, function(x) strsplit(x, split = '_')[[1]][3]))

annotated_peaks <- annotated_peaks[order(annotated_peaks$peakID), ]
annotated_promoters <- annotated_peaks[annotated_peaks$promoter_type_CT == 'promoter-like', ]
annotated_enhancers <- annotated_peaks[annotated_peaks$promoter_type_CT == 'enhancer-like', ]

tmp <- foreach(i = 1:dim(annotated_promoters)[1]) %dopar% {

  tryCatch({

    current_promoter <- annotated_promoters[i, ]

    # distance to every candidate enhancer (0 if overlapping, Inf if on a
    # different chromosome)
    enhancer_on_the_right_distances <- annotated_enhancers$start - current_promoter$end - 1
    enhancer_on_the_left_distances  <- current_promoter$start - annotated_enhancers$end - 1
    distances <- pmax(enhancer_on_the_right_distances, enhancer_on_the_left_distances)

    idx <- which(enhancer_on_the_left_distances * enhancer_on_the_right_distances > 0)
    if (length(idx) > 0) distances[idx] <- 0

    idx <- which(annotated_enhancers$seqnames != current_promoter$seqnames)
    if (length(idx) > 0) distances[idx] <- Inf

    idx <- which(distances <= ep_max_distance)
    if (length(idx) > 0) {
      current_peaks <- annotated_enhancers[idx, ]
      current_peaks$distance <- distances[idx]
    } else {
      current_peaks <- data.frame(matrix(NA, 0, ncol(annotated_enhancers) + 1))
      colnames(current_peaks) <- c(colnames(annotated_enhancers), 'distance')
    }

    if (use_hic) {

      # candidates within the "close enough regardless of Hi-C" range
      idx <- current_peaks$distance <= no_hic_filter_region
      current_peaks <- current_peaks[idx, ]
      current_peaks <- current_peaks[, c('peakID', 'distance')]

      # plus any enhancer with a significant Hi-C contact to this promoter
      current_peaks_1 <- hic[hic$aCRE_left  == current_promoter$peakID, c('aCRE_right', 'distance')]
      colnames(current_peaks_1) <- c('peakID', 'distance')
      current_peaks_2 <- hic[hic$aCRE_right == current_promoter$peakID, c('aCRE_left',  'distance')]
      colnames(current_peaks_2) <- c('peakID', 'distance')
      current_peaks <- rbind(current_peaks[, c('peakID', 'distance')],
                             unique(rbind(current_peaks_1, current_peaks_2)))
      if (dim(current_peaks)[1] == 0) stop('No peak within range')

      current_peaks <- current_peaks[current_peaks$peakID %in% annotated_enhancers$peakID, ]
      current_peaks <- merge(current_peaks, annotated_enhancers)
    }

    if (dim(current_peaks)[1] == 0) stop('No peak within range')

    to_add <- data.frame(source = current_peaks$peakID,
                         target = current_promoter$peakID,
                         source_type = 'Enhancer',
                         target_type = 'Promoter',
                         source_gene = current_peaks$geneName,
                         target_gene = current_promoter$geneName,
                         link_type = 'Enhancer_promoter',
                         link_score = current_peaks$distance)

    to_add

  }, error = function(x) paste0(print(x), ', error at iteration ', i))
}

errors <- which(sapply(tmp, class) == 'character')
if (length(errors) > 0) {
  sink(file.path(results_folder, 'error_messages_script_3.txt'))
  for (e in errors) print(tmp[[e]])
  sink()
  tmp <- tmp[-errors]
}

to_add <- do.call(rbind, tmp)
base_network <- rbind(base_network, to_add)

#### final formatting ####

base_network$source[base_network$source_type == 'TF'] <-
  toupper(base_network$source[base_network$source_type == 'TF'])
base_network$source_gene <- toupper(base_network$source_gene)
base_network$target_gene <- toupper(base_network$target_gene)

message('Base network edges: ', nrow(base_network),
        ' (', sum(base_network$link_type == 'TF_peak'), ' TF_peak, ',
        sum(base_network$link_type == 'Enhancer_promoter'), ' Enhancer_promoter)')

#### saving ####

variant_suffix <- paste0(ifelse(use_hic, '', '_no_hic'),
                         ifelse(use_footprinting, '', '_no_footprinting'))

filename <- file.path(results_folder, paste0('base_network', variant_suffix, '.rds'))
saveRDS(base_network, file = filename)

set.seed(12345)
ids <- sample(1:dim(base_network)[1], 100)
subset_filename <- file.path(results_folder, paste0('base_network', variant_suffix, '_subset'))
saveRDS(base_network[ids, ], file = paste0(subset_filename, '.rds'))
write.csv(base_network[ids, ], row.names = FALSE, file = paste0(subset_filename, '.csv'))

stopCluster(my.cluster)
