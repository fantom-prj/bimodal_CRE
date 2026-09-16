#### Regeneration script: enhancer tCRE topology and distance-binned edges ####
#
# Produces the shared computation behind FOUR panels:
#   Fig6b     edges per enhancer by distance bin (with-tCRE network)
#   Fig6c     new-route donut (tCRE triplets: new route vs shared)
#   ExtFig6c  edges per enhancer by distance bin (without-tCRE, equal-footing control)
#   ExtFig6d  top-20 TFs / top-20 target genes by new-route triplet count
#
# Asks what modelling an enhancer's own transcriptional activity (tCRE) adds
# to the regulatory network, beyond what chromatin accessibility (aCRE) alone
# already captures:
#   Topology: for every complete TF -> enhancer(tCRE) -> promoter triplet,
#     flags whether that TF/promoter pair is also reachable through the
#     aCRE-only network ("shared") or only appears once tCRE is modelled
#     ("new route").
#   Distance-binned edge density: for both the with-tCRE and without-tCRE
#     networks, the mean number of lasso-selected enhancer-promoter edges per
#     enhancer in 50 kb distance bins, comparing expressed vs non-expressed
#     enhancers -- pool-size corrected so every base-network enhancer in a
#     bin counts, including those with zero selected edges.
#   (The plots themselves live in the four panel scripts.)
#
# Filtering (deliberately NOT the shared filter_network() helper used
# elsewhere in Data_and_code/): adjpvalue <= 0.05, then per-stratum |estimate|
# p90 RECOMPUTED LIVE within EACH network via apply_estimate_filter(), then
# dev_ratio >= 0.3. The live recompute is intentional: this script compares
# two structurally different networks (with vs without enhancer tCRE
# measurements), and applying one frozen threshold set derived from the
# with-tCRE network to both would impose inconsistent fractional stringency
# given their scale differences. Every OTHER script in this repository reads
# the frozen p90 values from network_thresholds.txt -- this one is the
# documented exception, do not "fix" it to match the others.
#
# Outputs (all into Fig6/data/):
#   distance_bin_with_tCRE.csv     - source_class, dist_bin, n_enhancers, y_mean
#                                    (feeds Fig6b)
#   distance_bin_without_tCRE.csv  - same shape, without-tCRE network
#                                    (feeds ExtFig6c)
#   distance_bin_promoters.csv     - distinct promoters reached per enhancer
#                                    (double-chance control; not used by a
#                                    requested panel, kept because it is the
#                                    third sub-analysis of Theme 5 and costs
#                                    nothing extra to compute)
#   triplets_tCRE.csv              - one row per TF x enhancer(tCRE) x promoter
#                                    triplet: enhancer, TF, promoter,
#                                    promoter_measurement_type, tf_prom_id,
#                                    new_route, gene
#                                    (feeds Fig6c and ExtFig6d; the `gene`
#                                    column is resolved here via
#                                    promoter_gene_lookup so ExtFig6d does not
#                                    have to rebuild it)
#
# Cross-check (printed by this script): 3,929 total triplets, 90.3% new-route,
# 9.7% shared; 3,098 new TF->enhancer edges, 3,390 new enhancer->promoter edges.

rm(list = ls())
library(tidyverse)

# ── CONTROL PANEL ─────────────────────────────────────────────────────────────
source('../config.R')  # defines primary_data_folder

network_with_tCRE_file    <- file.path(primary_data_folder,
  '5_lasso_models/lasso_network_with_randomization.rds')
network_without_tCRE_file <- file.path(primary_data_folder,
  '5_lasso_models/lasso_network_with_randomization_no_enhancer_tCRE.rds')
dataset_file              <- file.path(primary_data_folder,
  '4_preparing_dataset/dataset.rds')
base_network_file         <- file.path(primary_data_folder,
  '3_preparing_base_network/base_network.rds')

thresholds_file <- '../network_loading/network_thresholds.txt'
out_folder      <- '../../Fig6/data'

apply_estimate_threshold <- TRUE
bin_width                <- 50000   # 50 kb distance bins
# ──────────────────────────────────────────────────────────────────────────────

# read_thresholds() supplies adjpvalue / dev_ratio; the |estimate| p90 values it
# also returns are deliberately unused here (see the filtering note above).
source('../network_loading/load_and_filter_network.R')
thresholds          <- read_thresholds(thresholds_file)
adjpvalue_threshold <- thresholds$adjpvalue_threshold
dev_ratio_threshold <- thresholds$dev_ratio_threshold
message('Thresholds loaded from: ', thresholds_file)

#### Helper: per-stratum p90 |estimate| filter, recomputed per network ####

apply_estimate_filter <- function(net) {
  all_keys <- ifelse(net$link_type == 'TF_peak',
                     paste0('TF_', net$target_type, '_', net$target_measurement_type),
                     paste0('EP_', net$source_measurement_type, '_', net$target_measurement_type))
  unique_keys <- unique(all_keys)
  thr <- sapply(unique_keys, function(k) {
    unname(quantile(abs(net$estimate[all_keys == k]), 0.90, na.rm = TRUE))
  })
  thr_per_edge <- thr[all_keys]
  net[!is.na(thr_per_edge) & abs(net$estimate) >= thr_per_edge, ]
}

#### Load networks ####

message('Loading networks...')
n_with    <- readRDS(network_with_tCRE_file)
n_without <- readRDS(network_without_tCRE_file)

n_with    <- n_with   [n_with$adjpvalue    <= adjpvalue_threshold, ]
n_without <- n_without[n_without$adjpvalue <= adjpvalue_threshold, ]

if (apply_estimate_threshold) {
  n_with    <- apply_estimate_filter(n_with)
  n_without <- apply_estimate_filter(n_without)
}
message('  With-tCRE edges after filtering: ',    nrow(n_with))
message('  Without-tCRE edges after filtering: ', nrow(n_without))

n_with_dr    <- n_with   [n_with$dev_ratio    >= dev_ratio_threshold, ]
n_without_dr <- n_without[n_without$dev_ratio >= dev_ratio_threshold, ]

message('Loading dataset and base network...')
dataset      <- readRDS(dataset_file)
base_network <- readRDS(base_network_file)

#### Topology: tCRE triplets and the new-route flag ####

message('Topology...')

tf_enh_with    <- n_with_dr   [n_with_dr$link_type    == 'TF_peak' & n_with_dr$target_type    == 'Enhancer', ]
tf_enh_without <- n_without_dr[n_without_dr$link_type == 'TF_peak' & n_without_dr$target_type == 'Enhancer', ]
ep_with        <- n_with_dr   [n_with_dr$link_type    == 'Enhancer_promoter', ]
ep_without     <- n_without_dr[n_without_dr$link_type == 'Enhancer_promoter', ]

# TF->enhancer(tCRE) and enhancer(tCRE)->promoter edges gained
tf_enh_tCRE_specific <- tf_enh_with[tf_enh_with$target_measurement_type == 'tCRE', ]
n_tf_enh_tCRE_gained <- nrow(tf_enh_tCRE_specific)
ep_tCRE_specific     <- ep_with[ep_with$source_measurement_type == 'tCRE', ]
n_ep_tCRE_gained     <- nrow(ep_tCRE_specific)

tf_to_enh_tCRE <- tf_enh_with[tf_enh_with$target_measurement_type == 'tCRE',
                              c('source', 'target')]
colnames(tf_to_enh_tCRE) <- c('TF', 'enhancer')

enh_tCRE_to_prom <- ep_with[ep_with$source_measurement_type == 'tCRE',
                            c('source', 'target', 'target_measurement_type')]
colnames(enh_tCRE_to_prom) <- c('enhancer', 'promoter', 'promoter_measurement_type')

# peak -> gene lookup (used to resolve the `gene` column below)
promoter_gene_lookup <- unique(ep_with[, c('target', 'target_gene')])
promoter_gene_lookup <- setNames(promoter_gene_lookup$target_gene,
                                 promoter_gene_lookup$target)

triplets_tCRE   <- unique(merge(tf_to_enh_tCRE, enh_tCRE_to_prom, by = 'enhancer'))
n_triplets_tCRE <- nrow(triplets_tCRE)

# new route = TF/promoter pair absent from the aCRE-only (without-tCRE) network
tf_to_enh_aCRE <- tf_enh_without[tf_enh_without$target_measurement_type == 'aCRE',
                                 c('source', 'target')]
colnames(tf_to_enh_aCRE) <- c('TF', 'enhancer')
enh_aCRE_to_prom <- ep_without[ep_without$source_measurement_type == 'aCRE',
                               c('source', 'target')]
colnames(enh_aCRE_to_prom) <- c('enhancer', 'promoter')
triplets_aCRE <- unique(merge(tf_to_enh_aCRE, enh_aCRE_to_prom,
                              by = 'enhancer')[, c('TF', 'promoter')])

triplets_tCRE$tf_prom_id <- paste0(triplets_tCRE$TF, '__', triplets_tCRE$promoter)
aCRE_tf_prom_ids         <- paste0(triplets_aCRE$TF, '__', triplets_aCRE$promoter)
triplets_tCRE$new_route  <- !(triplets_tCRE$tf_prom_id %in% aCRE_tf_prom_ids)
frac_new_route           <- mean(triplets_tCRE$new_route)

# gene names for the promoter end, used by ExtFig6d's target-gene ranking
triplets_tCRE$gene <- promoter_gene_lookup[triplets_tCRE$promoter]

message('  New TF->enhancer(tCRE) edges: ', n_tf_enh_tCRE_gained)
message('  New enhancer(tCRE)->promoter edges: ', n_ep_tCRE_gained)
message('  Complete tCRE triplets: ', n_triplets_tCRE,
        '  (new route ', round(100 * frac_new_route, 1), '%, shared ',
        round(100 * (1 - frac_new_route), 1), '%)')

#### Pool-size-corrected edges per enhancer in distance bins ####

message('Distance-binned edges per enhancer...')

all_peaks_in_dataset       <- sub('_aCRE$', '',
  rownames(dataset)[grepl('_aCRE$', rownames(dataset))])
expressed_peaks_in_dataset <- sub('_tCRE$', '',
  rownames(dataset)[grepl('_tCRE$', rownames(dataset))])

all_enhancers_in_base <- union(
  base_network$source[base_network$source_type == 'Enhancer'],
  base_network$target[base_network$target_type == 'Enhancer']
)

expressed_in_base <- expressed_peaks_in_dataset[
  expressed_peaks_in_dataset %in% all_enhancers_in_base]

# Shared denominator: every base-network enhancer, its median EP distance, and
# whether it is expressed. Distance binned at bin_width up to 500 kb.
all_enh_dist_base <- base_network %>%
  filter(link_type == 'Enhancer_promoter') %>%
  group_by(source) %>%
  summarise(median_dist = median(link_score), .groups = 'drop') %>%
  mutate(expressed_source = source %in% expressed_in_base,
         dist_bin = cut(median_dist,
                        breaks = seq(0, 500000, by = bin_width),
                        include.lowest = TRUE, right = FALSE))

# Data half of the distance-bin summary only (the ggplot half lives in the
# panel scripts). Enhancers absent from enh_counts_df contribute 0, which is
# the whole point of the pool-size correction.
make_bin_summary <- function(enh_counts_df, all_enh_base_df, y_var) {
  full_df <- all_enh_base_df %>%
    left_join(enh_counts_df %>% select(source, !!sym(y_var)), by = 'source') %>%
    mutate(!!y_var := ifelse(is.na(!!sym(y_var)), 0L, !!sym(y_var)))

  full_df %>%
    filter(!is.na(dist_bin)) %>%
    group_by(dist_bin, expressed_source) %>%
    summarise(n_enhancers = n(),
              y_mean      = mean(!!sym(y_var)),
              .groups = 'drop') %>%
    mutate(source_class = ifelse(expressed_source,
                                 'expressed enhancers',
                                 'non-expressed enhancers')) %>%
    select(source_class, dist_bin, n_enhancers, y_mean)
}

# All EP edges in the with-tCRE network  (Fig6b)
enh_ep_counts <- ep_with %>%
  group_by(source) %>%
  summarise(n_ep_edges = n(), median_dist = median(link_score), .groups = 'drop')
bin_summary_with <- make_bin_summary(enh_ep_counts, all_enh_dist_base, 'n_ep_edges')

# Distinct promoters reached (double-chance control; not a panel)
enh_prom_counts <- ep_with %>%
  group_by(source) %>%
  summarise(n_promoters_reached = n_distinct(target),
            median_dist         = median(link_score), .groups = 'drop')
bin_summary_prom <- make_bin_summary(enh_prom_counts, all_enh_dist_base,
                                     'n_promoters_reached')

# aCRE-only edges in the without-tCRE network  (ExtFig6c)
enh_ep_counts_wout <- ep_without %>%
  group_by(source) %>%
  summarise(n_ep_edges = n(), median_dist = median(link_score), .groups = 'drop')
bin_summary_wout <- make_bin_summary(enh_ep_counts_wout, all_enh_dist_base,
                                     'n_ep_edges')

#### Save ####

dir.create(out_folder, showWarnings = FALSE, recursive = TRUE)

write.csv(bin_summary_with, row.names = FALSE,
          file.path(out_folder, 'distance_bin_with_tCRE.csv'))
write.csv(bin_summary_wout, row.names = FALSE,
          file.path(out_folder, 'distance_bin_without_tCRE.csv'))
write.csv(bin_summary_prom, row.names = FALSE,
          file.path(out_folder, 'distance_bin_promoters.csv'))
write.csv(triplets_tCRE, row.names = FALSE,
          file.path(out_folder, 'triplets_tCRE.csv'))

message('Saved: distance_bin_with_tCRE.csv (', nrow(bin_summary_with), ' rows)')
message('Saved: distance_bin_without_tCRE.csv (', nrow(bin_summary_wout), ' rows)')
message('Saved: distance_bin_promoters.csv (', nrow(bin_summary_prom), ' rows)')
message('Saved: triplets_tCRE.csv (', nrow(triplets_tCRE), ' rows)')
message('Done.')
