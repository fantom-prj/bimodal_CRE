#### Regeneration script: per-enhancer centrality statistics ####
#
# Feeds (via copy_enhancer_stats.R -> Fig6/data/): Fig6d, Fig6e, ExtFig6e,
# ExtFig6f, and Data_and_code/screen_correlation/.
#
# For every enhancer in the filtered GRN, computes two centrality axes:
#   triplet_pairs           unique (TF, promoter) pairs bridged via the
#                           enhancer, i.e. distinct TF -> enhancer -> promoter
#                           paths it mediates
#   cum_abs_estimate_total  cumulative |estimate| over all edges touching the
#                           enhancer (outgoing + incoming)
# plus per-enhancer annotations, and flags the first Pareto front on the two
# axes (enhancers not dominated on both at once).
#
# Network filtering: adjpvalue, dev_ratio and per-stratum |estimate| p90
# (see ../network_loading/load_and_filter_network.R). With positive_only =
# TRUE, only activating edges (estimate > 0) are kept; the figures use the
# all-edges analysis (positive_only = FALSE).
#
# Reads (see ../config.R to set these paths):
#   [primary_data_folder]/5_lasso_models/lasso_network_with_randomization.rds
#   [primary_data_folder]/4_preparing_dataset/dataset.rds
#     (an enhancer counts as expressed if it has a tCRE measurement row)
#   ../network_loading/network_thresholds.txt
#
# Writes [primary_data_folder]/10_enhancer_centrality/<all_edges|positive_only>/
#   all_enhancers_stats.csv, one row per enhancer, sorted by triplet_pairs
#   (descending). Columns:
#     enhancer, triplet_pairs, cum_abs_estimate_out, cum_abs_estimate_in,
#     n_regulating_TFs, n_regulated_genes, n_regulated_promoters,
#     source_measurement_types (aCRE/tCRE modalities of the enhancer's
#       outgoing edges), is_expressed, is_mediating (on >= 1 complete
#       TF -> enhancer -> promoter path), cum_abs_estimate_total,
#     pareto (on the first Pareto front), label (enhancer id for the
#       enhancers labelled in centrality scatter plots, else empty)

rm(list = ls())
library(dplyr)

# ── CONTROL PANEL ─────────────────────────────────────────────────────────────
source('../config.R')  # defines primary_data_folder

positive_only            <- FALSE  # TRUE = activating edges only; FALSE = all edges
apply_estimate_threshold <- TRUE   # apply per-stratum p90 |estimate| filter
n_top_labels             <- 20     # enhancers given a non-empty label

network_file    <- file.path(primary_data_folder, '5_lasso_models/lasso_network_with_randomization.rds')
dataset_file    <- file.path(primary_data_folder, '4_preparing_dataset/dataset.rds')
thresholds_file <- '../network_loading/network_thresholds.txt'
out_folder      <- file.path(primary_data_folder, '10_enhancer_centrality',
                             if (positive_only) 'positive_only' else 'all_edges')
# ──────────────────────────────────────────────────────────────────────────────

dir.create(out_folder, showWarnings = FALSE, recursive = TRUE)

source('../network_loading/load_and_filter_network.R')
thresholds <- read_thresholds(thresholds_file)

#### Load and filter network ####

message('Loading network...')
n <- readRDS(network_file)
n <- filter_network(n, thresholds, apply_estimate_threshold = apply_estimate_threshold)
message('  Edges after filtering: ', nrow(n))

if (positive_only) {
  n <- n[n$estimate > 0, ]
  message('  Edges after positive-only filter: ', nrow(n))
}

#### Expressed and mediating enhancers ####

message('Loading dataset...')
dataset         <- readRDS(dataset_file)
tCRE_rows       <- rownames(dataset)[grepl('_tCRE$', rownames(dataset))]
expressed_peaks <- sub('_tCRE$', '', tCRE_rows)
rm(dataset)

all_enhancers <- union(
  n$source[n$source_type == 'Enhancer'],
  n$target[n$target_type == 'Enhancer']
)

# mediating: target of >= 1 TF_peak edge AND source of >= 1 Enhancer_promoter edge
tf_peak_targets     <- unique(n$target[n$link_type == 'TF_peak' & n$target_type == 'Enhancer'])
ep_sources          <- unique(n$source[n$link_type == 'Enhancer_promoter'])
mediating_enhancers <- intersect(tf_peak_targets, ep_sources)

message('  Enhancers in filtered network: ', length(all_enhancers))
message('  Mediating enhancers: ', length(mediating_enhancers))

#### Per-enhancer statistics ####

enhancer_stats <- function(net, enhancer_ids) {

  tf_to_enh <- net %>%
    filter(link_type == 'TF_peak', target %in% enhancer_ids) %>%
    select(TF = source, enhancer = target)

  enh_to_prom <- net %>%
    filter(link_type == 'Enhancer_promoter', source %in% enhancer_ids) %>%
    select(enhancer = source, promoter = target)

  triplet_counts <- inner_join(tf_to_enh, enh_to_prom, by = 'enhancer',
                               relationship = 'many-to-many') %>%
    group_by(enhancer) %>%
    summarise(triplet_pairs = n_distinct(paste(TF, promoter, sep = '__')),
              .groups = 'drop')

  abs_est_out <- net %>%
    filter(source %in% enhancer_ids) %>%
    group_by(source) %>%
    summarise(cum_abs_estimate_out = sum(abs(estimate), na.rm = TRUE), .groups = 'drop') %>%
    rename(enhancer = source)

  abs_est_in <- net %>%
    filter(target %in% enhancer_ids) %>%
    group_by(target) %>%
    summarise(cum_abs_estimate_in = sum(abs(estimate), na.rm = TRUE), .groups = 'drop') %>%
    rename(enhancer = target)

  n_tfs <- net %>%
    filter(target %in% enhancer_ids, link_type == 'TF_peak') %>%
    group_by(target) %>%
    summarise(n_regulating_TFs = n_distinct(source), .groups = 'drop') %>%
    rename(enhancer = target)

  n_genes <- net %>%
    filter(source %in% enhancer_ids, link_type == 'Enhancer_promoter') %>%
    group_by(source) %>%
    summarise(n_regulated_genes = n_distinct(target_gene), .groups = 'drop') %>%
    rename(enhancer = source)

  n_promoters <- net %>%
    filter(source %in% enhancer_ids, link_type == 'Enhancer_promoter') %>%
    group_by(source) %>%
    summarise(n_regulated_promoters = n_distinct(target), .groups = 'drop') %>%
    rename(enhancer = source)

  mtype_source <- net %>%
    filter(source %in% enhancer_ids, link_type == 'Enhancer_promoter') %>%
    group_by(source) %>%
    summarise(source_measurement_types = paste(sort(unique(source_measurement_type)),
                                               collapse = ','),
              .groups = 'drop') %>%
    rename(enhancer = source)

  flags <- data.frame(
    enhancer     = enhancer_ids,
    is_expressed = enhancer_ids %in% expressed_peaks,
    is_mediating = enhancer_ids %in% mediating_enhancers
  )

  stats <- data.frame(enhancer = enhancer_ids) %>%
    left_join(triplet_counts, by = 'enhancer') %>%
    left_join(abs_est_out,    by = 'enhancer') %>%
    left_join(abs_est_in,     by = 'enhancer') %>%
    left_join(n_tfs,          by = 'enhancer') %>%
    left_join(n_genes,        by = 'enhancer') %>%
    left_join(n_promoters,    by = 'enhancer') %>%
    left_join(mtype_source,   by = 'enhancer') %>%
    left_join(flags,          by = 'enhancer')

  num_cols <- c('triplet_pairs', 'cum_abs_estimate_out', 'cum_abs_estimate_in',
                'n_regulating_TFs', 'n_regulated_genes', 'n_regulated_promoters')
  stats[num_cols] <- lapply(stats[num_cols], function(x) ifelse(is.na(x), 0, x))
  stats$cum_abs_estimate_total <- stats$cum_abs_estimate_out + stats$cum_abs_estimate_in

  stats
}

# TRUE for enhancers not strictly dominated on both axes by another enhancer
pareto_front <- function(df) {
  n_rows    <- nrow(df)
  is_pareto <- rep(TRUE, n_rows)
  for (i in 1:n_rows) {
    if (!is_pareto[i]) next
    dominated <- df$triplet_pairs          < df$triplet_pairs[i] &
                 df$cum_abs_estimate_total < df$cum_abs_estimate_total[i]
    is_pareto[dominated] <- FALSE
  }
  is_pareto
}

# labels: every Pareto-front enhancer, topped up with the highest
# triplet_pairs enhancers if the front has fewer than n_top; if the front is
# larger, only its n_top highest on triplet_pairs + cum_abs_estimate_total
add_labels <- function(stats, n_top = n_top_labels) {
  label_ids <- stats$enhancer[stats$pareto]
  if (length(label_ids) < n_top) {
    extra <- stats %>%
      filter(!pareto) %>%
      arrange(desc(triplet_pairs)) %>%
      slice_head(n = n_top - length(label_ids)) %>%
      pull(enhancer)
    label_ids <- c(label_ids, extra)
  } else {
    label_ids <- stats %>%
      filter(pareto) %>%
      arrange(desc(triplet_pairs + cum_abs_estimate_total)) %>%
      slice_head(n = n_top) %>%
      pull(enhancer)
  }
  stats$label <- ifelse(stats$enhancer %in% label_ids, stats$enhancer, '')
  stats
}

message('Computing enhancer statistics...')
s_all <- enhancer_stats(n, all_enhancers)
s_all$pareto <- pareto_front(s_all)
s_all <- add_labels(s_all)

out_file <- file.path(out_folder, 'all_enhancers_stats.csv')
write.csv(s_all %>% arrange(desc(triplet_pairs)), row.names = FALSE, out_file)

message('Saved: ', out_file)
message('  Enhancers: ', format(nrow(s_all), big.mark = ','),
        '  expressed: ', sum(s_all$is_expressed),
        '  mediating: ', sum(s_all$is_mediating),
        '  Pareto front: ', sum(s_all$pareto))
message('Done.')
