#### Regeneration script: basic network node/edge composition across variants ####
#
# Feeds: Fig5/ext_fig5c_network_nodes.R, Fig5/ext_fig5c_network_edges.R
#
# For each of the 8 candidate network variants (base vs elastic-net-
# regularised, with/without TOBIAS footprinting, with/without Hi-C contacts),
# after applying the significance filter, counts distinct nodes (TFs,
# enhancers, promoters, split further by aCRE/tCRE measurement type) and
# edges (by link type: TF-to-peak, TF-to-promoter, TF-to-enhancer,
# enhancer-to-promoter, further split by measurement-type combination). This
# characterizes how much of the network's structure survives each
# combination of upstream modelling choices, independent of any downstream
# validation against ground truth.
#
# Reads (see ../config.R to set these paths):
#   [primary_data_folder]/3_preparing_base_network/base_network*.rds (4 files)
#   [primary_data_folder]/5_lasso_models/lasso_network*_with_randomization.rds (4 files)
#
# Writes to ../../Fig5/data/:
#   network_node_counts.csv  (network, metric, value -- num_tfs/num_enhancers/num_promoters)
#   network_edge_counts.csv  (network, metric, value -- the 5 edge-type metrics)
#   both already display-labelled (network as ordered factor levels, as text)
#   and metric-labelled, ready for direct ggplot use.

rm(list = ls())
library(dplyr)
library(tidyr)

# ── CONTROL PANEL ─────────────────────────────────────────────────────────────
source('../config.R')  # defines primary_data_folder

networks_to_compare <- c(
  'base', 'base_no_tobias', 'base_no_hic', 'base_no_hic_no_tobias',
  'lasso', 'lasso_no_tobias', 'lasso_no_hic', 'lasso_no_hic_no_tobias'
)
networks_files <- file.path(primary_data_folder, c(
  '3_preparing_base_network/base_network.rds',
  '3_preparing_base_network/base_network_no_footprinting.rds',
  '3_preparing_base_network/base_network_no_hic.rds',
  '3_preparing_base_network/base_network_no_hic_no_footprinting.rds',
  '5_lasso_models/lasso_network_with_randomization.rds',
  '5_lasso_models/lasso_network_no_footprinting_with_randomization.rds',
  '5_lasso_models/lasso_network_no_hic_with_randomization.rds',
  '5_lasso_models/lasso_network_no_hic_no_footprinting_with_randomization.rds'
))
networks_colors <- c(
  base                   = '#6BAED6',
  base_no_tobias         = '#FC8D59',
  base_no_hic            = '#74C476',
  base_no_hic_no_tobias  = '#BCBDDC',
  lasso                  = '#08519C',
  lasso_no_tobias        = '#B30000',
  lasso_no_hic           = '#238B45',
  lasso_no_hic_no_tobias = '#6A51A3'
)
significance_threshold <- 0.05
out_dir <- '../../Fig5/data'
# ──────────────────────────────────────────────────────────────────────────────

n_networks <- length(networks_files)

network_labels <- c(
  base                   = 'Base',
  base_no_tobias         = 'Base, no TOBIAS',
  base_no_hic            = 'Base, no Hi-C',
  base_no_hic_no_tobias  = 'Base, no Hi-C, no TOBIAS',
  lasso                  = 'Elastic net',
  lasso_no_tobias        = 'Elastic net, no TOBIAS',
  lasso_no_hic           = 'Elastic net, no Hi-C',
  lasso_no_hic_no_tobias = 'Elastic net, no Hi-C, no TOBIAS'
)
network_levels  <- unname(network_labels)
display_colors  <- setNames(networks_colors, network_labels[names(networks_colors)])

res <- data.frame(
  network                               = networks_to_compare,
  num_tfs                               = NA_integer_,
  num_enhancers                         = NA_integer_,
  num_promoters                         = NA_integer_,
  num_enhancers_aCRE                    = NA_integer_,
  num_enhancers_tCRE                    = NA_integer_,
  num_promoters_aCRE                    = NA_integer_,
  num_promoters_tCRE                    = NA_integer_,
  num_edges                             = NA_integer_,
  num_edges_TF_peak                     = NA_integer_,
  num_edges_TF_promoter                 = NA_integer_,
  num_edges_TF_promoter_aCRE            = NA_integer_,
  num_edges_TF_promoter_tCRE            = NA_integer_,
  num_edges_TF_enhancer                 = NA_integer_,
  num_edges_TF_enhancer_aCRE            = NA_integer_,
  num_edges_TF_enhancer_tCRE            = NA_integer_,
  num_edges_enhancer_promoter           = NA_integer_,
  num_edges_enhancer_aCRE_promoter_aCRE = NA_integer_,
  num_edges_enhancer_aCRE_promoter_tCRE = NA_integer_,
  num_edges_enhancer_tCRE_promoter_aCRE = NA_integer_,
  num_edges_enhancer_tCRE_promoter_tCRE = NA_integer_
)

for (i in seq_len(n_networks)) {

  message('Processing network: ', networks_to_compare[i])
  n <- readRDS(networks_files[i])

  if ('adjpvalue' %in% colnames(n)) {
    n <- n[n$adjpvalue <= significance_threshold, ]
  }

  res$num_tfs[i] <- length(unique(n$source[n$link_type == 'TF_peak']))
  res$num_enhancers[i] <- length(union(
    unique(n$source[n$source_type == 'Enhancer']),
    unique(n$target[n$target_type == 'Enhancer'])
  ))
  res$num_enhancers_aCRE[i] <- length(union(
    unique(n$source[n$source_type == 'Enhancer' & n$source_measurement_type == 'aCRE']),
    unique(n$target[n$target_type == 'Enhancer' & n$target_measurement_type == 'aCRE'])
  ))
  res$num_enhancers_tCRE[i] <- length(union(
    unique(n$source[n$source_type == 'Enhancer' & n$source_measurement_type == 'tCRE']),
    unique(n$target[n$target_type == 'Enhancer' & n$target_measurement_type == 'tCRE'])
  ))
  res$num_promoters[i] <- length(union(
    unique(n$source[n$source_type == 'Promoter']),
    unique(n$target[n$target_type == 'Promoter'])
  ))
  res$num_promoters_aCRE[i] <- length(union(
    unique(n$source[n$source_type == 'Promoter' & n$source_measurement_type == 'aCRE']),
    unique(n$target[n$target_type == 'Promoter' & n$target_measurement_type == 'aCRE'])
  ))
  res$num_promoters_tCRE[i] <- length(union(
    unique(n$source[n$source_type == 'Promoter' & n$source_measurement_type == 'tCRE']),
    unique(n$target[n$target_type == 'Promoter' & n$target_measurement_type == 'tCRE'])
  ))

  res$num_edges[i]                             <- nrow(n)
  res$num_edges_TF_peak[i]                     <- sum(n$link_type == 'TF_peak')
  res$num_edges_TF_promoter[i]                 <- sum(n$link_type == 'TF_peak' & n$target_type == 'Promoter')
  res$num_edges_TF_promoter_aCRE[i]            <- sum(n$link_type == 'TF_peak' & n$target_type == 'Promoter' & n$target_measurement_type == 'aCRE')
  res$num_edges_TF_promoter_tCRE[i]            <- sum(n$link_type == 'TF_peak' & n$target_type == 'Promoter' & n$target_measurement_type == 'tCRE')
  res$num_edges_TF_enhancer[i]                 <- sum(n$link_type == 'TF_peak' & n$target_type == 'Enhancer')
  res$num_edges_TF_enhancer_aCRE[i]            <- sum(n$link_type == 'TF_peak' & n$target_type == 'Enhancer' & n$target_measurement_type == 'aCRE')
  res$num_edges_TF_enhancer_tCRE[i]            <- sum(n$link_type == 'TF_peak' & n$target_type == 'Enhancer' & n$target_measurement_type == 'tCRE')
  res$num_edges_enhancer_promoter[i]           <- sum(n$link_type == 'Enhancer_promoter')
  res$num_edges_enhancer_aCRE_promoter_aCRE[i] <- sum(n$link_type == 'Enhancer_promoter' & n$source_measurement_type == 'aCRE' & n$target_measurement_type == 'aCRE')
  res$num_edges_enhancer_aCRE_promoter_tCRE[i] <- sum(n$link_type == 'Enhancer_promoter' & n$source_measurement_type == 'aCRE' & n$target_measurement_type == 'tCRE')
  res$num_edges_enhancer_tCRE_promoter_aCRE[i] <- sum(n$link_type == 'Enhancer_promoter' & n$source_measurement_type == 'tCRE' & n$target_measurement_type == 'aCRE')
  res$num_edges_enhancer_tCRE_promoter_tCRE[i] <- sum(n$link_type == 'Enhancer_promoter' & n$source_measurement_type == 'tCRE' & n$target_measurement_type == 'tCRE')
}

res_long <- res %>% pivot_longer(cols = -network, names_to = 'metric', values_to = 'value')

apply_network_factor <- function(df) {
  df$network <- factor(dplyr::recode(df$network, !!!network_labels), levels = network_levels)
  df
}

#### Figure 1 subset — node counts ####

node_metrics <- c('num_tfs', 'num_enhancers', 'num_promoters')
node_labels  <- c(
  num_tfs       = 'Number of TFs',
  num_enhancers = 'Number of enhancers',
  num_promoters = 'Number of promoters'
)

node_out <- res_long %>%
  filter(metric %in% node_metrics) %>%
  apply_network_factor() %>%
  mutate(metric = factor(dplyr::recode(metric, !!!node_labels), levels = unname(node_labels)),
         colour = display_colors[as.character(network)])

#### Figure 2 subset — edge counts ####

edge_metrics <- c('num_edges', 'num_edges_TF_peak',
                  'num_edges_TF_promoter', 'num_edges_TF_enhancer',
                  'num_edges_enhancer_promoter')
edge_labels  <- c(
  num_edges                   = 'Number of edges',
  num_edges_TF_peak           = 'Number of TF-to-peak edges',
  num_edges_TF_promoter       = 'Number of TF-to-promoter edges',
  num_edges_TF_enhancer       = 'Number of TF-to-enhancer edges',
  num_edges_enhancer_promoter = 'Number of enhancer-to-promoter edges'
)

edge_out <- res_long %>%
  filter(metric %in% edge_metrics) %>%
  apply_network_factor() %>%
  mutate(metric = factor(dplyr::recode(metric, !!!edge_labels), levels = unname(edge_labels)),
         colour = display_colors[as.character(network)])

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(node_out, row.names = FALSE, file.path(out_dir, 'network_node_counts.csv'))
write.csv(edge_out, row.names = FALSE, file.path(out_dir, 'network_edge_counts.csv'))
message('Saved: ', file.path(out_dir, 'network_node_counts.csv'), ' (', nrow(node_out), ' rows)')
message('Saved: ', file.path(out_dir, 'network_edge_counts.csv'), ' (', nrow(edge_out), ' rows)')
