#### Regeneration script: threshold sensitivity sweep for TF-set coupling ####
#
# Produces the per-node Jaccard category table used by SIX panels:
#   Fig5d   stacked bar, TF_promoter x estimate_pctl  (linear bar height)
#   ExtFig5d stacked bar, TF_promoter x adjpvalue     (linear bar height)
#   ExtFig5e stacked bar, TF_promoter x dev_ratio     (linear bar height)
#   Fig6a   stacked bar, TF_enhancer x estimate_pctl  (linear bar height)
#   ExtFig6a stacked bar, TF_enhancer x adjpvalue     (linear bar height)
#   ExtFig6b stacked bar, TF_enhancer x dev_ratio     (linear bar height)
#
# Tests how robust the "TF-set decoupling" result is to the network's
# significance thresholds. For each node (an expressed enhancer regulated at
# both chromatin and expression level -- TF_enhancer; or a promoter regulated
# at both chromatin and expression level -- TF_promoter), the two governing
# TF sets are compared by Jaccard index and bucketed as fully uncoupled
# (Jaccard=0), intermediate, or fully coupled (Jaccard=1). This is repeated
# across a grid of threshold values for each of the three filtering
# parameters in turn, holding the other two fixed at baseline, to check that
# the reported coupling proportions are not an artefact of one arbitrary
# cutoff choice.
#
# Sweeps, one-at-a-time, holding the other two thresholds at the
# network_thresholds.txt baseline (adjpvalue <= 0.05, dev_ratio >= 0.3,
# |estimate| >= per-stratum p90):
#   adjpvalue_threshold  : 0.02, 0.05*, 0.1, 0.2
#   dev_ratio_threshold  : 0.1, 0.2, 0.3*, 0.4, 0.5
#   |estimate| percentile: 75, 80, 90*, 95
# = 13 grid points; 13 x 2 themes x 3 categories = 78 category rows.
#
# IMPORTANT (do not "simplify"):
#   * Filter order is adjpvalue -> dev_ratio -> |estimate| percentile. The
#     percentile MUST be applied last, on edges that already passed both
#     other filters -- that is how the stored p90 values in
#     network_thresholds.txt were themselves derived. Swapping the last two
#     silently inflates every sample size, since the percentile would then
#     be computed over a different (larger) edge set than the one it is
#     meant to threshold.
#   * apply_estimate_filter_pctl() recomputes tapply()+quantile() per stratum
#     LIVE at every grid point, INCLUDING pctl = 90, instead of reading the
#     stored (3-decimal-rounded) p90 from network_thresholds.txt. This is by
#     design: the percentile is a swept axis here, so p90 is just one point
#     on it. It makes the pctl=90 point differ slightly from a fixed-
#     threshold computation using the rounded stored values (e.g.
#     TF_promoter n = 2457 vs 2466) -- the gap is entirely explained by that
#     3-decimal rounding. The published figures are sourced from this sweep,
#     all four percentile points, so this live recompute is the correct
#     behaviour to reproduce.
#
# Outputs (written to BOTH figure data folders, since the six panels that
# use them live in two different figure folders):
#   Fig5/data/sweep_category_table.csv
#   Fig6/data/sweep_category_table.csv
#     columns: sweep, sweep_value, theme, category, pct, n, n_nodes
#     (n_nodes is the sample size of the grid point x theme; it sets the bar
#      height in the linear-height stacked bars)

rm(list = ls())
library(dplyr)

# ── CONTROL PANEL ─────────────────────────────────────────────────────────────
source('../config.R')  # defines primary_data_folder

network_file      <- file.path(primary_data_folder,
                               '5_lasso_models/lasso_network_with_randomization.rds')
dataset_file      <- file.path(primary_data_folder,
                               '4_preparing_dataset/dataset.rds')
base_network_file <- file.path(primary_data_folder,
                               '3_preparing_base_network/base_network.rds')

thresholds_file   <- '../network_loading/network_thresholds.txt'
out_folders       <- c('../../Fig5/data', '../../Fig6/data')

adjpvalue_grid     <- c(0.02, 0.05, 0.1, 0.2)
dev_ratio_grid     <- c(0.1, 0.2, 0.3, 0.4, 0.5)
estimate_pctl_grid <- c(75, 80, 90, 95)
# ──────────────────────────────────────────────────────────────────────────────

# read_thresholds() -- only the baseline adjpvalue / dev_ratio are used here;
# the |estimate| threshold is recomputed per grid point (see header note).
source('../network_loading/load_and_filter_network.R')

thresholds         <- read_thresholds(thresholds_file)
baseline_adjpvalue <- thresholds$adjpvalue_threshold
baseline_dev_ratio <- thresholds$dev_ratio_threshold
baseline_pctl      <- 90

message('Baseline: adjpvalue <= ', baseline_adjpvalue,
        ', dev_ratio >= ', baseline_dev_ratio,
        ', |estimate| >= p', baseline_pctl, ' per stratum')

#### Load network, dataset, base network (once) ####

message('Loading network, dataset, base network...')
net_full     <- readRDS(network_file)
dataset      <- readRDS(dataset_file)
base_network <- readRDS(base_network_file)

expressed_peaks_in_dataset <- sub('_tCRE$', '',
  rownames(dataset)[grepl('_tCRE$', rownames(dataset))])
all_enhancers_in_base <- union(
  base_network$source[base_network$source_type == 'Enhancer'],
  base_network$target[base_network$target_type == 'Enhancer']
)
expressed_enhancers <- expressed_peaks_in_dataset[
  expressed_peaks_in_dataset %in% all_enhancers_in_base
]
message('  Expressed enhancers in base network: ', length(expressed_enhancers))

#### Helpers ####

jaccard <- function(a, b) {
  u <- length(union(a, b))
  if (u == 0L) return(NA_real_)
  length(intersect(a, b)) / u
}

# Per-stratum |estimate| filter at an arbitrary percentile, computed on the
# adjpvalue+dev_ratio-filtered network passed in.
apply_estimate_filter_pctl <- function(net, pctl) {
  stratum_key <- ifelse(net$link_type == 'TF_peak',
                        paste0('TF_', net$target_type, '_', net$target_measurement_type),
                        paste0('EP_', net$source_measurement_type, '_', net$target_measurement_type))
  thr_by_stratum <- tapply(abs(net$estimate), stratum_key,
                           function(x) unname(quantile(x, probs = pctl / 100)))
  thr_per_edge <- thr_by_stratum[stratum_key]
  net[!is.na(thr_per_edge) & abs(net$estimate) >= thr_per_edge, ]
}

# adjpvalue -> dev_ratio -> |estimate| percentile. Order is load-bearing.
filter_network_pctl <- function(net, adjp, dev_ratio, pctl) {
  net <- net[net$adjpvalue <= adjp, ]
  net <- net[net$dev_ratio >= dev_ratio, ]
  net <- apply_estimate_filter_pctl(net, pctl)
  net
}

# Core TF-set decoupling computation, for both node types.
run_decoupling <- function(net) {

  ## TF_enhancer: expressed enhancer, TF->aCRE vs TF->tCRE
  tf_enh      <- net[net$link_type == 'TF_peak' & net$target_type == 'Enhancer', ]
  tf_enh_aCRE <- tf_enh[tf_enh$target_measurement_type == 'aCRE', ]
  tf_enh_tCRE <- tf_enh[tf_enh$target_measurement_type == 'tCRE', ]
  enh_both    <- intersect(unique(tf_enh_aCRE$target), unique(tf_enh_tCRE$target))
  enh_both    <- enh_both[enh_both %in% expressed_enhancers]

  if (length(enh_both) > 0) {
    sets_a  <- split(tf_enh_aCRE$source, tf_enh_aCRE$target)[enh_both]
    sets_b  <- split(tf_enh_tCRE$source, tf_enh_tCRE$target)[enh_both]
    jac_enh <- mapply(jaccard, sets_a, sets_b)
  } else {
    jac_enh <- numeric(0)
  }

  ## TF_promoter: promoter, TF->aCRE vs TF->Expression
  tf_prom      <- net[net$link_type == 'TF_peak' & net$target_type == 'Promoter', ]
  tf_prom_aCRE <- tf_prom[tf_prom$target_measurement_type == 'aCRE', ]
  tf_prom_expr <- tf_prom[tf_prom$target_measurement_type == 'tCRE', ]
  prom_both    <- intersect(unique(tf_prom_aCRE$target), unique(tf_prom_expr$target))

  if (length(prom_both) > 0) {
    sets_a   <- split(tf_prom_aCRE$source, tf_prom_aCRE$target)[prom_both]
    sets_b   <- split(tf_prom_expr$source, tf_prom_expr$target)[prom_both]
    jac_prom <- mapply(jaccard, sets_a, sets_b)
  } else {
    jac_prom <- numeric(0)
  }

  list(
    n_nodes          = c(TF_enhancer = length(enh_both), TF_promoter = length(prom_both)),
    jaccard_by_theme = list(TF_enhancer = jac_enh, TF_promoter = jac_prom)
  )
}

run_grid_point <- function(sweep_name, sweep_value, adjp, dev_ratio, pctl) {
  message('  ', sweep_name, ' = ', sweep_value,
          '  (adjp<=', adjp, ', dev_ratio>=', dev_ratio, ', |estimate|>=p', pctl, ')')
  net_filt <- filter_network_pctl(net_full, adjp, dev_ratio, pctl)
  res <- run_decoupling(net_filt)
  res$sweep       <- sweep_name
  res$sweep_value <- sweep_value
  res
}

# Categorize a per-node Jaccard vector into 3 buckets and return %s.
categorize_jaccard <- function(jac) {
  jac <- jac[!is.na(jac)]
  n <- length(jac)
  if (n == 0) {
    return(data.frame(category = c('fully_uncoupled', 'intermediate', 'fully_coupled'),
                      pct = NA_real_, n = 0L))
  }
  data.frame(
    category = c('fully_uncoupled', 'intermediate', 'fully_coupled'),
    pct = c(mean(jac == 0), mean(jac > 0 & jac < 1), mean(jac == 1)) * 100,
    n   = c(sum(jac == 0), sum(jac > 0 & jac < 1), sum(jac == 1))
  )
}

#### Run sweeps (one-at-a-time, others held at baseline) ####

message('Sweeping adjpvalue_threshold...')
grid_adjp <- lapply(adjpvalue_grid, function(v)
  run_grid_point('adjpvalue_threshold', v, v, baseline_dev_ratio, baseline_pctl))

message('Sweeping dev_ratio_threshold...')
grid_devr <- lapply(dev_ratio_grid, function(v)
  run_grid_point('dev_ratio_threshold', v, baseline_adjpvalue, v, baseline_pctl))

message('Sweeping |estimate| stratum percentile...')
grid_pctl <- lapply(estimate_pctl_grid, function(v)
  run_grid_point('estimate_percentile', v, baseline_adjpvalue, baseline_dev_ratio, v))

#### Flatten into the category table ####

build_category_table <- function(grid_list) {
  do.call(rbind, lapply(grid_list, function(gp) {
    do.call(rbind, lapply(names(gp$jaccard_by_theme), function(th) {
      cat_df             <- categorize_jaccard(gp$jaccard_by_theme[[th]])
      cat_df$sweep       <- gp$sweep
      cat_df$sweep_value <- gp$sweep_value
      cat_df$theme       <- th
      cat_df$n_nodes     <- unname(gp$n_nodes[th])
      cat_df
    }))
  }))
}

category_table <- bind_rows(
  build_category_table(grid_adjp),
  build_category_table(grid_devr),
  build_category_table(grid_pctl)
) %>%
  relocate(sweep, sweep_value, theme, category, pct, n, n_nodes)

for (f in out_folders) {
  dir.create(f, showWarnings = FALSE, recursive = TRUE)
  write.csv(category_table, row.names = FALSE,
            file.path(f, 'sweep_category_table.csv'))
  message('Saved: ', file.path(f, 'sweep_category_table.csv'),
          ' (', nrow(category_table), ' rows)')
}

message('Done.')
