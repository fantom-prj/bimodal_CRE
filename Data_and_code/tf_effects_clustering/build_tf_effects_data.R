#### Regeneration script: TF cumulative-effect matrix + heatmap clustering ####
#
# Produces the shared computation behind TWO panels:
#   Fig5b  ComplexHeatmap of scaled TF effects across 4 target x modality columns
#   Fig5c  TF characteristics dot plot (t-tests of TF effects vs known TF roles)
#
# For every TF, sums its regulatory effect (signed model coefficient) across
# all edges into each of several target classes (all TF-peak edges, plus,
# separately, promoter and mediating-enhancer targets split by aCRE/tCRE
# measurement type), producing a cumulative-effect profile per TF. That
# profile is then:
#   - assembled into a 4-column matrix (mediating-enhancer aCRE/tCRE,
#     promoter aCRE/tCRE), column-scaled, and hierarchically clustered
#     (Ward's method on a sign-preserving sqrt-compressed transform) into
#     activator-like vs repressor-like groups, with activator-like TFs
#     further split into sub-clusters;
#   - ranked overall via a log rank-product across the 4 columns.
#   (The heatmap draw() call and the characteristics dot plot themselves live
#    in the panel scripts, Fig5/fig5b_heatmap.R and
#    Fig5/fig5c_tf_characteristics.R.)
#
# Control-panel values fixed to the variant used for the published figures:
#   use_abs_estimate         = FALSE  (signed cumulative effect, so activating
#                                      and repressive influence can be told
#                                      apart)
#   apply_estimate_threshold = TRUE   (per-stratum |estimate| p90 filter)
#   num_tf_to_show           = 25     (kept for provenance; not used by
#                                      either downstream panel directly)
#   target_k                 = 5      (Activator-like sub-clusters)
#
# Clustering notes (do not "simplify"):
#   * sqrt_signed(x) = sqrt(|x|) * sign(x) compresses outliers while keeping
#     sign and ordering; Euclidean distance + Ward.D2 is computed on it, while
#     the heatmap DISPLAYS the plain column-scaled matrix sh.
#   * Rows are ordered by descending log rank-product rk before splitting.
#   * row_group_fine's factor levels set the heatmap's top-to-bottom group
#     order; row_order is NOT compatible with row_split + cluster_rows = TRUE.
#   * The column dendrogram has its last two leaves swapped via
#     dendextend::rotate() so the promoter columns read aCRE-then-tCRE.
#
# Outputs:
#   Fig5/data/tf_effects_clustering.rds
#     a list with: h_wide, sh, sh_sqrt, rk, row_group, row_group_fine,
#     col_order, cdend, results_bar
#   Fig5/data/TF_grouping_by_chromvar_correlation_and_literature.tsv
#     verbatim copy of the pipeline's TF role annotation table (small enough
#     to ship with the figure); read by Fig5c.

rm(list = ls())
library(tidyverse)
library(dendextend)

# ── CONTROL PANEL ─────────────────────────────────────────────────────────────
source('../config.R')  # defines primary_data_folder, raw_data_folder
network_file <- file.path(primary_data_folder,
                          '5_lasso_models/lasso_network_with_randomization.rds')

tf_characteristics_file <- file.path(raw_data_folder,
                                     'TF_grouping_by_chromvar_correlation_and_literature.tsv')

thresholds_file <- '../network_loading/network_thresholds.txt'
out_folder      <- '../../Fig5/data'

use_abs_estimate         <- FALSE
apply_estimate_threshold <- TRUE
num_tf_to_show           <- 25
target_k                 <- 5

target_types <- c('All', 'Promoter', 'Enhancer', 'Expressed Enhancer', 'Mediating Enhancer')
# ──────────────────────────────────────────────────────────────────────────────

source('../network_loading/load_and_filter_network.R')

thresholds <- read_thresholds(thresholds_file)
message('Thresholds loaded from: ', thresholds_file)

#### Load and filter network ####

message('Loading network...')
n <- readRDS(network_file)
n <- filter_network(n, thresholds, apply_estimate_threshold = apply_estimate_threshold)

n$abs_estimate   <- abs(n$estimate)
n$signed_deltaR2 <- n$deltaR2 * sign(n$estimate)

message('  Edges after filtering: ', nrow(n))

#### Helper: log-signed transform ####

log_signed <- function(x) log2(abs(x) + 1) * sign(x)

#### Part 1 — per-target-type cumulative-effect summaries (results_bar) ####

message('Building per-target-type effect summaries...')

results_bar <- list()

for (tt in target_types) {

  idx <- switch(tt,
    'All'                = n$link_type == 'TF_peak',
    'Expressed Enhancer' = n$link_type == 'TF_peak' &
                           n$target_type == 'Enhancer' &
                           n$target_measurement_type == 'tCRE',
    'Mediating Enhancer' = n$link_type == 'TF_peak' &
                           n$target_type == 'Enhancer' &
                           n$target %in% n$source,
                           n$link_type == 'TF_peak' & n$target_type == tt
  )

  sub_n <- n[idx, c('source', 'estimate', 'target_measurement_type')]
  colnames(sub_n)[2] <- 'statistic'

  if (use_abs_estimate) sub_n$statistic <- abs(sub_n$statistic)

  to_plot <- sub_n %>%
    group_by(source, target_measurement_type) %>%
    summarise(
      total_effect       = sum(statistic),
      average_effect     = mean(statistic),
      log_total_effect   = log_signed(sum(statistic)),
      log_average_effect = log_signed(mean(statistic)),
      .groups = 'drop'
    ) %>%
    as.data.frame()

  results_bar[[tt]] <- to_plot
  message('  ', tt, ': ', length(unique(to_plot$source)), ' TFs')
}

#### Part 3 — heatmap matrix assembly and clustering ####

message('Assembling heatmap matrix and clustering...')

make_wide <- function(tt, label) {
  results_bar[[tt]] %>%
    select(source, target_measurement_type, log_total_effect) %>%
    mutate(col = paste0(label, '_', target_measurement_type)) %>%
    select(source, col, log_total_effect) %>%
    pivot_wider(names_from = col, values_from = log_total_effect)
}

h_wide <- full_join(
  make_wide('Mediating Enhancer', 'mediating_enhancer'),
  make_wide('Promoter',           'promoter'),
  by = 'source'
) %>%
  column_to_rownames('source') %>%
  as.matrix()

h_wide[is.na(h_wide)] <- 0

col_order <- c('mediating_enhancer_aCRE', 'mediating_enhancer_tCRE',
               'promoter_aCRE',           'promoter_tCRE')
h_wide <- h_wide[, col_order]

# scale columns (no centring)
sh <- scale(h_wide, center = FALSE, scale = TRUE)

sqrt_signed <- function(x) sqrt(abs(x)) * sign(x)
sh_sqrt <- sqrt_signed(sh)

# rank-product ordering (display ordering only)
rk <- log10(apply(apply(sh_sqrt, 2, rank), 1, prod))
sh_sqrt <- sh_sqrt[order(rk, decreasing = TRUE), ]
sh      <- sh     [order(rk, decreasing = TRUE), ]
rk      <- rk     [order(rk, decreasing = TRUE)]

# row split: sign of the column with the largest |sqrt-transformed| value
dominant_col  <- apply(sh_sqrt, 1, function(r) which.max(abs(r)))
dominant_sign <- sapply(seq_len(nrow(sh_sqrt)),
                        function(i) sign(sh_sqrt[i, dominant_col[i]]))
row_group <- ifelse(dominant_sign >= 0, 'Activator-like', 'Repressor-like')
names(row_group) <- rownames(sh_sqrt)

# Activator-like sub-clustering
act_idx <- which(row_group == 'Activator-like')
act_hc  <- hclust(dist(sh_sqrt[act_idx, ]), method = 'ward.D2')
act_sub <- cutree(act_hc, k = target_k)

message('  Activator-like clusters: ', max(act_sub),
        ' (sizes: ', paste(sort(table(act_sub), decreasing = TRUE), collapse = ', '), ')')

# remap raw cutree cluster numbers to the published labels
cluster_label_map <- c('1' = 'Activator-like 1',
                       '3' = 'Activator-like 2',
                       '2' = 'Mixed 1',
                       '4' = 'Mixed 2',
                       '5' = 'Mixed 3')
act_labels <- cluster_label_map[as.character(act_sub)]
names(act_labels) <- names(act_sub)

row_group_fine <- row_group
row_group_fine[names(act_sub)] <- act_labels

# factor levels control heatmap top-to-bottom group order
group_level_order <- c('Activator-like 1', 'Activator-like 2',
                       'Mixed 1', 'Mixed 2', 'Repressor-like', 'Mixed 3')
row_group_fine <- factor(row_group_fine, levels = group_level_order)

# column dendrogram — Euclidean on sqrt-signed values, last two leaves swapped
cdend <- as.dendrogram(hclust(dist(t(sh_sqrt)), method = 'ward.D2'))
col_leaf_order <- labels(cdend)
swapped_order  <- c(col_leaf_order[1:2], col_leaf_order[4], col_leaf_order[3])
cdend <- dendextend::rotate(cdend, swapped_order)

#### Save ####

dir.create(out_folder, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(
  h_wide         = h_wide,
  sh             = sh,
  sh_sqrt        = sh_sqrt,
  rk             = rk,
  row_group      = row_group,
  row_group_fine = row_group_fine,
  col_order      = col_order,
  cdend          = cdend,
  results_bar    = results_bar,
  num_tf_to_show = num_tf_to_show,
  stat_label     = if (use_abs_estimate) 'abs' else 'signed'
), file.path(out_folder, 'tf_effects_clustering.rds'))

message('Saved: ', file.path(out_folder, 'tf_effects_clustering.rds'),
        '  (', nrow(sh), ' TFs x ', ncol(sh), ' columns)')

# Ship the TF role annotation table alongside the RDS (Fig5c reads it).
file.copy(tf_characteristics_file,
          file.path(out_folder, 'TF_grouping_by_chromvar_correlation_and_literature.tsv'),
          overwrite = TRUE)
message('Saved: ', file.path(out_folder,
        'TF_grouping_by_chromvar_correlation_and_literature.tsv'))

message('Done.')
