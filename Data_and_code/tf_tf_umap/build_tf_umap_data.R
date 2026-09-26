#### Regeneration script: TF-TF regulatory UMAP with dominance classification ####
#
# Feeds: ExtFig7/ext_fig7b_tf_umap_dominance.R
#
# For every TF in the network, builds a fixed 2D layout (UMAP on a combined
# TF-TF regulatory-association matrix) and classifies each TF by whether its
# outgoing regulatory influence on other TFs is carried mostly through direct
# edges (TF -> promoter(TF)) or indirect ones (TF -> enhancer -> promoter(TF)):
#   - Network loading + filtering: adjpvalue, dev_ratio, per-stratum
#     |estimate| p90, positive-only Enhancer_promoter edges.
#   - TF co-expression (Spearman on pseudotime-smoothed expression) and three
#     TF-to-TF regulatory weight matrices: direct-only, indirect-only, and
#     combined (the stronger of the two per TF pair).
#   - A combined association matrix (co-expression scaled by combined
#     regulatory weight) reduced via PCA then UMAP into a fixed 2D layout.
#   - Per-TF dominance classification from the fraction of outgoing edge
#     weight that is direct: >2/3 direct-dominant, <1/3 indirect-dominant,
#     otherwise mixed.
#
# Reads (see ../config.R to set these paths):
#   [primary_data_folder]/5_lasso_models/lasso_network_with_randomization.rds
#   [primary_data_folder]/4_preparing_dataset/dataset.rds
#   [raw_data_folder]/pseudotime.scRNA.tsv
#   [raw_data_folder]/tableA_human_stage_TFs_markers_updated.tsv
#     (stage-marker TF annotation table, used to restrict to the TFs called
#      as differentiation-stage markers)
#   ../network_loading/network_thresholds.txt
#
# Writes: ../../ExtFig7/data/tf_umap_dominance.rds
#   list with: layout_all (ggraph layout: name, x, y, pr, dominance, timepoint),
#              edge_direct, edge_indirect, edge_combined (data frames: from,
#              to, value; edge_combined takes, per TF pair, whichever of the
#              direct/indirect weight is stronger),
#              node_names, num_label_per_cat, dom_colours

rm(list = ls())
library(dplyr)
library(ggraph)
library(tidygraph)
library(igraph)
library(uwot)
library(reshape)

# ── CONTROL PANEL ─────────────────────────────────────────────────────────────
source('../config.R')  # defines primary_data_folder, raw_data_folder

network_file      <- file.path(primary_data_folder, '5_lasso_models/lasso_network_with_randomization.rds')
dataset_file      <- file.path(primary_data_folder, '4_preparing_dataset/dataset.rds')
pseudotime_file   <- file.path(raw_data_folder, 'pseudotime.scRNA.tsv')
selected_tfs_file <- file.path(raw_data_folder, 'tableA_human_stage_TFs_markers_updated.tsv')
thresholds_file   <- '../network_loading/network_thresholds.txt'
out_file          <- '../../ExtFig7/data/tf_umap_dominance.rds'

apply_estimate_threshold <- TRUE
sdev_min  <- 0.01
num_tf_pr <- 30
num_label_per_cat <- 5
# ──────────────────────────────────────────────────────────────────────────────

source('../network_loading/load_and_filter_network.R')
thresholds <- read_thresholds(thresholds_file)

message('Loading network...')
network <- readRDS(network_file)
network <- filter_network(network, thresholds, apply_estimate_threshold = apply_estimate_threshold)

# keep only positive Enhancer_promoter edges (activating enhancer-promoter
# links only; repressive links are excluded from the indirect regulatory path)
to_remove <- which(network$link_type == 'Enhancer_promoter' & network$estimate < 0)
if (length(to_remove) > 0) network <- network[-to_remove, ]
message('  Edges after filtering: ', nrow(network))

message('Loading dataset and pseudotime...')
dataset    <- readRDS(dataset_file)
pseudotime <- read.table(pseudotime_file, header = TRUE)
rownames(pseudotime) <- pseudotime$SEACell

selected_tfs <- read.table(selected_tfs_file, header = TRUE, sep = '\t')

#### Identify TFs ####

tfs <- unique(network$source[network$link_type == 'TF_peak'])
selected_tfs <- selected_tfs[selected_tfs$gene %in% tfs & selected_tfs$gene_type == 'TF', ]

message('  TFs in network: ', length(tfs))

#### TF co-expression matrix + pseudotime-of-peak-expression ####

message('Computing TF co-expression correlations...')

tf_rows <- paste0(tfs, '_Expression')
tf_rows <- tf_rows[tf_rows %in% rownames(dataset)]
dataset <- t(dataset[tf_rows, ])
colnames(dataset) <- sub('_Expression$', '', tf_rows)
tfs <- colnames(dataset)

pseudotime <- pseudotime[rownames(dataset), ]

loess_dataset <- dataset
for (i in seq_len(ncol(dataset))) {
  lo   <- loess(dataset[, i] ~ pseudotime$pseudotime, span = 0.3, degree = 0)
  expr <- predict(lo)
  loess_dataset[, i] <- expr / max(expr) * 100
}

tf_max_timepoint <- apply(loess_dataset, 2, which.max) / nrow(dataset)
tf_max_timepoint <- data.frame(tf = names(tf_max_timepoint), timepoint = tf_max_timepoint)

tf_corrs <- cor(dataset, method = 'spearman')

#### TF-to-TF regulatory weight matrices ####

message('Building TF-to-TF regulatory weight matrices...')

fill_coeff_matrix <- function(base_mat, edge_df) {
  mat <- base_mat
  for (k in seq_len(nrow(edge_df))) {
    s <- edge_df$source[k]
    t <- edge_df$target_gene[k]
    v <- edge_df$estimate[k]
    if (!s %in% rownames(mat) || !t %in% rownames(mat)) next
    for (idx in list(c(s, t), c(t, s))) {
      cur <- mat[idx[1], idx[2]]
      mat[idx[1], idx[2]] <- if (cur == 0 || abs(v) > abs(cur)) v else cur
    }
  }
  mat
}

init_mat <- matrix(0, length(tfs), length(tfs), dimnames = list(tfs, tfs))

direct_edges <- network %>%
  filter(link_type == 'TF_peak', target_type == 'Promoter',
         target_measurement_type == 'tCRE',
         source %in% tfs, target_gene %in% tfs,
         source != target_gene) %>%
  select(source, target_gene, estimate) %>%
  as.data.frame()

tf_coeffs <- fill_coeff_matrix(init_mat, direct_edges)

tf_to_enh <- network %>%
  filter(link_type == 'TF_peak', target_type == 'Enhancer', source %in% tfs) %>%
  select(source, target, estimate) %>%
  dplyr::rename(enhancer = target, coeff_tf_enh = estimate)

enh_to_tf_prom <- network %>%
  filter(link_type == 'Enhancer_promoter', target_type == 'Promoter',
         target_measurement_type == 'tCRE', target_gene %in% tfs) %>%
  select(source, target_gene, estimate) %>%
  dplyr::rename(enhancer = source, coeff_enh_prom = estimate)

indirect_edges <- merge(tf_to_enh, enh_to_tf_prom, by = 'enhancer') %>%
  filter(source != target_gene) %>%
  mutate(estimate = ifelse(abs(coeff_tf_enh) < abs(coeff_enh_prom),
                           coeff_tf_enh, coeff_enh_prom)) %>%
  group_by(source, target_gene) %>%
  summarise(estimate = estimate[which.max(abs(estimate))], .groups = 'drop') %>%
  as.data.frame()

tf_coeffs_enh <- fill_coeff_matrix(init_mat, indirect_edges)

tf_coeffs_all <- tf_coeffs
stronger <- abs(tf_coeffs_enh) > abs(tf_coeffs)
tf_coeffs_all[stronger] <- tf_coeffs_enh[stronger]

message('  Non-zero entries — direct: ', sum(tf_coeffs != 0),
        '  indirect: ',  sum(tf_coeffs_enh != 0),
        '  combined: ',  sum(tf_coeffs_all != 0))

# combined association matrix: co-expression correlation scaled up by
# regulatory weight, so TF pairs that are both co-expressed and strongly
# linked in the network sit closer together in the UMAP embedding
tf_assocs_all <- tf_corrs * sqrt(abs(tf_coeffs_all) + 1)

#### Section 3 — fixed-layout UMAP + outgoing-edge dominance ####

message('Building fixed-layout comparison UMAP...')

pc_all <- irlba::prcomp_irlba(tf_assocs_all, n = 50)
to_keep_all <- max(which(pc_all$sdev^2 / sum(pc_all$sdev^2) > sdev_min)) + 1
to_keep_all <- min(to_keep_all, ncol(pc_all$x))
x_all <- pc_all$x[, seq_len(to_keep_all)]
rownames(x_all) <- rownames(tf_assocs_all)
message('  PCs retained: ', to_keep_all)

set.seed(12345)
umap_fixed <- uwot::umap(x_all)
rownames(umap_fixed) <- rownames(x_all)

diag(tf_coeffs)     <- 0
diag(tf_coeffs_enh) <- 0
sum_direct   <- rowSums(abs(tf_coeffs))
sum_indirect <- rowSums(abs(tf_coeffs_enh))
total        <- sum_direct + sum_indirect
frac_direct  <- ifelse(total > 0, sum_direct / total, NA_real_)

dominance <- setNames(dplyr::case_when(
  is.na(frac_direct) ~ 'No outgoing edges',
  frac_direct > 2/3   ~ 'Direct-dominant',
  frac_direct < 1/3   ~ 'Indirect-dominant',
  TRUE                ~ 'Mixed'
), names(frac_direct))

dom_colours <- c(
  'Direct-dominant'   = '#1f78b4',
  'Indirect-dominant' = '#e31a1c',
  'Mixed'             = '#b2df8a',
  'No outgoing edges' = 'grey80'
)

coeff_plot <- tf_coeffs_all
diag(coeff_plot) <- 0
coeff_plot[upper.tri(coeff_plot)] <- 0
suppressWarnings(edge_all <- reshape::melt(coeff_plot))
colnames(edge_all) <- c('from', 'to', 'value')
edge_all <- edge_all[edge_all$value != 0, ]
edge_all$value <- abs(edge_all$value)

graph_all    <- igraph::graph_from_data_frame(edge_all)
nodes_all    <- names(V(graph_all))
umap_sub_all <- umap_fixed[nodes_all, ]
layout_all   <- ggraph::create_layout(graph_all, layout = umap_sub_all)

pr_all     <- igraph::page_rank(igraph::reverse_edges(graph_all))
node_names <- as.character(layout_all$name)
layout_all$pr        <- pr_all$vector[node_names]
layout_all$dominance <- dominance[node_names]
layout_all$timepoint <- tf_max_timepoint[match(node_names, tf_max_timepoint$tf), 'timepoint']

message('  Dominance breakdown: ', paste(names(table(layout_all$dominance)),
        table(layout_all$dominance), sep = '=', collapse = ', '))

make_edge_table <- function(mat) {
  m <- mat; diag(m) <- 0; m[upper.tri(m)] <- 0
  suppressWarnings(df <- reshape::melt(m))
  colnames(df) <- c('from', 'to', 'value')
  df <- df[df$value != 0, ]
  df$value <- abs(df$value)
  df
}
edge_direct   <- make_edge_table(tf_coeffs)
edge_indirect <- make_edge_table(tf_coeffs_enh)
# Same combined matrix already used to build layout_all above (tf_coeffs_all:
# per TF pair, whichever of the direct/indirect weight is stronger), just
# run through the same make_edge_table() helper as the other two tables for
# a consistent from/to/value structure.
edge_combined <- make_edge_table(tf_coeffs_all)

#### Save ####

dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(list(
  layout_all         = layout_all,
  edge_direct        = edge_direct,
  edge_indirect      = edge_indirect,
  edge_combined      = edge_combined,
  node_names         = node_names,
  umap_fixed         = umap_fixed,
  num_label_per_cat  = num_label_per_cat,
  dom_colours        = dom_colours
), out_file)

message('Saved: ', out_file)
