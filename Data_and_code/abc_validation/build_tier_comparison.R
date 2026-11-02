#### Regeneration script: network-tier vs ABC concordance summary ####
#
# Feeds (via copy_tier_comparison.R -> Fig5/data/): ExtFig5d.
#
# Validates the GRN's Enhancer_promoter edges against ABC (Activity-by-
# Contact) enhancer-gene scores, computed independently from chromatin
# accessibility, H3K27ac and Hi-C (see ../Fig4.ABC/). Three network tiers
# are compared:
#   base         every proximity/Hi-C-eligible candidate pair in the base
#                network, before any model fitting (the no-selection floor)
#   elastic_net  edges kept by the elastic-net model, before the standard
#                significance/effect-size filters
#   filtered     the final high-confidence GRN
# against four ABC variants: each cell type (iPSC, NSC, Neuron) and the
# per-pair maximum across cell types.
#
# Eligibility. ABC's full score matrix has a row for every (enhancer, gene)
# pair ABC evaluated, with the actual score (including 0), and no row for
# pairs it never evaluated; a per-cell-type NA means the pair was not
# evaluated in that cell type. For each variant the comparison universe is
# the base-network candidate pairs that ABC evaluated for that variant, so
# "never evaluated" is never counted as an ABC negative. ABC peaks that are
# promoter-like in the aCRE catalogue (Reg_Ele CGI / TATA_box / Both) are
# removed first, so that only genuine enhancer-gene pairs are compared.
# A pair is ABC-positive if its score exceeds abc_score_threshold (0.02,
# ABC's own reporting threshold).
#
# Metrics per (tier, variant): Jaccard of network vs ABC-positive pairs;
# precision, recall, F1 and Fisher's exact test within the eligible
# universe; ABC base rate; Spearman/Pearson correlation between |estimate|
# and ABC score over network pairs in the universe (tiers with estimates
# only); AUC of the ABC score for discriminating network edges. The base
# tier calls every pair in its universe, so its Fisher test is undefined
# (NA) and its AUC is NA.
#
# Reads (see ../config.R to set these paths):
#   [raw_data_folder]/ABCscore_matrix_3cells.tsv.gz
#     full ABC score matrix: peakID, TargetGene, iPSC, NSC, Neuron
#     (decompressed with gzip, so no extra R package is needed for .gz;
#      an uncompressed .tsv also works)
#   [raw_data_folder]/all_aCRE.PE.final.tsv  (Reg_Ele column)
#   [primary_data_folder]/3_preparing_base_network/base_network.rds
#   [primary_data_folder]/5_lasso_models/lasso_network_with_randomization.rds
#   [primary_data_folder]/5_lasso_models/lasso_network_filtered.rds
#
# Writes [primary_data_folder]/19_abc_validation_full_matrix/tier_comparison.csv
#   12 rows (3 tiers x 4 ABC variants). Columns: tier, abc_variant,
#   n_network_edges, n_abc_edges, n_intersect, n_union, jaccard,
#   n_eligible_universe, abc_base_rate, precision, recall, f1, fisher_OR,
#   fisher_ci_low, fisher_ci_high, fisher_p, n_both_for_corr, spearman_rho,
#   spearman_p, pearson_r, pearson_p, auc
#
# Memory: the matrix has ~24 million rows (~73 million in long format);
# expect a peak of tens of GB of RAM.

rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(pROC)
})

# ── CONTROL PANEL ─────────────────────────────────────────────────────────────
source('../config.R')  # defines primary_data_folder, raw_data_folder

abc_matrix_file     <- file.path(raw_data_folder, 'ABCscore_matrix_3cells.tsv.gz')
acre_catalogue_file <- file.path(raw_data_folder, 'all_aCRE.PE.final.tsv')
base_network_file   <- file.path(primary_data_folder, '3_preparing_base_network/base_network.rds')
network_files       <- c(
  elastic_net = file.path(primary_data_folder, '5_lasso_models/lasso_network_with_randomization.rds'),
  filtered    = file.path(primary_data_folder, '5_lasso_models/lasso_network_filtered.rds')
)
out_folder <- file.path(primary_data_folder, '19_abc_validation_full_matrix')

abc_score_threshold <- 0.02
cell_types          <- c('iPSC', 'NSC', 'Neuron')
tier_order          <- c('base', 'elastic_net', 'filtered')
# ──────────────────────────────────────────────────────────────────────────────

dir.create(out_folder, showWarnings = FALSE, recursive = TRUE)

#### ABC score matrix ####

message('Loading ABC score matrix (large)...')
abc_mat <- if (grepl('\\.gz$', abc_matrix_file)) {
  fread(cmd = paste('gzip -dc', shQuote(abc_matrix_file)), data.table = FALSE)
} else {
  fread(abc_matrix_file, data.table = FALSE)
}
message('  ', nrow(abc_mat), ' rows, ', length(unique(abc_mat$peakID)), ' enhancer peaks, ',
        length(unique(abc_mat$TargetGene)), ' target genes')

# drop ABC "enhancers" that are promoter-like in the aCRE catalogue
acre_raw <- fread(acre_catalogue_file, sep = '\t', quote = '', data.table = FALSE,
                  stringsAsFactors = FALSE)
promoter_like_peaks <- unique(acre_raw$peakID[acre_raw$Reg_Ele %in% c('CGI', 'TATA_box', 'Both')])
rm(acre_raw)

n_abc_peaks_before <- length(unique(abc_mat$peakID))
abc_mat <- abc_mat[!(abc_mat$peakID %in% promoter_like_peaks), ]
message('  Removed ', n_abc_peaks_before - length(unique(abc_mat$peakID)),
        ' promoter-like ABC peaks')

# long format: one row per (enhancer, gene, cell_type)
abc_long <- abc_mat %>%
  rename(enhancer = peakID, gene = TargetGene) %>%
  pivot_longer(cols = all_of(cell_types), names_to = 'cell_type', values_to = 'ABC.Score') %>%
  filter(!is.na(gene), gene != '')
rm(abc_mat); invisible(gc())

# per cell type: pairs evaluated (non-NA score) in that cell type;
# max_across_celltypes: pairs evaluated in >= 1 cell type, max score over them
abc_by_celltype <- setNames(
  lapply(cell_types, function(ct) abc_long %>% filter(cell_type == ct, !is.na(ABC.Score))),
  cell_types
)
abc_max <- abc_long %>%
  filter(!is.na(ABC.Score)) %>%
  group_by(enhancer, gene) %>%
  summarise(ABC.Score = max(ABC.Score), .groups = 'drop')
abc_variants <- c(abc_by_celltype, list(max_across_celltypes = abc_max))
rm(abc_long); invisible(gc())

abc_eligible_key_by_variant <- lapply(abc_variants, function(tab) {
  unique(paste(tab$enhancer, tab$gene, sep = '::'))
})

#### Base-network candidate pairs ####

message('Loading base network candidate pairs...')
base_net <- readRDS(base_network_file)
base_ep  <- base_net[base_net$link_type == 'Enhancer_promoter', ]
base_ep  <- base_ep[!is.na(base_ep$target_gene) & base_ep$target_gene != '', ]
candidate_df <- base_ep %>%
  transmute(enhancer = source, gene = target_gene) %>%
  group_by(enhancer, gene) %>%
  summarise(.groups = 'drop') %>%
  mutate(key = paste(enhancer, gene, sep = '::'))
rm(base_net, base_ep); invisible(gc())
message('  ', nrow(candidate_df), ' candidate (enhancer, gene) pairs')

#### Comparison ####

extract_network_ep_pairs <- function(net) {
  ep <- net[net$link_type == 'Enhancer_promoter', ]
  ep <- ep[!is.na(ep$target_gene) & ep$target_gene != '', ]
  ep %>%
    transmute(enhancer = source, gene = target_gene, abs_estimate = abs(estimate)) %>%
    group_by(enhancer, gene) %>%
    summarise(abs_estimate = max(abs_estimate), .groups = 'drop')
}

compare_one <- function(net_pairs, abc_tab, abc_variant_name, candidate_keys) {
  universe <- intersect(candidate_keys, abc_eligible_key_by_variant[[abc_variant_name]])

  net_key     <- paste(net_pairs$enhancer, net_pairs$gene, sep = '::')
  abc_pos_key <- abc_tab %>% filter(ABC.Score > abc_score_threshold) %>%
    mutate(key = paste(enhancer, gene, sep = '::')) %>% pull(key)

  n_inter <- length(intersect(net_key, abc_pos_key))
  n_union <- length(union(net_key, abc_pos_key))

  in_net <- universe %in% net_key
  in_abc <- universe %in% abc_pos_key

  fisher_res <- list(odds_ratio = NA, ci_low = NA, ci_high = NA, p_value = NA)
  if (length(unique(in_net)) == 2) {
    fisher_res <- tryCatch({
      ft <- fisher.test(table(factor(in_net, c(FALSE, TRUE)), factor(in_abc, c(FALSE, TRUE))))
      list(odds_ratio = unname(ft$estimate), ci_low = ft$conf.int[1],
           ci_high = ft$conf.int[2], p_value = ft$p.value)
    }, error = function(e) fisher_res)
  }

  tp <- sum(in_net & in_abc); fp <- sum(in_net & !in_abc); fn <- sum(!in_net & in_abc)
  precision <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  recall    <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  f1        <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0)
                 2 * precision * recall / (precision + recall) else NA_real_

  abc_score_lookup <- abc_tab %>% mutate(key = paste(enhancer, gene, sep = '::')) %>%
    select(key, ABC.Score) %>% distinct(key, .keep_all = TRUE)

  # |estimate| vs ABC score, over network pairs inside the universe
  n_both <- 0
  spear <- pear <- NULL
  if ('abs_estimate' %in% colnames(net_pairs)) {
    cont_df <- data.frame(key = intersect(net_key, universe)) %>%
      left_join(net_pairs %>% mutate(key = paste(enhancer, gene, sep = '::')) %>%
                  select(key, abs_estimate), by = 'key') %>%
      left_join(abc_score_lookup, by = 'key') %>%
      mutate(ABC.Score = ifelse(is.na(ABC.Score), 0, ABC.Score))
    n_both <- nrow(cont_df)
    if (n_both >= 3) {
      spear <- cor.test(cont_df$abs_estimate, cont_df$ABC.Score, method = 'spearman')
      pear  <- cor.test(cont_df$abs_estimate, cont_df$ABC.Score, method = 'pearson')
    }
  }

  # AUC of the ABC score for network-edge membership within the universe
  auc_df <- data.frame(key = universe, is_network_edge = in_net) %>%
    left_join(abc_score_lookup, by = 'key') %>%
    mutate(ABC.Score = ifelse(is.na(ABC.Score), 0, ABC.Score))
  auc_val <- NA_real_
  if (length(unique(auc_df$is_network_edge)) == 2) {
    auc_val <- as.numeric(pROC::auc(pROC::roc(auc_df$is_network_edge, auc_df$ABC.Score,
                                              quiet = TRUE)))
  }

  data.frame(
    n_network_edges = length(net_key), n_abc_edges = length(abc_pos_key),
    n_intersect = n_inter, n_union = n_union,
    jaccard = if (n_union == 0) NA_real_ else n_inter / n_union,
    n_eligible_universe = length(universe), abc_base_rate = mean(in_abc),
    precision = precision, recall = recall, f1 = f1,
    fisher_OR = fisher_res$odds_ratio, fisher_ci_low = fisher_res$ci_low,
    fisher_ci_high = fisher_res$ci_high, fisher_p = fisher_res$p_value,
    n_both_for_corr = n_both,
    spearman_rho = if (!is.null(spear)) unname(spear$estimate) else NA_real_,
    spearman_p   = if (!is.null(spear)) spear$p.value else NA_real_,
    pearson_r    = if (!is.null(pear)) unname(pear$estimate) else NA_real_,
    pearson_p    = if (!is.null(pear)) pear$p.value else NA_real_,
    auc          = auc_val
  )
}

net_pairs_by_tier <- list(
  base        = candidate_df %>% select(enhancer, gene),
  elastic_net = extract_network_ep_pairs(readRDS(network_files[['elastic_net']])),
  filtered    = extract_network_ep_pairs(readRDS(network_files[['filtered']]))
)

rows <- list()
for (tier in tier_order) {
  for (abc_variant in names(abc_variants)) {
    message('Tier ', tier, ' x ABC variant ', abc_variant)
    rows[[paste(tier, abc_variant)]] <- data.frame(
      tier = tier, abc_variant = abc_variant,
      compare_one(net_pairs_by_tier[[tier]], abc_variants[[abc_variant]],
                  abc_variant, candidate_df$key)
    )
  }
}
tier_comparison <- bind_rows(rows)

out_file <- file.path(out_folder, 'tier_comparison.csv')
write.csv(tier_comparison, row.names = FALSE, out_file)

message('Saved: ', out_file, ' (', nrow(tier_comparison), ' rows)')
f1_max <- tier_comparison$f1[tier_comparison$abc_variant == 'max_across_celltypes']
message('  F1 (max_across_celltypes), ', paste(tier_order, collapse = ' -> '), ': ',
        paste(round(f1_max, 3), collapse = ' -> '))
message('Done.')
