#### Regeneration script: enhancer-promoter epigenomic feature enrichment ####
#
# Feeds: ExtFig7/ext_fig7a_ep_enrichment_collapsed.R
#
# Tests whether pairs of epigenomic features (chromatin state, stage-specific
# super-enhancer membership, CpG islands, TATA boxes) co-occur more or less
# often than expected between an enhancer and the promoter(s) it regulates in
# the GRN. For each of 9x9 enhancer/promoter feature combinations, builds a
# 2x2 (or larger) contingency table over enhancer-promoter edges and tests
# for association (chi-squared test), reporting odds ratio and FDR-adjusted
# p-value. A combination where either feature takes a single value across the
# promoters cannot be tested and is reported with NA odds ratio and p-value.
# Repeats this across 3 strata of enhancer-promoter edges: all edges, edges
# from expressed (tCRE) enhancers only, and edges from "mediating"
# enhancers only (enhancers that are themselves under TF control, i.e. also a
# target of a TF-to-peak edge in the network) -- mediating status is resolved
# directly from the network's edge structure.
#
# Chromatin states come from the per-cell-type state calls (Active, Primed,
# Repressed, Bivalent, Null; one state per peak and cell type). Each state is
# used as a feature that is 'Yes' when the peak has that state in at least one
# of iPSC, NSC and Neuron; Null is the implicit remainder and is not tested.
#
# Reads (see ../config.R to set these paths):
#   [primary_data_folder]/5_lasso_models/lasso_network_with_randomization.rds
#   [primary_data_folder]/4_preparing_dataset/dataset.rds  (not used by this
#     collapsed-panel computation once the mediating-enhancer set is resolved
#     from the network alone, but loaded for consistency with the other
#     enhancer-role analyses run against the same inputs)
#   [raw_data_folder]/all_aCRE.TSS.access.markalone_cuttag.tsv.gz
#     (peak feature table: chromatin states per cell type, SE, CGI, TATA box)
#   ../network_loading/network_thresholds.txt
#
# Writes: ../../ExtFig7/data/ep_enrichment_collapsed.csv
#   columns: enhancer, promoter, stratum, odds_ratio, pvalue, adjpvalue
#   (enhancer/promoter columns here are FEATURE names, e.g. 'any.Active_enhancer',
#    not peak ids -- one row per feature x feature x stratum combination)

rm(list = ls())
library(dplyr)
library(stringr)

# ── CONTROL PANEL ─────────────────────────────────────────────────────────────
source('../config.R')  # defines primary_data_folder, raw_data_folder

network_file              <- file.path(primary_data_folder, '5_lasso_models/lasso_network_with_randomization.rds')
dataset_file              <- file.path(primary_data_folder, '4_preparing_dataset/dataset.rds')
peak_characteristics_file <- file.path(raw_data_folder, 'all_aCRE.TSS.access.markalone_cuttag.tsv.gz')
thresholds_file           <- '../network_loading/network_thresholds.txt'
out_dir                   <- '../../ExtFig7/data'

apply_estimate_threshold <- TRUE
significance_threshold   <- 0.05

chromatin_states <- c('Active', 'Primed', 'Repressed', 'Bivalent')
characteristics_enhancer <- c(paste0('any.', chromatin_states),
                               'iPSC_SE', 'NSC_SE', 'Neuron_SE',
                               'CGI', 'TATA_box')
characteristics_promoter <- characteristics_enhancer
# ──────────────────────────────────────────────────────────────────────────────

source('../network_loading/load_and_filter_network.R')
thresholds <- read_thresholds(thresholds_file)
adjpvalue_threshold <- thresholds$adjpvalue_threshold
dev_ratio_threshold <- thresholds$dev_ratio_threshold
abs_estimate_thresholds_ep <- thresholds$abs_estimate_thresholds[
  c('EP_aCRE_aCRE', 'EP_aCRE_tCRE', 'EP_tCRE_aCRE', 'EP_tCRE_tCRE')]
names(abs_estimate_thresholds_ep) <- c('aCRE_aCRE', 'aCRE_tCRE', 'tCRE_aCRE', 'tCRE_tCRE')

message('Loading network...')
n <- readRDS(network_file)
n <- n[n$adjpvalue <= adjpvalue_threshold, ]
n <- n[n$dev_ratio >= dev_ratio_threshold, ]
message('  Edges after adjpvalue + dev_ratio filter: ', nrow(n))

message('Loading dataset...')
dataset <- readRDS(dataset_file)

#### Part 1 (partial) — enhancer role classification: 'mediating' only ####

message('Classifying enhancers (mediating set only)...')

tmp <- unique(n$target[n$target_type == 'Enhancer'])
mediating_enhancers <- tmp[tmp %in% n$source]
message('  Mediating enhancers: ', length(mediating_enhancers))

#### Part 3 (collapsed only) — EP enrichment tests ####

message('Loading peak feature table...')
# tab-separated, with spaces inside some fields
peak_data <- read.table(peak_characteristics_file, header = TRUE, sep = '\t',
                        quote = '', comment.char = '', stringsAsFactors = FALSE)
peak_data$tCRE   <- str_to_title(peak_data$tCRE)
peak_data$peakID <- gsub('-', '_', peak_data$peakID)
colnames(peak_data) <- sub('^(iPSC|NSC|Neuron)\\.SE$', '\\1_SE', colnames(peak_data))
# a peak without a state call in a cell type counts as not having the state there
for (state in chromatin_states) {
  peak_data[[paste0('any.', state)]] <- ifelse(
    peak_data[[paste0('iPSC.', state)]] %in% 'Yes' |
      peak_data[[paste0('NSC.', state)]] %in% 'Yes' |
      peak_data[[paste0('Neuron.', state)]] %in% 'Yes',
    'Yes', 'No'
  )
}
stopifnot(all(characteristics_enhancer %in% colnames(peak_data)))

ep <- n[n$link_type == 'Enhancer_promoter', ]
ep$edge_id <- paste0(ep$source, '_', ep$target_gene)

if (apply_estimate_threshold) {
  stratum_key <- paste0(ep$source_measurement_type, '_', ep$target_measurement_type)
  abs_thr     <- abs_estimate_thresholds_ep[stratum_key]
  ep <- ep[!is.na(abs_thr) & abs(ep$estimate) >= abs_thr, ]
}
message('  EP edges after all filters: ', nrow(ep))

measurement_strata <- list(
  all       = list(source = c('aCRE', 'tCRE'), target = c('aCRE', 'tCRE'), mediating_only = FALSE),
  tCRE_all  = list(source = 'tCRE',            target = c('aCRE', 'tCRE'), mediating_only = FALSE),
  mediating = list(source = c('aCRE', 'tCRE'), target = c('aCRE', 'tCRE'), mediating_only = TRUE)
)

compute_or_p <- function(lk, raw_enh, raw_prom) {
  # a feature with a single level across the promoters (e.g. no promoter whose
  # enhancers all carry a rare state) leaves nothing to test: report NA
  if (length(unique(raw_enh)) < 2 || length(unique(raw_prom)) < 2) {
    return(list(or = NA_real_, pv = NA_real_))
  }
  or <- (lk$weight[4] / lk$weight[3]) / (lk$weight[2] / lk$weight[1])
  pv <- tryCatch(chisq.test(raw_enh, raw_prom)$p.value, error = function(e) 1)
  list(or = or, pv = pv)
}

res_cols      <- c('enhancer', 'promoter', 'stratum', 'odds_ratio', 'pvalue', 'adjpvalue')
res_collapsed <- data.frame(matrix(NA, 0, length(res_cols)), stringsAsFactors = FALSE)
colnames(res_collapsed) <- res_cols

message('Running EP enrichment tests (collapsed only)...')

for (k in seq_along(measurement_strata)) {

  mt_name <- names(measurement_strata)[k]
  idx     <- ep$source_measurement_type %in% measurement_strata[[k]]$source &
             ep$target_measurement_type %in% measurement_strata[[k]]$target
  cur <- ep[idx, ]
  if (mt_name == 'tCRE_all') cur <- cur[!(cur$edge_id %in% ep$edge_id[!idx]), ]
  if (measurement_strata[[k]]$mediating_only)
    cur <- cur[cur$source %in% mediating_enhancers, ]

  cur_base <- cur %>% dplyr::select(source, target, source_type, target_type) %>% unique()

  for (feat_enh in characteristics_enhancer) {
    for (feat_prom in characteristics_promoter) {

      to_plot <- cur_base
      feat_enh_col  <- paste0(feat_enh,  '_enhancer')
      feat_prom_col <- paste0(feat_prom, '_promoter')

      if (feat_enh %in% colnames(peak_data)) {
        tmp <- peak_data[, c('peakID', feat_enh)]
        colnames(tmp)[2] <- feat_enh_col
        to_plot <- merge(to_plot, tmp, by.x = 'source', by.y = 'peakID')
      }
      if (feat_prom %in% colnames(peak_data)) {
        tmp <- peak_data[, c('peakID', feat_prom)]
        colnames(tmp)[2] <- feat_prom_col
        to_plot <- merge(to_plot, tmp, by.x = 'target', by.y = 'peakID')
      }
      to_plot[is.na(to_plot)] <- 'NA'

      to_plot_col <- to_plot %>%
        group_by(target, !!sym(feat_prom_col)) %>%
        summarise(
          source      = 'all',
          source_type = 'Enhancer',
          target_type = 'Promoter',
          !!feat_enh_col := ifelse(all(!!sym(feat_enh_col) == 'Yes'), 'Yes', 'No'),
          .groups = 'drop'
        ) %>%
        as.data.frame()

      lk_col <- to_plot_col %>% count(!!sym(feat_enh_col), !!sym(feat_prom_col), name = 'weight')
      colnames(lk_col) <- c('enhancer', 'promoter', 'weight')

      ep_res_c <- compute_or_p(lk_col, to_plot_col[[feat_enh_col]], to_plot_col[[feat_prom_col]])

      res_collapsed <- rbind(res_collapsed, data.frame(
        enhancer   = feat_enh_col, promoter = feat_prom_col,
        stratum    = mt_name,
        odds_ratio = ep_res_c$or, pvalue = ep_res_c$pv, adjpvalue = NA,
        stringsAsFactors = FALSE
      ))
    }
  }

  idx_mt <- res_collapsed$stratum == mt_name
  res_collapsed$adjpvalue[idx_mt] <- p.adjust(res_collapsed$pvalue[idx_mt], method = 'fdr')
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(res_collapsed, row.names = FALSE, file.path(out_dir, 'ep_enrichment_collapsed.csv'))
message('Saved: ', file.path(out_dir, 'ep_enrichment_collapsed.csv'), ' (', nrow(res_collapsed), ' rows)')
