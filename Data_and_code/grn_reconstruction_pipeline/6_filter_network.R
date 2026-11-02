#### Step 6: edge thresholds and the published high-confidence GRN ####
#
# Derives the edge-filtering thresholds from the published network variant
# (lasso_network_with_randomization.rds, from step 5a), writes them to
# ../network_loading/network_thresholds.txt, and applies them to write
# lasso_network_filtered.rds, the network shipped in this repository:
#   adjpvalue <= adjpvalue_threshold, then dev_ratio >= dev_ratio_threshold,
#   then |estimate| >= the per-stratum threshold.
# The per-stratum |estimate| threshold is the 90th percentile of |estimate|
# within each stratum (TF_peak: target type x target measurement type;
# Enhancer_promoter: source x target measurement type), computed on the
# edges that pass the first two filters and rounded to 3 decimals.
#
# Every downstream script reads network_thresholds.txt through
# read_thresholds() in ../network_loading/load_and_filter_network.R, so to
# change a cutoff, edit the control panel below and rerun this script.
#
# Reads:
#   [primary_data_folder]/5_lasso_models/lasso_network_with_randomization.rds
# Writes:
#   ../network_loading/network_thresholds.txt
#   [primary_data_folder]/5_lasso_models/lasso_network_filtered.rds
#
# Cross-check (published network): 89,371 TF_peak + 52,362 Enhancer_promoter
# = 141,733 edges.

rm(list = ls())

# ── CONTROL PANEL ─────────────────────────────────────────────────────────────
source('../config.R')  # defines primary_data_folder
network_file        <- file.path(primary_data_folder,
                                 '5_lasso_models/lasso_network_with_randomization.rds')
thresholds_file     <- '../network_loading/network_thresholds.txt'
out_folder          <- file.path(primary_data_folder, '5_lasso_models')
adjpvalue_threshold <- 0.05
dev_ratio_threshold <- 0.3
estimate_percentile <- 0.9
# ──────────────────────────────────────────────────────────────────────────────

source('../network_loading/load_and_filter_network.R')

message('Loading network...')
network <- readRDS(network_file)
message('  Edges: ', nrow(network))

#### Per-stratum |estimate| thresholds ####

net_pd  <- network[network$adjpvalue <= adjpvalue_threshold, ]
net_pd  <- net_pd[net_pd$dev_ratio >= dev_ratio_threshold, ]
stratum <- ifelse(net_pd$link_type == 'TF_peak',
                  paste0('TF_peak_', net_pd$target_type, '_', net_pd$target_measurement_type),
                  paste0('Enhancer_promoter_', net_pd$source_measurement_type, '_',
                         net_pd$target_measurement_type))
p90 <- tapply(abs(net_pd$estimate), stratum,
              function(x) unname(quantile(x, probs = estimate_percentile)))

ep_strata <- paste0('Enhancer_promoter_', c('aCRE_aCRE', 'aCRE_tCRE', 'tCRE_aCRE', 'tCRE_tCRE'))
tf_strata <- paste0('TF_peak_', c('Enhancer_aCRE', 'Enhancer_tCRE', 'Promoter_aCRE', 'Promoter_tCRE'))
stopifnot(all(c(ep_strata, tf_strata) %in% names(p90)))
threshold_line <- function(s) sprintf('abs_estimate_%s = %.3f', s, p90[[s]])

#### Write network_thresholds.txt ####

pctl_label <- paste0('p', round(100 * estimate_percentile))
writeLines(c(
  '# Network edge filtering thresholds',
  '# ============================================================',
  '# Written by grn_reconstruction_pipeline/6_filter_network.R from',
  '# lasso_network_with_randomization.rds -- to change a cutoff, edit that',
  '# script\'s control panel and rerun it rather than editing this file.',
  '# Load these values with read_thresholds() from load_and_filter_network.R.',
  '#',
  '#   adjpvalue_threshold  — permutation-test adjusted p-value from lasso randomisation',
  '#   dev_ratio_threshold  — proportion of deviance explained by the edge in the lasso model',
  paste0('#   abs_estimate_thresholds — ', round(100 * estimate_percentile),
         'th percentile of |estimate| within each edge-type stratum,'),
  '#                             computed AFTER applying adjpvalue and dev_ratio filters.',
  '#                             Stratum keys: <source_measurement_type>_<target_measurement_type>',
  '#                             for Enhancer_promoter edges, and',
  '#                             TF_<target_type>_<target_measurement_type> for TF_peak edges.',
  '# ============================================================',
  '',
  '# ── Global filters (applied to all lasso network edges) ──────────────────────',
  paste0('adjpvalue_threshold = ', adjpvalue_threshold),
  paste0('dev_ratio_threshold = ', dev_ratio_threshold),
  '',
  paste0('# ── Per-stratum |estimate| thresholds (', pctl_label,
         ', Enhancer_promoter edges) ─────────'),
  '# Stratum: enhancer_source_modality → promoter_target_modality',
  sapply(ep_strata, threshold_line),
  '',
  paste0('# ── Per-stratum |estimate| thresholds (', pctl_label,
         ', TF_peak edges) ───────────────────'),
  '# Stratum: TF → target_type (Enhancer/Promoter) × target_measurement_type',
  sapply(tf_strata, threshold_line),
  '',
  '# ── Notes ─────────────────────────────────────────────────────────────────────',
  '# - chipseq_eqtl_validation/ sweeps dev_ratio from 0 to 0.8 to draw performance',
  '#   curves; dev_ratio_threshold above marks the operating point highlighted in',
  '#   its summary figures.',
  '# - The abs_estimate thresholds were cross-validated against a collaborator\'s',
  '#   independent elbow-point analysis on a previous network version; the p90',
  '#   criterion reproduces their values within ±0.05 for all strata.'
), thresholds_file)
message('Saved: ', thresholds_file)

#### Filter with the thresholds just written, and save ####

# Filtering through the written file (rather than the unrounded values above)
# keeps the shipped network exactly reproducible by every downstream script.
thresholds <- read_thresholds(thresholds_file)
filtered   <- filter_network(network, thresholds)
saveRDS(filtered, file.path(out_folder, 'lasso_network_filtered.rds'))

counts <- table(filtered$link_type)
message('Saved: lasso_network_filtered.rds (',
        paste(names(counts), counts, sep = ' = ', collapse = ', '),
        '; total ', nrow(filtered), ')')
