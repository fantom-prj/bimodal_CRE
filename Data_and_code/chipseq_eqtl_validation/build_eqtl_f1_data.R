#### Regeneration script: eQTL validation F1 curve data ####
#
# For enhancer-to-promoter GRN edges, checks how well each network variant's
# inferred links recover known enhancer-to-gene regulatory relationships from
# fine-mapped eQTLs (credible-set SNPs in the enhancer, linked to their eGene).
# An inferred edge counts as a true positive if it connects the same
# peak/gene pair as an eQTL association; precision, recall, and F1 are
# computed across a range of deviance-ratio thresholds, separately for each
# of 8 candidate network variants (base vs elastic-net-regularised,
# with/without TOBIAS footprinting, with/without Hi-C contacts) and for 5
# source/target measurement-type pairs (aCRE/tCRE combinations).
#
# This produces the data behind Extended Figure 5b; see
# ../../Fig5/ext_fig5b_eqtl_f1.R for the plot itself, which uses only the
# pooled 'All interactions' slice.
#
# Reads (see ../config.R to set these paths):
#   [raw_data_folder]/aCRE2egene_eqtl_LD_pip3_all_method_pairbase.tsv
#   [primary_data_folder]/3_preparing_base_network/base_network*.rds (4 files)
#   [primary_data_folder]/5_lasso_models/lasso_network*_with_randomization.rds (4 files)
#
# Writes: ../../Fig5/data/eqtl_f1_all_interactions.csv
#   columns: network, dev_ratio, measurement_type, metric, value
#   (pre-filtered to metric == 'f1', measurement_type == 'All interactions')

rm(list = ls())
library(org.Hs.eg.db)

# ── CONTROL PANEL ─────────────────────────────────────────────────────────────
source('../config.R')  # defines primary_data_folder, raw_data_folder
eqtl_file <- file.path(raw_data_folder, 'aCRE2egene_eqtl_LD_pip3_all_method_pairbase.tsv')
out_file  <- '../../Fig5/data/eqtl_f1_all_interactions.csv'

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

thresholds <- c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
measurement_types <- list(
  all       = list(c('aCRE', 'tCRE'), c('aCRE', 'tCRE')),
  aCRE_aCRE = list('aCRE', 'aCRE'),
  aCRE_tCRE = list('aCRE', 'tCRE'),
  tCRE_tCRE = list('tCRE', 'tCRE'),
  tCRE_aCRE = list('tCRE', 'aCRE')
)
significance_threshold <- 0.05
# ──────────────────────────────────────────────────────────────────────────────

n_networks <- length(networks_files)

#### Load and format eQTL ground truth ####

eqtl_raw <- read.table(eqtl_file, header = TRUE)
eqtl_raw <- eqtl_raw[eqtl_raw$promoter_type2 == 'enhancer-like', ]
eqtl_raw <- eqtl_raw[, c('peakID', 'eqtl_gene_id')]
colnames(eqtl_raw) <- c('peak_id', 'ensembl_gene_id')
eqtl_raw$peak_id <- gsub('-', '_', eqtl_raw$peak_id)

sym_map <- AnnotationDbi::select(org.Hs.eg.db,
                         keys    = unique(eqtl_raw$ensembl_gene_id),
                         keytype = 'ENSEMBL',
                         columns = 'SYMBOL')
eqtl <- merge(eqtl_raw, sym_map, by.x = 'ensembl_gene_id', by.y = 'ENSEMBL')
eqtl$ensembl_gene_id <- NULL
eqtl <- unique(eqtl)
eqtl$id <- paste0(eqtl$peak_id, '_', eqtl$SYMBOL)

message('eQTL ground-truth pairs: ', length(unique(eqtl$id)))

#### Compute precision / recall / F1 ####

col_names <- c('network', 'dev_ratio', 'measurement_type', 'metric', 'value')
results   <- data.frame(matrix(NA, 0, length(col_names)), stringsAsFactors = FALSE)
colnames(results) <- col_names

tp_col_labels <- sapply(measurement_types, function(x) paste(x, collapse = '_'))

for (i in seq_len(n_networks)) {

  message('Processing network: ', networks_to_compare[i])
  network <- readRDS(networks_files[i])

  if (!('dev_ratio' %in% colnames(network))) network$dev_ratio <- 1
  if ('adjpvalue' %in% colnames(network)) {
    network <- network[network$adjpvalue <= significance_threshold, ]
  }

  network <- network[network$link_type == 'Enhancer_promoter', ]
  network$id <- paste0(network$source, '_', network$target_gene)

  for (j in seq_along(thresholds)) {

    thr     <- thresholds[j]
    net_thr <- network[network$dev_ratio >= thr, ]

    for (k in seq_along(measurement_types)) {

      if ('source_measurement_type' %in% colnames(net_thr)) {
        idx <- net_thr$source_measurement_type %in% measurement_types[[k]][[1]] &
               net_thr$target_measurement_type %in% measurement_types[[k]][[2]]
        cur <- net_thr[idx, ]
        if (names(measurement_types)[k] != 'all') {
          cur <- cur[!(cur$id %in% net_thr$id[!idx]), ]
        }
      } else {
        cur <- net_thr
      }

      actual    <- unique(eqtl$id)
      predicted <- unique(cur$id)
      TP        <- sum(actual %in% predicted)
      FP        <- length(predicted) - TP
      FN        <- length(actual)    - TP
      precision <- if ((TP + FP) > 0) TP / (TP + FP) else 0
      recall    <- if ((TP + FN) > 0) TP / (TP + FN) else 0
      F1        <- if ((precision + recall) > 0) 2 * precision * recall / (precision + recall) else 0

      to_add <- data.frame(
        network          = networks_to_compare[i],
        dev_ratio        = thr,
        measurement_type = names(measurement_types)[k],
        metric           = c('precision', 'recall', 'f1', 'tp'),
        value            = c(precision, recall, F1, TP),
        stringsAsFactors = FALSE
      )
      results <- rbind(results, to_add)
    }
  }
}

results$value[is.na(results$value)] <- 0

mt_labels <- c(
  all       = 'All interactions',
  aCRE_aCRE = 'aCRE / aCRE interactions',
  aCRE_tCRE = 'aCRE / tCRE interactions',
  tCRE_tCRE = 'tCRE / tCRE interactions',
  tCRE_aCRE = 'tCRE / aCRE interactions'
)
results$measurement_type <- mt_labels[results$measurement_type]

results_f1_all <- results[results$metric == 'f1' &
                           results$measurement_type == 'All interactions', ]

dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
write.csv(results_f1_all, row.names = FALSE, file = out_file)
message('Saved: ', out_file, ' (', nrow(results_f1_all), ' rows)')
