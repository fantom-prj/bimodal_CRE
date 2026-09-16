#### Validate TF-to-peak GRN edges against ChIP-seq ground truth ####
#
# For two transcription factors with available ChIP-seq data (TEAD4,
# ONECUT2), checks how well each network variant's inferred TF-to-peak edges
# recover the TF's true ChIP-seq binding sites, across a range of deviance-
# ratio thresholds (the minimum model fit quality required to keep an edge).
# An edge counts as recovered whether it reaches the target directly
# (TF -> peak) or indirectly through an intermediate enhancer
# (TF -> enhancer -> peak); this script computes precision, recall, and F1
# for the "either direct or indirect" case across all 8 candidate network
# variants (base vs elastic-net-regularised, with/without TOBIAS footprinting,
# with/without Hi-C contacts).
#
# This produces the data behind Extended Figure 5a; see
# ../../Fig5/ext_fig5a_chipseq_f1.R for the plot itself.
#
# Reads (see ../config.R to set these paths):
#   [raw_data_folder]/tobias_chip_seq.merged.TEAD.ONECUT.tsv
#   [primary_data_folder]/3_preparing_base_network/base_network*.rds (4 files)
#   [primary_data_folder]/5_lasso_models/lasso_network*_with_randomization.rds (4 files)
#
# Writes: ../../Fig5/data/chipseq_f1_either_all.csv
#   columns: network, dev_ratio, TF, connection_type, target_type, metric, value
#   (kept only for connection_type == 'either', target_type == 'all',
#    metric == 'f1_score' -- the slice plotted in Extended Figure 5a)

rm(list = ls())
library(dplyr)

# ── CONTROL PANEL ─────────────────────────────────────────────────────────────
source('../config.R')  # defines primary_data_folder, raw_data_folder
chip_file <- file.path(raw_data_folder, 'tobias_chip_seq.merged.TEAD.ONECUT.tsv')
out_file  <- '../../Fig5/data/chipseq_f1_either_all.csv'

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

TFs                    <- c('TEAD4', 'ONECUT2')
thresholds             <- c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
significance_threshold <- 0.05
# ──────────────────────────────────────────────────────────────────────────────

n_networks <- length(networks_files)

compute_metrics <- function(retrieved, relevant) {
  retrieved <- unique(retrieved)
  relevant  <- unique(relevant)
  TP <- sum(retrieved %in% relevant)
  FP <- sum(!(retrieved %in% relevant))
  FN <- sum(!(relevant %in% retrieved))
  precision <- if ((TP + FP) > 0) TP / (TP + FP) else 0
  recall    <- if ((TP + FN) > 0) TP / (TP + FN) else 0
  f1_score  <- if ((precision + recall) > 0) 2 * precision * recall / (precision + recall) else 0
  c(true_positives = TP, false_positives = FP, false_negatives = FN,
    precision = precision, recall = recall, f1_score = f1_score)
}

df_colnames <- c('source', 'target', 'target_type', 'target_measurement_type',
                 'target_dev_ratio', 'direct_connection',
                 'enhancer', 'enhancer_measurement_type', 'enhancer_dev_ratio')

df_list <- vector('list', n_networks)
names(df_list) <- networks_to_compare

for (i in seq_len(n_networks)) {

  message('Loading network: ', networks_to_compare[i])
  network <- readRDS(networks_files[i])
  network$network_type <- networks_to_compare[i]

  if (!('dev_ratio' %in% colnames(network)))               network$dev_ratio <- 1
  if (!('target_measurement_type' %in% colnames(network))) network$target_measurement_type <- NA_character_

  if ('adjpvalue' %in% colnames(network)) {
    network <- network[network$adjpvalue <= significance_threshold, ]
  }

  df <- data.frame(matrix(NA_character_, 0, length(df_colnames)), stringsAsFactors = FALSE)
  colnames(df) <- df_colnames

  direct <- network %>%
    filter(source %in% TFs, link_type == 'TF_peak') %>%
    select(source, target, target_type, target_measurement_type, dev_ratio) %>%
    rename(target_dev_ratio = dev_ratio) %>%
    mutate(direct_connection = TRUE, enhancer = NA_character_,
           enhancer_measurement_type = NA_character_, enhancer_dev_ratio = 1)
  df <- rbind(df, direct)

  tf_to_enh <- network %>%
    filter(source %in% TFs, link_type == 'TF_peak', target_type == 'Enhancer') %>%
    select(source, target, target_measurement_type, dev_ratio) %>%
    rename(enhancer = target,
           enhancer_measurement_type = target_measurement_type,
           enhancer_dev_ratio = dev_ratio)

  enh_to_prom <- network %>%
    filter(source %in% tf_to_enh$enhancer, link_type == 'Enhancer_promoter',
           target_type == 'Promoter') %>%
    select(source, target, target_type, target_measurement_type, dev_ratio) %>%
    rename(enhancer = source, target_dev_ratio = dev_ratio)

  indirect <- merge(enh_to_prom, tf_to_enh, by = 'enhancer', all = FALSE)
  indirect$direct_connection <- FALSE
  df <- rbind(df, indirect[, df_colnames])

  df$network   <- networks_to_compare[i]
  df_list[[i]] <- df
}

df <- do.call(rbind, df_list)

#### Compute precision / recall / F1 ####

chip_data <- read.delim(chip_file)
chip_data <- chip_data[, c('TF', 'V4')]
colnames(chip_data) <- c('TF', 'target')

res_colnames <- c('network', 'dev_ratio', 'TF', 'connection_type', 'target_type',
                  'metric', 'value')
res <- data.frame(matrix(NA, 0, length(res_colnames)), stringsAsFactors = FALSE)
colnames(res) <- res_colnames

for (i in seq_len(n_networks)) {

  network_name <- networks_to_compare[i]
  network_df   <- df[df$network == network_name, ]

  for (k in seq_along(thresholds)) {

    threshold <- thresholds[k]
    idx <- network_df$target_dev_ratio >= threshold &
           network_df$enhancer_dev_ratio >= threshold
    nd  <- network_df[idx, ]

    for (j in seq_along(TFs)) {

      TF       <- TFs[j]
      tf_nd    <- nd[nd$source == TF, ]
      chip_tf  <- chip_data[chip_data$TF == TF, ]
      relevant <- paste0(chip_tf$TF, '_', chip_tf$target)

      add_row <- function(sub_df, conn_type, tgt_type, use_enhancer = FALSE) {
        ids <- if (use_enhancer) paste0(sub_df$source, '_', sub_df$enhancer)
               else              paste0(sub_df$source, '_', sub_df$target)
        m <- compute_metrics(ids, relevant)
        res <<- rbind(res, data.frame(
          network = network_name, dev_ratio = threshold, TF = TF,
          connection_type = conn_type, target_type = tgt_type,
          metric = names(m), value = m, stringsAsFactors = FALSE
        ))
      }

      # 'either' = edge counts whether recovered via a direct TF->peak link
      # or an indirect TF->enhancer->peak path; 'all' = all target types
      # pooled together (promoters and enhancers alike).
      add_row(tf_nd, 'either', 'all')
    }
  }
}

rownames(res) <- NULL
message('Total rows computed: ', nrow(res))

res_f1_either_all <- res[res$metric == 'f1_score' &
                          res$connection_type == 'either' &
                          res$target_type == 'all', ]

dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
write.csv(res_f1_either_all, row.names = FALSE, file = out_file)
message('Saved: ', out_file, ' (', nrow(res_f1_either_all), ' rows)')
