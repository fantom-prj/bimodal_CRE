#### Calibrate a TOBIAS footprinting-score threshold from ChIP-seq ####
#
# Step 1 of the GRN reconstruction pipeline. For each of two transcription
# factors with available ChIP-seq ground truth (TEAD4, ONECUT2), scores every
# accessible peak by its maximum TOBIAS footprinting score across the TF's
# motif(s), then compares each candidate score threshold against the ChIP-seq
# binding calls to compute precision, recall, and F1. The resulting
# precision-recall curve is used to choose the footprinting-score quantile
# threshold applied when building TF-to-peak edges in
# 2_tf_2_peak.R (see that script's `quantile_score_threshold` control-panel
# value, chosen as the lowest threshold at which both TFs' recall exceeds
# 0.85).
#
# Reads (see ../config.R to set raw_data_folder):
#   [raw_data_folder]/BINDetect_compare_output/<TF-motif-folder>/*_overview.txt
#     (TOBIAS BINDetect per-peak footprinting scores, one folder per motif)
#   [raw_data_folder]/tobias_chip_seq.merged.TEAD.ONECUT.tsv
#     (ChIP-seq peak calls for TEAD4 and ONECUT2, used as ground truth)
#
# Writes, one set per TF, into output/1_tf_2_peak_set_threshold/:
#   <TF>.png   precision-recall curve across score quantile thresholds
#   <TF>.rds / <TF>.csv
#     columns: precision, recall, false_positive_rate, false_negative_rate,
#     f1, quantile, threshold (one row per score_quantile_thresholds value)

rm(list = ls())
library(dplyr)
library(ggplot2)

# ── CONTROL PANEL ─────────────────────────────────────────────────────────────
source('../config.R')  # defines raw_data_folder, primary_data_folder

score_quantile_thresholds <- seq(from = 0.75, to = 0.99, by = 0.01)
tfs <- c('TEAD4', 'ONECUT2')
tobias_root_folder <- file.path(raw_data_folder, 'BINDetect_compare_output')
chip_res_file       <- file.path(raw_data_folder, 'tobias_chip_seq.merged.TEAD.ONECUT.tsv')
out_dir             <- file.path(primary_data_folder, '1_tf_2_peak_set_threshold')
# ──────────────────────────────────────────────────────────────────────────────

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ChIP-seq ground truth: one row per confirmed TF-binding peak
all_chip_res <- read.table(chip_res_file, header = TRUE)
names(all_chip_res) <- c('peak_id', 'promoter_type', 'TF')
all_chip_res$seqnames <- sapply(all_chip_res$peak_id, function(x) strsplit(x, '_')[[1]][1])
all_chip_res$start    <- sapply(all_chip_res$peak_id, function(x) strsplit(x, '_')[[1]][2])
all_chip_res$end      <- sapply(all_chip_res$peak_id, function(x) strsplit(x, '_')[[1]][3])

for (tf in tfs) {

  chip_res <- all_chip_res[all_chip_res$TF == tf, ]

  # load TOBIAS BINDetect per-peak footprinting scores for this TF's motif(s)
  tobias_folders <- dir(tobias_root_folder)
  current_tobias_folder <- tobias_folders[grep(tf, tobias_folders)]
  t_res <- read.table(file.path(tobias_root_folder,
                                current_tobias_folder,
                                paste0(current_tobias_folder, '_overview.txt')),
                      header = TRUE, sep = '\t')

  # per-peak maximum footprinting score across all score columns
  score_values <- t_res[, grep('score', colnames(t_res))]
  score_values$TFBS_score <- NULL
  score_values <- as.matrix(score_values)
  t_res$max_score_values <- apply(score_values, 1, max, na.rm = TRUE)

  tfbs <- unique(t_res[, c('peak_chr', 'peak_start', 'peak_end', 'max_score_values')])
  colnames(tfbs) <- c('seqnames', 'start', 'end', 'max_score')
  tfbs <- tfbs %>%
    group_by(seqnames, start, end) %>%
    summarise(max_score = max(max_score)) %>%
    as.data.frame()
  tfbs$start   <- tfbs$start + 1  # BED half-open start -> 1-based
  tfbs$peak_id <- paste0(tfbs$seqnames, '_', tfbs$start, '_', tfbs$end)

  tfbs$chip_seq_confirmed <- 0
  tfbs$chip_seq_confirmed[tfbs$peak_id %in% chip_res$peak_id] <- 1

  # precision/recall/F1 at each candidate score-quantile threshold
  thresholds <- quantile(tfbs$max_score[tfbs$max_score > 0], score_quantile_thresholds)
  names(thresholds) <- score_quantile_thresholds
  res <- data.frame(
    precision = rep(NA, length(thresholds)),
    recall    = rep(NA, length(thresholds)),
    false_positive_rate = rep(NA, length(thresholds)),
    false_negative_rate = rep(NA, length(thresholds)),
    f1        = rep(NA, length(thresholds)),
    quantile  = 1 - score_quantile_thresholds,
    threshold = thresholds
  )
  for (i in seq_along(thresholds)) {
    threshold <- thresholds[i]
    TP <- sum(tfbs$max_score >= threshold & tfbs$chip_seq_confirmed == 1)
    FP <- sum(tfbs$max_score >= threshold & tfbs$chip_seq_confirmed == 0)
    FN <- sum(tfbs$max_score <  threshold & tfbs$chip_seq_confirmed == 1)
    TN <- sum(tfbs$max_score <  threshold & tfbs$chip_seq_confirmed == 0)
    res$precision[i] <- TP / (TP + FP)
    res$recall[i]    <- TP / (TP + FN)
    res$false_positive_rate[i] <- FP / (TN + FP)
    res$false_negative_rate[i] <- FN / (TP + FN)
    res$f1[i] <- 2 * (res$precision[i] * res$recall[i]) / (res$precision[i] + res$recall[i])
  }

  tmp <- which.max(res$f1)[1]
  max_f1_precision <- res$precision[tmp]
  max_f1_recall    <- res$recall[tmp]
  res$label <- paste0('Q: ', round(res$quantile, 2), ', S: ', round(res$threshold, 2))

  p <- ggplot(res, aes(x = recall, y = precision, label = label)) +
    geom_line() +
    geom_point() +
    geom_text(nudge_x = 0.07) +
    geom_vline(xintercept = max_f1_recall, color = 'darkred', linetype = 'dashed') +
    geom_hline(yintercept = max_f1_precision, color = 'darkred', linetype = 'dashed') +
    scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.05)) +
    scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.05)) +
    theme_bw()

  png(filename = file.path(out_dir, paste0(tf, '.png')), width = 3300, height = 3300, res = 300)
  plot(p)
  dev.off()

  saveRDS(res, file.path(out_dir, paste0(tf, '.rds')))
  write.csv(res, row.names = FALSE, file = file.path(out_dir, paste0(tf, '.csv')))
  message('Saved threshold curve for ', tf, ' (best F1 at quantile ',
          round(res$quantile[tmp], 2), ', threshold ', round(res$threshold[tmp], 3), ')')
}
