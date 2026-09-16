#### Call TF-to-peak edges from TOBIAS footprinting scores ####
#
# Step 2 of the GRN reconstruction pipeline. For every transcription factor
# with a known motif, pools its TOBIAS BINDetect footprinting scores across
# all peaks and calls a peak "bound" by that TF if its footprinting score
# clears a threshold. Two modes, controlled by `use_footprinting`:
#   TRUE  — threshold is the per-TF `quantile_score_threshold` quantile of
#           that TF's own score distribution (footprinting-informed).
#   FALSE — threshold is the minimum observed score, i.e. every peak with a
#           motif match for the TF is called bound (motif-only, no
#           footprinting evidence required). This is the variant used
#           throughout the published network (see 3_preparing_base_network.R).
# `quantile_score_threshold = 0.82` is the operating point chosen in
# 1_tf_2_peak_set_threshold.R (the lowest score quantile at which both
# validated TFs, TEAD4 and ONECUT2, reach recall > 0.85 against ChIP-seq).
#
# Reads (see ../config.R to set raw_data_folder):
#   [raw_data_folder]/expressed_400motif.tsv
#     motif table mapping each JASPAR-style motif label to its TF gene
#     symbol(s) (two columns for heterodimeric motifs)
#   [raw_data_folder]/BINDetect_compare_output/<motif-label>/*_overview.txt
#     TOBIAS BINDetect per-peak footprinting scores, one folder per motif
#
# Writes into output/2_tf_2_peak/:
#   tf2peak.rds / tf2peak.csv               (use_footprinting = TRUE)
#   tf2peak_no_footprinting.rds / .csv      (use_footprinting = FALSE)
#     columns: tf, seqnames, start, end, max_score
#     (one row per TF-peak pair called bound; peaks with several motifs
#      belonging to the same TF are collapsed by keeping the max score)

rm(list = ls())
library(dplyr)
library(foreach)
library(doParallel)

# ── CONTROL PANEL ─────────────────────────────────────────────────────────────
source('../config.R')  # defines raw_data_folder, primary_data_folder

use_footprinting <- FALSE
quantile_score_threshold <- 0.82
num_cpus <- 10
tfs_file <- file.path(raw_data_folder, 'expressed_400motif.tsv')
tobias_res_location <- file.path(raw_data_folder, 'BINDetect_compare_output')
out_dir  <- file.path(primary_data_folder, '2_tf_2_peak')
# ──────────────────────────────────────────────────────────────────────────────

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
registerDoParallel(num_cpus)

# motif table -> unique TF gene symbols (heterodimers list two TF columns)
motifs_table <- read.table(tfs_file, header = TRUE)
tfs <- unique(c(motifs_table$Approved.symbol.TF1, motifs_table$Approved.symbol.TF2))
tfs <- tfs[!is.na(tfs)]

tobias_folders <- dir(tobias_res_location)

tf2peak <- foreach(i = 1:length(tfs), .packages = c('dplyr')) %dopar% {

  print(paste0('tf ', i, ' out of ', length(tfs)))
  tf <- tfs[i]

  # every motif folder whose TF1 or TF2 column matches this TF
  current_tobias_folders <- motifs_table$label[motifs_table$Approved.symbol.TF1 == tf]
  current_tobias_folders <- c(current_tobias_folders,
                              motifs_table$label[!is.na(motifs_table$Approved.symbol.TF2) &
                                                   motifs_table$Approved.symbol.TF2 == tf])
  current_tobias_folders <- unique(current_tobias_folders)

  if (length(current_tobias_folders) == 0) {
    return(NULL)
  }

  tmp_res <- data.frame(tf = c(), seqnames = c(), start = c(), end = c(), max_score = c())

  for (j in 1:length(current_tobias_folders)) {

    tmp <- gsub('::', '', current_tobias_folders[j])
    overview_file <- file.path(tobias_res_location, tmp, paste0(tmp, '_overview.txt'))
    if (!file.exists(overview_file)) next()
    t_res <- read.table(overview_file, header = TRUE, sep = '\t')

    score_values <- t_res[, grep('score', colnames(t_res))]
    score_values$TFBS_score <- NULL
    score_values <- as.matrix(score_values)
    t_res$max_score <- apply(score_values, 1, max, na.rm = TRUE)

    t_res <- unique(t_res[, c('peak_chr', 'peak_start', 'peak_end', 'max_score')])
    # a peak can appear more than once with different max_scores; keep the max
    t_res <- t_res %>%
      group_by(peak_chr, peak_start, peak_end) %>%
      summarise(max_score = max(max_score)) %>%
      as.data.frame()

    if (use_footprinting) {
      threshold <- quantile(t_res$max_score[t_res$max_score > 0], quantile_score_threshold)
    } else {
      threshold <- min(t_res$max_score, na.rm = TRUE)
    }

    is_it_bound <- t_res$max_score >= threshold
    bound_peaks <- unique(t_res[which(is_it_bound), ])

    colnames(bound_peaks) <- c('seqnames', 'start', 'end', 'max_score')
    bound_peaks <- cbind(tf = tf, bound_peaks)
    bound_peaks <- unique(bound_peaks)
    tmp_res <- rbind(tmp_res, bound_peaks)
  }

  return(tmp_res)
}

names(tf2peak) <- sapply(tf2peak, function(x) unique(x$tf))
tf2peak <- tf2peak[sapply(tf2peak, nrow) > 0]
tf2peak <- do.call('rbind', tf2peak)
rownames(tf2peak) <- 1:dim(tf2peak)[1]
tf2peak <- unique(tf2peak)

# collapse duplicate (tf, peak) pairs across motif folders, keeping the max score
tf2peak <- tf2peak %>%
  group_by(tf, seqnames, start, end) %>%
  summarise(max_score = max(max_score)) %>%
  as.data.frame()

message('TF-to-peak edges called: ', nrow(tf2peak), ' across ', length(unique(tf2peak$tf)), ' TFs')

out_suffix <- ifelse(use_footprinting, '', '_no_footprinting')
saveRDS(tf2peak, file = file.path(out_dir, paste0('tf2peak', out_suffix, '.rds')))
write.csv(tf2peak, row.names = FALSE, file = file.path(out_dir, paste0('tf2peak', out_suffix, '.csv')))
