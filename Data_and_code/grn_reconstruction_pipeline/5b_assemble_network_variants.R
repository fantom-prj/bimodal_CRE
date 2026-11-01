#### Assemble the lasso network variants from separately fitted parts ####
#
# Companion to 5a_lasso_models.R. Script 5 fits two independent families of
# models: TF-to-peak models (TFs as the only predictors) and
# enhancer-to-promoter models (enhancers as the only predictors). The two
# network flags each affect only one family:
#   - use_hic changes only the enhancer-to-promoter candidate edges;
#   - use_footprinting changes only the TF-to-peak candidate edges;
#   - use_enhancer_tCRE = FALSE drops the enhancer tCRE target models from the
#     TF-to-peak part, and changes the predictor set of the
#     enhancer-to-promoter models.
# Every model's input values, target list, target order and permutation seeds
# are the same in every variant that shares a part, and p-values are
# BH-adjusted within each model. So a part fitted for one variant is identical
# to the same part of any other variant that shares it, and only the
# following parts need fitting beyond the full published network:
#
#   network (output file)                                  TF_peak part            Enhancer_promoter part
#   lasso_network_with_randomization                       full run of script 5 (both parts)
#   lasso_network_with_randomization_no_enhancer_tCRE      published, minus        fit: use_enhancer_tCRE = FALSE
#                                                          enhancer tCRE targets
#   lasso_network_no_hic_with_randomization                published               fit: use_hic = FALSE
#   lasso_network_no_footprinting_with_randomization       fit: use_footprinting   published
#                                                          = FALSE
#   lasso_network_no_hic_no_footprinting_with_randomization   as no_footprinting   as no_hic
#
# The parts are fitted with script 5's fit_edge_classes control-panel option,
# which writes them as ..._TF_peak_part.rds / ..._Enhancer_promoter_part.rds.
# Each network is rbind(TF_peak part, Enhancer_promoter part), the same order
# script 5 uses.
#
# If any of the four inputs below was fitted in target blocks (script 5's
# target_blocks option), it is first rebuilt from its
# <name>_block<b>of<B>.rds files: the per-target results are put back in
# target order and bound exactly as script 5 binds a single run's. An input
# is rebuilt when it is missing or older than its newest block. The variants
# are assembled once all four inputs exist, so the published network can be
# rebuilt and checked on its own first.
#
# Reads from output/5_lasso_models/:
#   lasso_network_with_randomization.rds
#   lasso_network_with_randomization_no_enhancer_tCRE_Enhancer_promoter_part.rds
#   lasso_network_no_hic_with_randomization_Enhancer_promoter_part.rds
#   lasso_network_no_footprinting_with_randomization_TF_peak_part.rds
#
# Writes into output/5_lasso_models/:
#   lasso_network_with_randomization_no_enhancer_tCRE.rds
#   lasso_network_no_hic_with_randomization.rds
#   lasso_network_no_footprinting_with_randomization.rds
#   lasso_network_no_hic_no_footprinting_with_randomization.rds

rm(list = ls())

# ── CONTROL PANEL ─────────────────────────────────────────────────────────────
source('../config.R')  # defines primary_data_folder

results_folder <- file.path(primary_data_folder, '5_lasso_models')
# ──────────────────────────────────────────────────────────────────────────────

read_network <- function(name, link_type = NULL) {
  network <- readRDS(file.path(results_folder, paste0(name, '.rds')))
  if (!is.null(link_type)) {
    stopifnot(all(network$link_type == link_type))
  }
  network
}

assemble <- function(tf_peak_part, enhancer_promoter_part) {
  stopifnot(all(tf_peak_part$link_type == 'TF_peak'),
            all(enhancer_promoter_part$link_type == 'Enhancer_promoter'))
  network <- rbind(tf_peak_part, enhancer_promoter_part)
  rownames(network) <- NULL
  network
}

combine_blocks <- function(name) {
  pattern <- paste0('^', name, '_block([0-9]+)of([0-9]+)\\.rds$')
  files <- list.files(results_folder, pattern)
  output_file <- file.path(results_folder, paste0(name, '.rds'))
  if (length(files) == 0 ||
      (file.exists(output_file) &&
       file.mtime(output_file) > max(file.mtime(file.path(results_folder, files))))) {
    return(invisible(NULL))
  }
  n_blocks <- as.integer(unique(sub(pattern, '\\2', files)))
  if (length(n_blocks) != 1 ||
      !setequal(as.integer(sub(pattern, '\\1', files)), seq_len(n_blocks))) {
    stop(name, ': expected blocks 1..B of a single split, found ', paste(files, collapse = ', '))
  }
  blocks <- lapply(file.path(results_folder, files), readRDS)
  parts <- lapply(c('TF_peak', 'Enhancer_promoter'), function(edge_class) {
    tmp <- do.call('c', lapply(blocks, `[[`, edge_class))
    tmp <- unname(tmp[order(as.integer(names(tmp)))])
    # as in script 5: drop the targets without a model, then bind
    errors <- which(sapply(tmp, class) == 'character')
    if (length(errors) > 0) {
      tmp <- tmp[-errors]
    }
    do.call('rbind', tmp)
  })
  network <- rbind(parts[[1]], parts[[2]])
  RNGkind("L'Ecuyer-CMRG")
  set.seed(42)
  idx <- sample(1:dim(network)[1], min(1000, dim(network)[1]))
  write.csv(network[idx, ], row.names = FALSE,
            file.path(results_folder, paste0(name, '_sample.csv')))
  saveRDS(network, output_file)
  counts <- table(network$link_type)
  message(name, ' (rebuilt from ', n_blocks, ' blocks): ',
          paste(names(counts), counts, sep = ' ', collapse = ', '))
}

inputs <- c('lasso_network_with_randomization',
            'lasso_network_with_randomization_no_enhancer_tCRE_Enhancer_promoter_part',
            'lasso_network_no_hic_with_randomization_Enhancer_promoter_part',
            'lasso_network_no_footprinting_with_randomization_TF_peak_part')
for (name in inputs) {
  combine_blocks(name)
}
missing <- inputs[!file.exists(file.path(results_folder, paste0(inputs, '.rds')))]
if (length(missing) > 0) {
  message('Variants not assembled yet, missing: ', paste(missing, collapse = ', '))
} else {
  published <- read_network('lasso_network_with_randomization')
  published_tf_peak <- published[published$link_type == 'TF_peak', ]
  published_enhancer_promoter <- published[published$link_type == 'Enhancer_promoter', ]

  no_enhancer_tCRE_enhancer_promoter <- read_network(
    'lasso_network_with_randomization_no_enhancer_tCRE_Enhancer_promoter_part', 'Enhancer_promoter')
  no_hic_enhancer_promoter <- read_network(
    'lasso_network_no_hic_with_randomization_Enhancer_promoter_part', 'Enhancer_promoter')
  no_footprinting_tf_peak <- read_network(
    'lasso_network_no_footprinting_with_randomization_TF_peak_part', 'TF_peak')

  networks <- list(
    lasso_network_with_randomization_no_enhancer_tCRE = assemble(
      published_tf_peak[!(published_tf_peak$target_type == 'Enhancer' &
                            published_tf_peak$target_measurement_type == 'tCRE'), ],
      no_enhancer_tCRE_enhancer_promoter),
    lasso_network_no_hic_with_randomization = assemble(
      published_tf_peak, no_hic_enhancer_promoter),
    lasso_network_no_footprinting_with_randomization = assemble(
      no_footprinting_tf_peak, published_enhancer_promoter),
    lasso_network_no_hic_no_footprinting_with_randomization = assemble(
      no_footprinting_tf_peak, no_hic_enhancer_promoter))

  for (name in names(networks)) {
    saveRDS(networks[[name]], file.path(results_folder, paste0(name, '.rds')))
    counts <- table(networks[[name]]$link_type)
    message(name, ': ', paste(names(counts), counts, sep = ' ', collapse = ', '))
  }
}
