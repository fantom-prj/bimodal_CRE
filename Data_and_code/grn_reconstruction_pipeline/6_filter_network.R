#### Step 6: filter the lasso network down to the published high-confidence GRN ####
#
# Applies the three confidence filters used by every figure to the published
# network variant (lasso_network_with_randomization.rds, from step 5a) and
# writes lasso_network_filtered.rds, the network shipped in this repository:
#   adjpvalue <= adjpvalue_threshold, then dev_ratio >= dev_ratio_threshold,
#   then |estimate| >= the stored per-stratum 90th percentile.
# Thresholds are read from ../network_loading/network_thresholds.txt and
# applied by filter_network() in ../network_loading/load_and_filter_network.R.
#
# The per-stratum |estimate| thresholds stored in network_thresholds.txt are
# the 90th percentile of |estimate| within each stratum, computed on this
# network after the adjpvalue and dev_ratio filters and rounded to 3
# decimals. This script recomputes them and writes them, next to the stored
# values, to abs_estimate_p90_recomputed.tsv; it warns if any differ. After
# refitting step 5a, copy the recomputed values into network_thresholds.txt
# and rerun this script. The thresholds file is never edited automatically.
#
# Reads:
#   [primary_data_folder]/5_lasso_models/lasso_network_with_randomization.rds
#   ../network_loading/network_thresholds.txt
# Writes:
#   [primary_data_folder]/5_lasso_models/lasso_network_filtered.rds
#   [primary_data_folder]/5_lasso_models/abs_estimate_p90_recomputed.tsv
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
estimate_percentile <- 0.9
# ──────────────────────────────────────────────────────────────────────────────

source('../network_loading/load_and_filter_network.R')
thresholds <- read_thresholds(thresholds_file)

message('Loading network...')
network <- readRDS(network_file)
message('  Edges: ', nrow(network))

#### Recompute the per-stratum |estimate| thresholds ####

# Same order as filter_network(): the percentile is taken on the edges that
# pass the adjpvalue and dev_ratio filters.
net_pd <- network[network$adjpvalue <= thresholds$adjpvalue_threshold, ]
net_pd <- net_pd[net_pd$dev_ratio >= thresholds$dev_ratio_threshold, ]
stratum <- ifelse(net_pd$link_type == 'TF_peak',
                  paste0('TF_', net_pd$target_type, '_', net_pd$target_measurement_type),
                  paste0('EP_', net_pd$source_measurement_type, '_', net_pd$target_measurement_type))
p90 <- tapply(abs(net_pd$estimate), stratum,
              function(x) unname(quantile(x, probs = estimate_percentile)))

strata    <- names(thresholds$abs_estimate_thresholds)
p90_table <- data.frame(
  stratum    = strata,
  n_edges    = as.integer(table(stratum)[strata]),
  stored     = unname(thresholds$abs_estimate_thresholds),
  recomputed = round(unname(p90[strata]), 3)
)
write.table(p90_table, file.path(out_folder, 'abs_estimate_p90_recomputed.tsv'),
            sep = '\t', quote = FALSE, row.names = FALSE)
print(p90_table, row.names = FALSE)

differs <- is.na(p90_table$recomputed) |
           abs(p90_table$stored - p90_table$recomputed) > 1e-9
if (any(differs)) {
  warning('Recomputed |estimate| thresholds differ from network_thresholds.txt for: ',
          paste(p90_table$stratum[differs], collapse = ', '),
          '. The filtered network below uses the stored values.')
}

#### Filter and save ####

filtered <- filter_network(network, thresholds)
saveRDS(filtered, file.path(out_folder, 'lasso_network_filtered.rds'))

counts <- table(filtered$link_type)
message('Saved: lasso_network_filtered.rds (',
        paste(names(counts), counts, sep = ' = ', collapse = ', '),
        '; total ', nrow(filtered), ')')
