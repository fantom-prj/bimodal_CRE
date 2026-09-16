#### Shared helper: load network_thresholds.txt and filter the lasso GRN ####
#
# Sourced by other Data_and_code/ regeneration scripts. Defines:
#   read_thresholds(thresholds_file)
#   filter_network(net, thresholds, apply_estimate_threshold = TRUE)
#
# Filter order (adjpvalue -> dev_ratio -> |estimate| percentile) matches how
# the stored p90 values in network_thresholds.txt were themselves derived;
# applying the percentile filter before dev_ratio silently inflates sample
# sizes, since it changes which edges the per-stratum quantile is computed
# over.

read_thresholds <- function(thresholds_file) {
  thr_lines  <- readLines(thresholds_file)
  thr_lines  <- thr_lines[!grepl('^\\s*#', thr_lines) & grepl('=', thr_lines)]
  thr_parsed <- setNames(
    as.numeric(trimws(sub('.*=', '', thr_lines))),
    trimws(sub('=.*', '', thr_lines))
  )
  list(
    adjpvalue_threshold = unname(thr_parsed['adjpvalue_threshold']),
    dev_ratio_threshold = unname(thr_parsed['dev_ratio_threshold']),
    abs_estimate_thresholds = c(
      TF_Enhancer_aCRE = unname(thr_parsed['abs_estimate_TF_peak_Enhancer_aCRE']),
      TF_Enhancer_tCRE = unname(thr_parsed['abs_estimate_TF_peak_Enhancer_tCRE']),
      TF_Promoter_aCRE = unname(thr_parsed['abs_estimate_TF_peak_Promoter_aCRE']),
      TF_Promoter_tCRE = unname(thr_parsed['abs_estimate_TF_peak_Promoter_tCRE']),
      EP_aCRE_aCRE      = unname(thr_parsed['abs_estimate_Enhancer_promoter_aCRE_aCRE']),
      EP_aCRE_tCRE      = unname(thr_parsed['abs_estimate_Enhancer_promoter_aCRE_tCRE']),
      EP_tCRE_aCRE      = unname(thr_parsed['abs_estimate_Enhancer_promoter_tCRE_aCRE']),
      EP_tCRE_tCRE      = unname(thr_parsed['abs_estimate_Enhancer_promoter_tCRE_tCRE'])
    )
  )
}

filter_network <- function(net, thresholds, apply_estimate_threshold = TRUE) {
  net <- net[net$adjpvalue <= thresholds$adjpvalue_threshold, ]
  net <- net[net$dev_ratio >= thresholds$dev_ratio_threshold, ]

  if (apply_estimate_threshold) {
    tf_key <- paste0('TF_', net$target_type, '_', net$target_measurement_type)
    tf_thr <- thresholds$abs_estimate_thresholds[tf_key]
    ep_key <- paste0('EP_', net$source_measurement_type, '_', net$target_measurement_type)
    ep_thr <- thresholds$abs_estimate_thresholds[ep_key]
    thr    <- ifelse(net$link_type == 'TF_peak', tf_thr, ep_thr)
    net    <- net[!is.na(thr) & abs(net$estimate) >= thr, ]
  }
  net
}
