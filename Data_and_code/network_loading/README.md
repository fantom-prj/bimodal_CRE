# Shared network loading and filtering

Common helper sourced by every other `Data_and_code/` regeneration script that
needs the lasso GRN. Not a figure-specific analysis on its own.

* `load_and_filter_network.R` -> defines:
  - `read_thresholds(thresholds_file)` — parses `network_thresholds.txt`
    into `adjpvalue_threshold`, `dev_ratio_threshold`, and the named vector
    `abs_estimate_thresholds` (per-stratum |estimate| p90).
  - `filter_network(net, thresholds, apply_estimate_threshold = TRUE)` —
    applies, in order: `adjpvalue <= adjpvalue_threshold`,
    `dev_ratio >= dev_ratio_threshold`, then (if requested) the per-stratum
    |estimate| p90 filter. Order matters: the stored p90 values in
    `network_thresholds.txt` were themselves derived AFTER the adjpvalue and
    dev_ratio filters, so the percentile filter must be applied last to
    reproduce them exactly.
* `network_thresholds.txt` -> a copy of the canonical per-stratum threshold
  file (adjpvalue, dev_ratio, and the |estimate| p90 cutoffs used to call an
  edge significant), so every script that filters the network to the
  published edge set uses the same values.

## Inputs (external, not committed — see path note below)

* `lasso_network_with_randomization.rds` — full lasso GRN with permutation
  p-values (~123 MB).
* `dataset.rds` — metacell expression/accessibility matrix used by several
  downstream analyses (~73 MB).
* `base_network.rds` — pre-lasso candidate edge set (~67 MB).

These files are large intermediate outputs of the GRN reconstruction pipeline
(gene regulatory network inference from paired scRNA-seq / snATAC-seq / Hi-C
data) and are not included in this repository because of their size. Every
regeneration script in `Data_and_code/` that needs them reads a
`primary_data_folder` path set at the top of the script (edit it to point at
your local copy of these files — see `Data_and_code/config.R`).
