# TF cumulative-effect matrix + heatmap clustering

Regenerates the TF-level cumulative-effect summaries and heatmap clustering
shared by Fig5b (heatmap) and Fig5c (TF-characteristics dot plot).

* `build_tf_effects_data.R` — for every TF, sums its regulatory effect
  (signed model coefficient) across all edges into each of several target
  classes (per-target-type cumulative-effect summaries, `results_bar`), then
  assembles a 4-column matrix (mediating-enhancer aCRE/tCRE, promoter
  aCRE/tCRE) and hierarchically clusters it into activator-like vs
  repressor-like TF groups (`h_wide`, `sh`, `sh_sqrt`, `rk`, `row_group`,
  `row_group_fine`, `col_order`, `cdend`). Reads
  `lasso_network_with_randomization.rds` via a `primary_data_folder`
  placeholder and `../network_loading/network_thresholds.txt`.

Writes `Fig5/data/tf_effects_clustering.rds` (a single list with all of the
above fields plus `results_bar`, `num_tf_to_show`, `stat_label`).

Fixed control-panel choice: `use_abs_estimate = FALSE` (signed cumulative
effect), `apply_estimate_threshold = TRUE` — the variant used in the
published figures.
