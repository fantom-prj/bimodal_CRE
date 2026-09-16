# Threshold sensitivity sweep (TF-set coupling)

Regenerates the per-node Jaccard category table, and holds the shared
plotting logic, for the six stacked-bar panels: Fig5d, ExtFig5d, ExtFig5e,
Fig6a, ExtFig6a, ExtFig6b.

* `build_sweep_data.R` — tests how robust the TF-set decoupling result is to
  the network's significance thresholds. `run_decoupling()` computes, for
  each node (an expressed enhancer regulated at both chromatin and
  expression level = `TF_enhancer`; or a promoter regulated at both levels =
  `TF_promoter`), the Jaccard index between its two governing TF sets, then
  `categorize_jaccard()` buckets nodes into fully uncoupled / intermediate /
  fully coupled. `run_grid_point()` repeats this across three one-at-a-time
  sweeps — adjpvalue 0.02/0.05*/0.1/0.2, dev_ratio 0.1/0.2/0.3*/0.4/0.5,
  |estimate| percentile 75/80/90*/95 (* = baseline), 13 grid points total.
  Reads the lasso network, dataset and base network via a
  `primary_data_folder` placeholder and
  `../network_loading/network_thresholds.txt` (baseline adjpvalue/dev_ratio
  only).

  Two behaviours are load-bearing and must not be "simplified":
  filter order is adjpvalue -> dev_ratio -> |estimate| percentile (the
  percentile is computed on edges that already passed both other filters,
  matching how the stored p90 values were derived); and the percentile is
  recomputed live per stratum at every grid point including 90, rather than
  read from the thresholds file, because the percentile is the swept axis.

  Writes `sweep_category_table.csv` (78 rows = 13 grid points x 2 themes x 3
  Jaccard categories; columns sweep, sweep_value, theme, category, pct, n,
  n_nodes) into **both** `Fig5/data/` and `Fig6/data/`, since the six panels
  that read it live in two figure folders.

* `plot_stacked_bar_height.R` — draws the linear-height stacked bar (total
  bar height = raw n_nodes, so shrinking sample size at stringent thresholds
  is visible in the figure itself). All six panel scripts `source()` this
  and call `make_stacked_bar_linheight(category_table, theme_name,
  sweep_name)`; they differ only in those two arguments.

Cross-check printed by the script (e.g. estimate_percentile=90, TF_promoter:
fully_uncoupled 61.17%, fully_coupled 10.50%, n_nodes 2457).
