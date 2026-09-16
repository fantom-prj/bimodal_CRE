# Per-enhancer centrality statistics (provenance)

Supplies the per-enhancer centrality table and the two screen-overlap tables
read by four Fig6 panels: Fig6d, Fig6e, ExtFig6e and ExtFig6f.

* `copy_enhancer_stats.R` — performs no computation. `all_enhancers_stats.csv`
  is already a finished output of the GRN reconstruction pipeline's
  enhancer-centrality analysis (control panel: `positive_only = FALSE`, i.e.
  the `all_edges` analysis; per-stratum |estimate| p90 filter applied).
  Regenerating it from scratch means rerunning that analysis against the full
  lasso network. This script records that provenance and copies the file —
  plus `screen_overlap_enhancers.csv` and `screen_overlap_promoters.csv`,
  recording how CRISPRi neurogenesis screen peaks overlap the network's
  enhancer and promoter peaks — into `Fig6/data/` so the panel scripts read
  one documented location. Reads them via a `primary_data_folder` placeholder.

Columns used downstream: `enhancer`; `triplet_pairs` (centrality axis 1 —
unique (TF, promoter) pairs bridged via the enhancer); `cum_abs_estimate_total`
(axis 2 — cumulative |estimate| over all edges touching the enhancer, in +
out); and the `is_expressed` / `is_mediating` / `n_regulating_TFs` /
`n_regulated_genes` / `n_regulated_promoters` annotations. In the
screen-overlap tables each row is one (screen peak x network peak x contrast)
triple, with `hit_score` NA for screen peaks that have no quantified score in
the screen results.

Cross-check printed by the script: 52,334 enhancers in the `all_edges`
analysis, of which 9,268 expressed and 10,813 mediating.
