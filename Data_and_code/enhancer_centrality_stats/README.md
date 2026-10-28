# Per-enhancer centrality statistics and screen overlap

Supplies the per-enhancer centrality table and the two screen-overlap tables
read by four Fig6 panels: Fig6d, Fig6e, ExtFig6e and ExtFig6f (and by
`../screen_correlation/`).

Run in this order, from inside this folder:

1. `build_enhancer_stats.R` — for every enhancer in the filtered GRN
   (adjpvalue, dev_ratio, per-stratum |estimate| p90), computes two
   centrality axes — `triplet_pairs` (unique (TF, promoter) pairs bridged
   via the enhancer) and `cum_abs_estimate_total` (cumulative |estimate|
   over all edges touching it, in + out) — plus per-enhancer annotations and
   the first Pareto front on the two axes. Control panel `positive_only =
   FALSE` (the `all_edges` analysis, used by the figures); `TRUE` keeps
   activating edges only. Reads the lasso network and `dataset.rds` via
   `primary_data_folder`. Writes
   `[primary_data_folder]/10_enhancer_centrality/all_edges/all_enhancers_stats.csv`.
2. `build_screen_overlap.R` — overlaps the CRISPRi neurogenesis screen's
   peak library with the network's promoter and enhancer peaks by genomic
   coordinates (full, unfiltered network: the overlap is purely positional).
   Reads `all_promoter_coordinate.bed` (screen peak library) and
   `overall.promoter.base.result.osc` (screen results) via `raw_data_folder`.
   Writes `screen_overlap_enhancers.csv` and `screen_overlap_promoters.csv`
   into `[primary_data_folder]/10_enhancer_centrality/`.
3. `copy_enhancer_stats.R` — performs no computation; copies the three
   files into `Fig6/data/` so the panel scripts read one documented location.

Columns used downstream: `enhancer`; `triplet_pairs`; `cum_abs_estimate_total`;
and the `is_expressed` (the enhancer has a tCRE measurement) / `is_mediating`
(on >= 1 complete TF -> enhancer -> promoter path) / `n_regulating_TFs` /
`n_regulated_genes` / `n_regulated_promoters` annotations. In the
screen-overlap tables each row is one (screen peak x network peak x contrast)
triple, with `hit_score` NA for screen peaks that have no quantified score in
the screen results.

Cross-check (published network): 52,334 enhancers in the `all_edges`
analysis, of which 9,268 expressed and 10,813 mediating; 16,615 promoter and
5,280 enhancer screen-overlap rows.

If `build_screen_overlap.R` stops in `findOverlaps` with
`build_NCList: memory allocation failed` — even on a handful of ranges —
the cause is not memory: the installed IRanges / S4Vectors / GenomicRanges
come from different Bioconductor releases. Reinstall them from a single
release (`BiocManager::valid()` lists the mismatched packages).
