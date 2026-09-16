# ChIP-seq and eQTL validation F1 curve data

Two independent regeneration scripts, each validating the inferred GRN edges
against an orthogonal ground truth and producing the F1-curve data for one
Ext Fig 5 panel.

* `build_chipseq_f1_data.R` — for two transcription factors with available
  ChIP-seq data (TEAD4, ONECUT2), checks how well each network variant's
  TF-to-peak edges recover the TF's true ChIP-seq binding sites (direct
  TF->peak or indirect TF->enhancer->peak), computing precision/recall/F1
  across dev_ratio thresholds x 8 network variants. Writes
  `Fig5/data/chipseq_f1_either_all.csv`. Feeds `Fig5/ext_fig5a_chipseq_f1.R`.

* `build_eqtl_f1_data.R` — checks how well each network variant's
  enhancer-to-promoter edges recover known enhancer-to-gene links from
  fine-mapped eQTLs (gene symbols resolved via org.Hs.eg.db), computing
  precision/recall/F1 across dev_ratio thresholds x 8 network variants x 5
  source/target measurement-type pairs, pooled into an 'All interactions'
  slice. Writes `Fig5/data/eqtl_f1_all_interactions.csv`. Feeds
  `Fig5/ext_fig5b_eqtl_f1.R`.

Both read the 8 base+lasso network RDS files via `primary_data_folder`, and a
small external ground-truth table (ChIP-seq or eQTL) via
`raw_data_folder` (see `../config.R`).
