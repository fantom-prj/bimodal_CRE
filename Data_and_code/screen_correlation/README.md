# CRISPRi screen x enhancer-centrality analyses

Regenerates the statistics behind three panels: Fig6e (Spearman correlation
heatmap), ExtFig6e (GSEA enrichment curves) and ExtFig6f (hit-score sign vs
triplet participation).

* `build_screen_analyses_data.R` — relates each enhancer's two network-
  centrality measures (triplet-pair count and cumulative |estimate|) to its
  CRISPRi neurogenesis screen behaviour, for the `all_edges` analysis:
  Part 1 GSEA of centrality rankings against the expressed-enhancer and
  mediating-enhancer gene sets; Part 2a Spearman correlation of mean screen
  hit score against each centrality axis per differentiation-step contrast;
  Part 3 Wilcoxon/Fisher tests of centrality against the sign of the mean hit
  score.

  Three behaviours are deliberately built in: the `mediating_TRUE` x
  `triplet_pairs` GSEA run is **skipped** (trivially determined by set
  membership, since non-mediating enhancers have `triplet_pairs = 0`) and
  excluded from the BH pool, leaving 3 meaningful runs; `contrast_order` is
  the three adjacent differentiation-step transitions only (iPSvsNSC_P0,
  NSC_P0vsNSC_P2, NSC_P2vsNeuron), which also fixes the FDR pool sizes; and
  Part 3's Wilcoxon and Fisher p-values are BH-corrected as two separate
  families of 3, never pooled into one family of 6.

  The published GSEA curves are not `fgsea::plotEnrichment()` output — the
  running enrichment score (a step up at each gene-set member, down
  elsewhere) is not a summary statistic, so it is computed directly here and
  written out (subsampled to ~4000 points per panel; the rug hit positions
  are written in full).

Reads `all_enhancers_stats.csv` and `screen_overlap_enhancers.csv` from
`Fig6/data/`, where `../enhancer_centrality_stats/copy_enhancer_stats.R`
placed them (that script documents their provenance). Needs no large RDS
file, so it has no `primary_data_folder` placeholder.

Writes into `Fig6/data/`: `gsea_results.csv`, `gsea_running_curves.csv`,
`gsea_hit_positions.csv`, `screen_correlation.csv`,
`hitscore_sign_tests.csv`, `hitscore_sign_enhancers.csv`,
`bimodality_diagnostics.csv`.

Cross-check printed by the script: 145 scored enhancers; iPSvsNSC_P0 x
triplet_pairs rho = -0.235, FDR = 0.024.
