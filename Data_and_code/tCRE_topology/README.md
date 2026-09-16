# Enhancer tCRE topology and distance-binned edge counts

Regenerates the shared computation behind four Fig6 panels: Fig6b, Fig6c,
ExtFig6c and ExtFig6d.

* `build_tCRE_topology_data.R` — asks what modelling an enhancer's own
  transcriptional activity (tCRE) adds to the regulatory network beyond
  chromatin accessibility (aCRE) alone. Computes complete TF ->
  enhancer(tCRE) -> promoter triplets and the `new_route` flag (a triplet is
  a new route when its TF/promoter pair is absent from the aCRE-only
  network), plus pool-size-corrected edges-per-enhancer in 50 kb
  enhancer-promoter distance bins, for the with-tCRE network, the
  distinct-promoters double-chance control, and the without-tCRE
  equal-footing control. The plots themselves live in the four panel
  scripts; only the distance-bin summary's data half is computed here.

  The `gene` column on the triplet table is resolved here, via a
  `promoter_gene_lookup` (built from `ep_with$target` / `target_gene`), so
  ExtFig6d does not rebuild the lookup.

  **Filtering is deliberately not the shared `filter_network()` helper.**
  This script applies adjpvalue <= 0.05, then a per-stratum |estimate| p90
  recomputed *live within each network separately*, then dev_ratio >= 0.3.
  The live recompute is intentional: this is the one analysis comparing two
  structurally different networks (with vs without enhancer tCRE
  measurements), where imposing one frozen threshold set derived from the
  with-tCRE network would apply inconsistent fractional stringency. Every
  other regeneration script here reads the frozen p90 values from
  `network_thresholds.txt`; this one is the documented exception.

Reads `lasso_network_with_randomization.rds`,
`lasso_network_with_randomization_no_enhancer_tCRE.rds`, `dataset.rds` and
`base_network.rds` via a `primary_data_folder` placeholder, plus
`../network_loading/network_thresholds.txt`.

Writes into `Fig6/data/`:
`distance_bin_with_tCRE.csv` and `distance_bin_without_tCRE.csv` (columns
source_class, dist_bin, n_enhancers, y_mean), `distance_bin_promoters.csv`
(the double-chance control, not used by a requested panel but part of the
same Theme 5 computation), and `triplets_tCRE.csv` (3,929 rows: enhancer, TF,
promoter, promoter_measurement_type, tf_prom_id, new_route, gene).

Cross-check printed by the script and matching the published figure: 3,929
complete tCRE triplets, 90.3% new route / 9.7% shared, 3,098 new
TF->enhancer(tCRE) edges, 3,390 new enhancer(tCRE)->promoter edges.
