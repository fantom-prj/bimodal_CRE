# Enhancer-promoter epigenomic feature enrichment (collapsed)

Regenerates the promoter-collapsed odds-ratio/p-value enrichment table used
by ExtFig7a.

* `build_ep_enrichment_data.R` — tests whether pairs of epigenomic features
  (H3K27ac, H3K27me3, stage-specific super-enhancer membership, CpG islands,
  TATA boxes) co-occur more or less often than expected between an enhancer
  and the promoter(s) it regulates in the GRN, across 7 enhancer features x 7
  promoter features x 3 edge strata (all edges / expressed-enhancer edges /
  mediating-enhancer edges), FDR-corrected within stratum. Mediating-enhancer
  status (an enhancer that is itself under TF control) is resolved directly
  from the network's edge structure.

Reads the full lasso network + dataset (via `primary_data_folder`) and the
peak epigenomic feature table `all_aCRE.PE.final.select.tsv` (via
`raw_data_folder` — see `../config.R`).

Writes `ExtFig7/data/ep_enrichment_collapsed.csv` (147 rows = 7x7x3).
