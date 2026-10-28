# Network-tier vs ABC concordance summary

Supplies the tier-comparison table read by one panel: ExtFig5d.

Independently validates the GRN's Enhancer_promoter edges against ABC
(Activity-by-Contact), an orthogonal enhancer-gene regulatory map built from
chromatin accessibility, H3K27ac and Hi-C contact frequency (see
`../Fig4.ABC/`), asking whether elastic-net selection and the standard
confidence filters increase agreement with this independent method relative
to a no-selection baseline.

Run in this order, from inside this folder:

1. `build_tier_comparison.R` — compares three network tiers (`base`: every
   proximity/Hi-C-eligible candidate pair before model fitting;
   `elastic_net`: edges kept by the elastic-net model before the standard
   filters; `filtered`: the final published high-confidence GRN) against
   four ABC score variants (three cell types plus their per-pair maximum).
   Each comparison is restricted to the candidate pairs ABC actually
   evaluated for that variant: ABC's full score matrix has a row (with the
   real score, including 0) for every pair it evaluated and none for pairs
   it did not, so "not evaluated" is never counted as an ABC negative. ABC
   peaks that are promoter-like in the aCRE catalogue are removed first. A
   pair is ABC-positive if its score exceeds 0.02, ABC's own reporting
   threshold. Reads `ABCscore_matrix_3cells.tsv.gz` and
   `all_aCRE.PE.final.tsv` via `raw_data_folder`, and the base, lasso and
   filtered networks via `primary_data_folder`. Writes
   `[primary_data_folder]/19_abc_validation_full_matrix/tier_comparison.csv`.
   Memory-heavy: the matrix has ~24 million rows (expect a peak of tens of
   GB of RAM).
2. `copy_tier_comparison.R` — performs no computation; copies the file into
   `Fig5/data/` so the panel script reads one documented location.

Columns used downstream: `tier` (`base` / `elastic_net` / `filtered`),
`abc_variant` (`iPSC` / `NSC` / `Neuron` / `max_across_celltypes`), `f1` (F1
between network-edge membership and ABC positivity over the pairs ABC
evaluated for that variant). The file also carries Jaccard, precision,
recall, Fisher odds ratio, AUC, and Spearman/Pearson correlation columns,
not used by this panel.

Cross-check (published network): 3 tiers, 4 ABC variants (12 rows total).
F1 for `max_across_celltypes` rises monotonically across tiers: 0.084
(`base`) -> 0.124 (`elastic_net`) -> 0.170 (`filtered`).
