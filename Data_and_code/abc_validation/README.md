# Network-tier vs ABC concordance summary (provenance)

Supplies the tier-comparison table read by one panel: ExtFig5d.

Independently validates the GRN's Enhancer_promoter edges against ABC
(Activity-by-Contact), an orthogonal enhancer-gene regulatory map built from
chromatin accessibility and Hi-C contact frequency, asking whether
elastic-net selection and the standard confidence filters increase agreement
with this independent method relative to a no-selection baseline.

* `copy_tier_comparison.R` — performs no computation. `tier_comparison.csv`
  is already a finished output of the GRN reconstruction pipeline's ABC
  validation analysis, which compares three network tiers (`base`: every
  proximity/Hi-C-eligible candidate pair before model fitting;
  `elastic_net`: edges kept by the elastic-net model before the standard
  filters; `filtered`: the final published high-confidence GRN) against four
  ABC score variants (three cell types plus their per-pair maximum), each
  restricted to the candidate pairs ABC actually evaluated for that variant.
  Regenerating it from scratch means re-running that analysis against ABC's
  full, unfiltered genome-wide score matrix (tens of millions of candidate
  pairs). This script records that provenance and copies the file into
  `Fig5/data/` so the panel script reads one documented location. Reads it
  via a `primary_data_folder` placeholder.

Columns used downstream: `tier` (`base` / `elastic_net` / `filtered`),
`abc_variant` (`iPSC` / `NSC` / `Neuron` / `max_across_celltypes`), `f1` (F1
between network-edge membership and ABC positivity — ABC score above its own
0.02 reporting threshold — over the pairs ABC evaluated for that variant).
The file also carries precision, recall, Fisher odds ratio, AUC, and
Spearman/Pearson correlation columns, not used by this panel.

Cross-check printed by the script: 3 tiers, 4 ABC variants (12 rows total).
F1 for `max_across_celltypes` rises monotonically across tiers: 0.084
(`base`) -> 0.124 (`elastic_net`) -> 0.170 (`filtered`).
