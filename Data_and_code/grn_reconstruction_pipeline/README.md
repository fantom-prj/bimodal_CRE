# GRN reconstruction pipeline

Reconstructs the gene regulatory network (GRN) itself, from raw TOBIAS
footprinting / Hi-C / paired scRNA-seq-snATAC-seq data through to the final
elastic-net-pruned network. Every other `Data_and_code/` subfolder, and
`Fig5/fig5a_sankey.R`, reads this pipeline's output (via `primary_data_folder`
in `../config.R`) rather than the raw data directly — this is the one place
in the repository where that output is actually produced.

Run the five scripts in order, from inside this folder:

1. **`1_tf_2_peak_set_threshold.R`** — calibrates a TOBIAS footprinting-score
   threshold against ChIP-seq ground truth for two validated transcription
   factors (TEAD4, ONECUT2), by scanning a grid of score-quantile thresholds
   and reporting precision/recall/F1 at each. Informs (but does not directly
   set) the quantile used in step 2.
2. **`2_tf_2_peak.R`** — calls TF-to-peak edges genome-wide from TOBIAS
   BINDetect footprinting scores, for every TF with a motif in the
   reference set. Two modes: footprinting-informed (a per-TF score quantile
   threshold) or motif-only (every peak with a matching motif is called
   bound). The published network uses the motif-only variant.
3. **`3_preparing_base_network.R`** — combines the TF-to-peak calls with
   genomic-proximity and Hi-C chromatin-contact evidence into a single
   candidate edge table (TF->Promoter, TF->Enhancer, Enhancer->Promoter),
   before any expression/accessibility-based pruning. Two control-panel
   flags (`use_hic`, `use_footprinting`) each toggle one evidence source;
   run all four combinations to reproduce the network-variant comparisons in
   `../network_basic_stats/` and `../chipseq_eqtl_validation/`.
4. **`4_preparing_dataset.R`** — restricts the paired scRNA-seq / snATAC-seq
   metacell object to exactly the measurements the base network needs as
   predictors or targets, normalises each modality, and assembles them into
   one matrix ordered by pseudotime.
5. **`5_lasso_models.R`** — the pruning step: fits an elastic-net regression
   of each target's expression/accessibility against its candidate
   regulators, keeping only regulators with a non-zero coefficient as edges
   in the final network. Each retained edge carries a permutation-test
   p-value, a deviance-ratio model-fit-quality score, and a post-hoc OLS
   delta-R^2 — the three quantities every downstream figure filters the
   network on (see `../network_loading/`). By far the most computationally
   expensive step; expect substantial wall time even with heavy
   parallelisation.

## One-time setup

Set `raw_data_folder` in `../config.R` to your copy of the raw pipeline
inputs (TOBIAS BINDetect output, Hi-C contact tables, the metacell object,
ChIP-seq ground truth, motif and pseudotime reference tables — see each
script's header for its exact file list). `primary_data_folder` in the same
file already defaults to `grn_reconstruction_pipeline/output/` here, so no
further configuration is needed before running these five scripts — their
output lands exactly where every other `Data_and_code/` script expects to find
it.

## Already included: the final, published network

`output/5_lasso_models/lasso_network_filtered.rds` ships pre-built in this
repository (127,938 edges) — it is `lasso_network_with_randomization.rds`
(step 5's raw output) after the three confidence filters every figure in
this paper applies (adjpvalue <= 0.05, dev_ratio >= 0.3, per-stratum
|estimate| >= 90th percentile; see `../network_loading/`). If you only want
to explore the final GRN or rerun a downstream figure script, you do not
need `raw_data_folder` or any of the five scripts above — this file already
has everything they would produce, filtered down to the edges the paper
reports.

## Producing every network variant used downstream

Steps 3-5 each have control-panel flags governing which network variant they
build. To reproduce the full 8-network comparison used in
`../network_basic_stats/` and `../chipseq_eqtl_validation/` (4 base variants
x un-pruned/lasso-pruned), rerun steps 3-4 for each combination of
`use_hic`/`use_footprinting`, then rerun step 5 for each of those four
`(base_network, dataset)` pairs — additionally toggling `use_enhancer_tCRE`
on the published (`use_hic = TRUE, use_footprinting = TRUE`) combination
produces the `no_enhancer_tCRE` variant used in `../tCRE_topology/`.
