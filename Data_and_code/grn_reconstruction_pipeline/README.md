# GRN reconstruction pipeline

Reconstructs the gene regulatory network (GRN) itself, from raw TOBIAS
footprinting / Hi-C / paired scRNA-seq-snATAC-seq data through to the final
elastic-net-pruned network. Every other `Data_and_code/` subfolder, and
`Fig5/fig5a_sankey.R`, reads this pipeline's output (via `primary_data_folder`
in `../config.R`) rather than the raw data directly — this is the one place
in the repository where that output is actually produced.

Run the scripts in order, from inside this folder:

1. **`1_tf_2_peak_set_threshold.R`** — calibrates a TOBIAS footprinting-score
   threshold against ChIP-seq ground truth for two validated transcription
   factors (TEAD4, ONECUT2), by scanning a grid of score-quantile thresholds
   and reporting precision/recall/F1 at each. Informs (but does not directly
   set) the quantile used in step 2.
2. **`2_tf_2_peak.R`** — calls TF-to-peak edges genome-wide from TOBIAS
   BINDetect footprinting scores, for every TF with a motif in the
   reference set. Two modes: footprinting-informed (a per-TF score quantile
   threshold) or motif-only (every peak with a matching motif is called
   bound). The published network uses the footprinting-informed variant;
   the motif-only one feeds the `no_footprinting` network variants.
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
5. **`5a_lasso_models.R`** — the pruning step: fits an elastic-net regression
   of each target's expression/accessibility against its candidate
   regulators, keeping only regulators with a non-zero coefficient as edges
   in the final network. Each retained edge carries a permutation-test
   p-value, a deviance-ratio model-fit-quality score, and a post-hoc OLS
   delta-R^2 — the three quantities every downstream figure filters the
   network on (see `../network_loading/`). By far the most computationally
   expensive step; expect substantial wall time even with heavy
   parallelisation.
6. **`5b_assemble_network_variants.R`** — builds the four network variants
   other than the published one from parts fitted separately by step 5a, and
   first rebuilds any step 5a output that was run in target blocks (see
   "Producing every network variant used downstream" below). No model
   fitting.

## One-time setup

Set `raw_data_folder` in `../config.R` to your copy of the raw pipeline
inputs (TOBIAS BINDetect output, Hi-C contact tables, the metacell object,
ChIP-seq ground truth, motif and pseudotime reference tables — see each
script's header for its exact file list). `primary_data_folder` in the same
file already defaults to `grn_reconstruction_pipeline/output/` here, so no
further configuration is needed before running these scripts — their
output lands exactly where every other `Data_and_code/` script expects to find
it.

## Already included: the final, published network

`output/5_lasso_models/lasso_network_filtered.rds` ships pre-built in this
repository (141,733 edges) — it is `lasso_network_with_randomization.rds`
(step 5a's raw output) after the three confidence filters every figure in
this paper applies (adjpvalue <= 0.05, dev_ratio >= 0.3, per-stratum
|estimate| >= 90th percentile; see `../network_loading/`). If you only want
to explore the final GRN or rerun a downstream figure script, you do not
need `raw_data_folder` or any of the scripts above — this file already
has everything they would produce, filtered down to the edges the paper
reports.

## Producing every network variant used downstream

Steps 3-5a each have control-panel flags governing which network variant they
build. To reproduce the full 8-network comparison used in
`../network_basic_stats/` and `../chipseq_eqtl_validation/` (4 base variants
x un-pruned/lasso-pruned), rerun steps 3-4 for each combination of
`use_hic`/`use_footprinting`. Toggling `use_enhancer_tCRE` in step 5a on the
published (`use_hic = TRUE, use_footprinting = TRUE`) combination produces the
`no_enhancer_tCRE` variant used in `../tCRE_topology/`.

Step 5a does not need a full run for every lasso-pruned variant. It fits two
independent families of models: TF-to-peak models, with TFs as the only
predictors, and enhancer-to-promoter models, with enhancers as the only
predictors. `use_hic` changes only the enhancer-to-promoter candidate edges and
`use_footprinting` only the TF-to-peak ones. `use_enhancer_tCRE = FALSE` drops
the enhancer tCRE target models from the TF-to-peak part and changes the
predictors of the enhancer-to-promoter models. Within a part, every model's
input values, target list, target order and permutation seeds are the same in
all variants that share that part, and p-values are BH-adjusted within each
model, so a shared part is identical whichever variant it was fitted for.
Setting `fit_edge_classes` in step 5a to a single class fits only that part and
writes it as `..._TF_peak_part.rds` or `..._Enhancer_promoter_part.rds`, never
under a full-network name. The variants are then built as follows:

| Network | Step 5a runs needed |
|---|---|
| `lasso_network_with_randomization` (published) | full run (`use_hic = TRUE`, `use_footprinting = TRUE`) |
| `lasso_network_with_randomization_no_enhancer_tCRE` | `fit_edge_classes = 'Enhancer_promoter'`, `use_enhancer_tCRE = FALSE` |
| `lasso_network_no_hic_with_randomization` | `fit_edge_classes = 'Enhancer_promoter'`, `use_hic = FALSE` |
| `lasso_network_no_footprinting_with_randomization` | `fit_edge_classes = 'TF_peak'`, `use_footprinting = FALSE` |
| `lasso_network_no_hic_no_footprinting_with_randomization` | none (TF-to-peak part of no_footprinting + enhancer-to-promoter part of no_hic) |

After these four step 5a runs, `5b_assemble_network_variants.R` writes the four
variants under their standard names. The TF-to-peak part without footprinting
has about five times as many candidate edges as the published one and is the
longest of these runs.

### Splitting step 5a across jobs

On a cluster, any step 5a run can be split into blocks of targets that run as
separate jobs: set `target_blocks` to the number of blocks and `target_block`
to 1, 2, ... in each job. Block b fits targets b, b + B, b + 2B, ..., and
every target keeps the permutation seeds it has in a single run. Each block
writes `<name>_block<b>of<B>.rds`; `5b_assemble_network_variants.R` rebuilds
`<name>.rds` from the complete set of blocks, identical to a single run, and
assembles the variants once all four inputs above exist. Running 5b after
the published network's blocks alone rebuilds just that network.
