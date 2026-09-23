# GRN reconstruction and shared data regeneration for Fig5, Fig6, ExtFig7

This document covers the subset of `Data_and_code/` subfolders behind
Figures 5-6 and Extended Data Figure 7: `grn_reconstruction_pipeline/`,
`tf_effects_clustering/`, `threshold_sensitivity_sweep/`,
`chipseq_eqtl_validation/`, `network_basic_stats/`, `tCRE_topology/`,
`enhancer_centrality_stats/`, `screen_correlation/`, `ep_enrichment/`,
`tf_tf_umap/`, `abc_validation/`, and `network_loading/` — plus `config.R`
at the top of this folder. None of these subfolders hold any figure panels
themselves:

- `grn_reconstruction_pipeline/` — the pipeline that builds the gene
  regulatory network (GRN) itself, from raw data through to the final
  elastic-net-pruned network; every other subfolder listed above reads its
  output.
- Everything else listed above: scripts that regenerate the small,
  pre-extracted data files consumed by the panel scripts in `../Fig5/`,
  `../Fig6/`, and `../ExtFig7/` — anything computed from the GRN pipeline's
  output that either (a) is shared by more than one panel, or (b) is cheap to
  isolate into its own regeneration step for traceability. One panel
  (`Fig5/fig5a_sankey.R`) is a documented exception that reads the GRN
  pipeline's output directly instead — see the note at the bottom.

## One-time setup

Edit `config.R` in this folder: set `raw_data_folder` to your local copy of
the raw pipeline inputs (see `grn_reconstruction_pipeline/README.md` for the
exact file list). `primary_data_folder` — the GRN pipeline's own output —
already defaults to `grn_reconstruction_pipeline/output/`, so it needs no
edit unless you already have those files elsewhere. Every regeneration
script below sources this one file — you only need to set it once.

## Running order

Run `grn_reconstruction_pipeline/`'s five scripts first (in numeric order),
to produce the GRN itself. After that, there is no cross-dependency between
the remaining subfolders — each can be run independently, in any order. Run
a subfolder's `build_*.R` script (from inside that subfolder, so its own
relative paths resolve) before running the panel script(s) listed as its
consumers.

| Subfolder | Regenerates | Feeds panel(s) |
|---|---|---|
| `grn_reconstruction_pipeline/` | The gene regulatory network itself: TF-to-peak calls, the candidate base network, the expression/accessibility dataset, and the final lasso-pruned network(s) | (all subfolders below, which read its output as `primary_data_folder`) |
| `tf_effects_clustering/` | TF cumulative-effect matrix + heatmap clustering (`tf_effects_clustering.rds`) | Fig5b, Fig5c |
| `threshold_sensitivity_sweep/` | 13-grid-point threshold sweep, both themes (`sweep_category_table.csv`) + shared plotting function | Fig5d, Ext5e/f (combined, one script/image), Fig6a, Ext6a, Ext6b |
| `chipseq_eqtl_validation/` | ChIP-seq F1 curve data + eQTL F1 curve data | Ext5a, Ext5b |
| `network_basic_stats/` | Node/edge counts across all 8 network variants | Ext5c (both node and edge panels) |
| `abc_validation/` | Network-tier vs ABC concordance F1 table (`tier_comparison.csv`) | Ext5d |
| `tCRE_topology/` | New-route triplet table + distance-binning tables (with/without tCRE) | Fig6b, Fig6c, Ext6c, Ext6d |
| `enhancer_centrality_stats/` | Copy of the pipeline's `all_enhancers_stats.csv` + screen overlap tables | Fig6d, Fig6e, Ext6e, Ext6f |
| `screen_correlation/` | Screen hit-score correlation, GSEA, and hit-score-sign test tables | Fig6e, Ext6e, Ext6f |
| `ep_enrichment/` | Promoter-collapsed epigenomic feature enrichment table | Ext7a |
| `tf_tf_umap/` | TF coefficient/association matrices + fixed UMAP layout + dominance classification | Ext7b |
| `network_loading/` | Shared helper (`load_and_filter_network.R`) + canonical `network_thresholds.txt` — sourced by most of the above, not a data-producer itself | (infrastructure only) |

Each subfolder's own `README.md` documents exactly what analysis it performs
and the exact output file(s) it writes.

## Documented exception: Fig5a

`Fig5/fig5a_sankey.R` reads `lasso_network_with_randomization.rds` and
`dataset.rds` directly (from `primary_data_folder`) and is not fed by any
`Data_and_code/` subfolder — its custom Sankey ribbon geometry needs dozens of
fine-grained per-modality-pair edge counts that don't decompose into a
reusable summary table, so splitting counting from drawing would have been a
leaky abstraction. It still sources `Data_and_code/config.R` for the data-root
path and `Data_and_code/network_loading/plot_theme_and_save.R` for the save
helpers.

## Composing the final figures

Each figure's panels are saved as standalone PNG/PDF files in that figure's
own `out/` folder, styled consistently via the shared `bimodal_theme`
(Arial, 5-7pt). There is no in-repo assembly script that composites a
figure's panels into one page — multi-panel figures are hand-composed later
outside R, from these individually-polished panel files.
