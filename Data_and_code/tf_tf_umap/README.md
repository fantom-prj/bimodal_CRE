# TF-TF regulatory UMAP with dominance classification

Regenerates the fixed-layout UMAP and direct/indirect edge tables used by
ExtFig7b.

* `build_tf_umap_data.R` — builds a fixed 2D layout for every TF in the
  network (network loading/filtering by adjpvalue, dev_ratio, per-stratum
  |estimate| p90, positive-only EP edges; a TF co-expression matrix from
  pseudotime-smoothed expression; the three TF-TF regulatory coefficient
  matrices — direct, indirect, and combined; a combined association matrix
  reduced via PCA then UMAP), then classifies each TF by whether its
  outgoing regulatory influence on other TFs runs mostly through direct
  edges (TF -> promoter(TF)) or indirect ones (TF -> enhancer -> promoter(TF)).

Reads the full lasso network + dataset (via `primary_data_folder`), a
pseudotime table, and a stage-marker TF annotation table (both via
`raw_data_folder` — see `../config.R`).

Writes `ExtFig7/data/tf_umap_dominance.rds`: layout_all, edge_direct,
edge_indirect, node_names, umap_fixed, num_label_per_cat, dom_colours.

Cross-check (this run): 327 TFs, 115,634 filtered edges, 1,084 direct /
838 indirect / 1,824 combined non-zero TF-TF coefficient entries; dominance
breakdown Direct-dominant=133, Indirect-dominant=24, Mixed=106.
