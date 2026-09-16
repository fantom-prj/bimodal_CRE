#### Data-root configuration ####
#
# Every regeneration script under Data_and_code/, plus Fig5/fig5a_sankey.R,
# sources this file to locate the data behind the GRN reconstruction pipeline
# (gene regulatory network inference from paired scRNA-seq / snATAC-seq /
# Hi-C data across iPSC -> NSC -> Neuron differentiation).

# Folder containing the raw and small reference inputs to the pipeline —
# not produced by any script in this repository, so you need to supply your
# own copy. Expected contents:
#   TOBIAS BINDetect output (per-TF footprinting scores)
#   Hi-C contact tables (iPSC / NSC / Neuron, one RDS per cell type)
#   the paired scRNA-seq / snATAC-seq metacell object (Seurat/Signac)
#   ChIP-seq ground truth for the two validated TFs (TEAD4, ONECUT2)
#   motif and TF-grouping reference tables
#   pseudotime, eQTL, and disease/GWAS reference tables
# See Data_and_code/grn_reconstruction_pipeline/README.md for the exact file
# list each script needs.
raw_data_folder <- '[raw_data_folder]'

# Folder containing the pipeline's own outputs — the GRN itself and every
# intermediate built along the way (TF-to-peak calls, the base network, the
# expression/accessibility dataset, and the final lasso-regularised network).
#
# Defaults to Data_and_code/grn_reconstruction_pipeline/output/ (resolved
# relative to this config file, not to whichever script sourced it), so
# running the scripts in Data_and_code/grn_reconstruction_pipeline/ in numeric
# order, from raw_data_folder, reproduces everything every other script below
# needs — no external results folder required. If you already have these
# files from elsewhere, replace the line below with an absolute path to that
# folder and skip re-running the pipeline scripts.
#   1_tf_2_peak_set_threshold/<TF>.rds
#   2_tf_2_peak/tf2peak.rds
#   3_preparing_base_network/base_network.rds (+ variants)
#   4_preparing_dataset/dataset.rds (+ variants)
#   5_lasso_models/lasso_network_with_randomization.rds (+ variants)
#   10_enhancer_centrality/all_edges/all_enhancers_stats.csv
#   10_enhancer_centrality/screen_overlap_enhancers.csv
#   10_enhancer_centrality/screen_overlap_promoters.csv
#
# One file ships pre-built, so most scripts below work out of the box even
# before running the pipeline:
#   5_lasso_models/lasso_network_filtered.rds — the final, published
#     high-confidence GRN (127,938 edges): lasso_network_with_randomization.rds
#     after adjpvalue <= 0.05, dev_ratio >= 0.3, and the per-stratum
#     |estimate| p90 filter (see network_loading/load_and_filter_network.R).
this_config_dir     <- dirname(sys.frames()[[1]]$ofile)
primary_data_folder <- file.path(this_config_dir, 'grn_reconstruction_pipeline/output')
