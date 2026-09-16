#### Regeneration script: per-enhancer centrality statistics (provenance) ####
#
# Supplies the per-enhancer centrality table used by FOUR panels:
#   Fig6d     screen-coverage scatter
#   Fig6e     screen correlation heatmap
#   ExtFig6e  GSEA enrichment plots
#   ExtFig6f  hit-score sign vs triplet participation bar
#
# Unlike the other Data_and_code/ regeneration scripts, this one performs no
# computation: `all_enhancers_stats.csv` is already a finished output of the
# GRN reconstruction pipeline's enhancer-centrality analysis (control panel:
# positive_only = FALSE, i.e. the `all_edges` analysis; per-stratum
# |estimate| p90 filter applied). Regenerating it from scratch means rerunning
# that analysis against the full lasso network. This script exists to record
# that provenance and to copy the file (plus the two screen-overlap tables it
# is analysed against) into Fig6/data/ so the panel scripts read one
# documented location.
#
# `all_enhancers_stats.csv` columns used downstream:
#   enhancer                 peak id
#   triplet_pairs            centrality axis 1: unique (TF, promoter) pairs
#                            bridged via this enhancer
#   cum_abs_estimate_total   centrality axis 2: cumulative |estimate| over all
#                            edges touching the enhancer (in + out)
#   is_expressed, is_mediating, n_regulating_TFs, n_regulated_genes,
#   n_regulated_promoters    per-enhancer annotations
# 52,334 enhancers in the all_edges analysis.
#
# The two screen-overlap tables record how each CRISPRi neurogenesis screen
# peak overlaps the network's enhancer and promoter peaks. Each row is one
# (screen peak x network peak x contrast) triple; `hit_score` is NA for
# screen peaks with no quantified score in the screen results.
#
# Outputs (into Fig6/data/):
#   all_enhancers_stats.csv
#   screen_overlap_enhancers.csv
#   screen_overlap_promoters.csv

rm(list = ls())

# ── CONTROL PANEL ─────────────────────────────────────────────────────────────
source('../config.R')  # defines primary_data_folder

analysis <- 'all_edges'   # 'all_edges' or 'positive_only'; the panels use all_edges

centrality_folder <- file.path(primary_data_folder, '10_enhancer_centrality')
stats_file        <- file.path(centrality_folder, analysis, 'all_enhancers_stats.csv')
enh_overlap_file  <- file.path(centrality_folder, 'screen_overlap_enhancers.csv')
prom_overlap_file <- file.path(centrality_folder, 'screen_overlap_promoters.csv')

out_folder <- '../../Fig6/data'
# ──────────────────────────────────────────────────────────────────────────────

dir.create(out_folder, showWarnings = FALSE, recursive = TRUE)

for (f in c(stats_file, enh_overlap_file, prom_overlap_file)) {
  if (!file.exists(f)) stop('Missing input: ', f)
  file.copy(f, file.path(out_folder, basename(f)), overwrite = TRUE)
  message('Copied: ', basename(f))
}

stats <- read.csv(file.path(out_folder, 'all_enhancers_stats.csv'),
                  stringsAsFactors = FALSE)
message('Enhancers in ', analysis, ': ', format(nrow(stats), big.mark = ','))
message('  expressed: ', sum(stats$is_expressed),
        '  mediating: ', sum(stats$is_mediating))
message('Done.')
