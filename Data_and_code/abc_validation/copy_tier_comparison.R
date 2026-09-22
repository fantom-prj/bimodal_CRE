#### Regeneration script: network-tier vs ABC concordance summary (provenance) ####
#
# Supplies the tier-comparison table used by ONE panel: ExtFig5d.
#
# Independently validates the GRN's Enhancer_promoter edges against ABC
# (Activity-by-Contact), an orthogonal enhancer-gene regulatory map built
# from chromatin accessibility and Hi-C contact frequency. Three network
# tiers are compared, each restricted to the (enhancer, gene) candidate
# pairs ABC actually evaluated for the relevant comparison:
#   base         every proximity/Hi-C-eligible candidate pair, before any
#                model fitting — the no-selection floor.
#   elastic_net  edges kept by the elastic-net model, before the standard
#                significance/effect-size filters.
#   filtered     the final published high-confidence GRN (selection + the
#                standard filters).
# For each tier, F1 is computed between "network calls this pair an edge"
# and "ABC score exceeds its own reporting threshold (0.02)", separately
# for four ABC score variants (three cell types plus their per-pair maximum)
# — quantifying whether elastic-net selection and the standard confidence
# filters increase agreement with this independent method, versus a
# no-selection baseline.
#
# Unlike most Data_and_code/ regeneration scripts, this one performs no
# computation: `tier_comparison.csv` is already a finished output of the GRN
# reconstruction pipeline's ABC validation analysis. Regenerating it from
# scratch means re-running that analysis against ABC's full, unfiltered
# genome-wide score matrix (tens of millions of candidate pairs). This
# script exists to record that provenance and to copy the file into
# Fig5/data/ so the panel script reads one documented location.
#
# `tier_comparison.csv` columns used downstream:
#   tier          base | elastic_net | filtered
#   abc_variant   iPSC | NSC | Neuron | max_across_celltypes
#   f1            F1 between network-edge membership and ABC positivity,
#                 restricted to the pairs ABC evaluated for that variant
# (also carries precision, recall, Fisher OR, AUC, and Spearman/Pearson
# correlation columns, not used by this panel)
#
# Output (into Fig5/data/): tier_comparison.csv

rm(list = ls())

# ── CONTROL PANEL ─────────────────────────────────────────────────────────────
source('../config.R')  # defines primary_data_folder

tier_comparison_file <- file.path(primary_data_folder,
                                  '19_abc_validation_full_matrix', 'tier_comparison.csv')
out_folder <- '../../Fig5/data'
# ──────────────────────────────────────────────────────────────────────────────

dir.create(out_folder, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(tier_comparison_file)) stop('Missing input: ', tier_comparison_file)
file.copy(tier_comparison_file, file.path(out_folder, 'tier_comparison.csv'), overwrite = TRUE)
message('Copied: tier_comparison.csv')

tc <- read.csv(file.path(out_folder, 'tier_comparison.csv'), stringsAsFactors = FALSE)
message('Tiers: ', paste(unique(tc$tier), collapse = ', '))
message('ABC variants: ', paste(unique(tc$abc_variant), collapse = ', '))
message('Done.')
