#### Regeneration script: CRISPRi screen x enhancer-centrality analyses ####
#
# Produces the statistics behind THREE panels:
#   Fig6e     screen correlation heatmap (Spearman rho per contrast x axis)
#   ExtFig6e  GSEA enrichment curves (expressed / mediating vs centrality)
#   ExtFig6f  hit-score sign vs triplet participation bar
#
# Relates each enhancer's two network-centrality measures (triplet-pair
# count and cumulative |estimate|) to its CRISPRi neurogenesis screen
# behaviour, testing whether topologically central enhancers -- ones that
# bridge many TF/promoter pairs or carry strong cumulative regulatory
# weight -- are more likely to matter functionally in the screen:
#   Part 1  GSEA of centrality rankings against the expressed-enhancer and
#           mediating-enhancer gene sets                    -> ExtFig6e
#   Part 2a Spearman correlation of mean screen hit score against each
#           centrality axis, per differentiation-step contrast -> Fig6e
#   Part 3  Wilcoxon/Fisher tests of centrality against the SIGN of the mean
#           hit score (does a positive vs negative screen phenotype track
#           with centrality or triplet participation?)       -> ExtFig6f
#   Analysis variant: all_edges.
#
# contrast_order is restricted to the 3 adjacent differentiation-step
# transitions (rather than all pairwise screen contrasts); this also sets the
# FDR pool size (6 tests for Part 2a, 3 per test family for Part 3). Do not
# widen it without re-checking every downstream number.
#
# Part 1 notes:
#   * 2 gene sets (expressed_TRUE, mediating_TRUE) x 2 centrality rankings.
#   * mediating_TRUE x triplet_pairs is SKIPPED: non-mediating enhancers have
#     triplet_pairs = 0 by definition, so the test is trivially determined by
#     set membership (and fgsea takes ~10 min on the degenerate ranking). Its
#     row carries NES/pval/fdr = NA and meaningful = FALSE; BH is applied to
#     the 3 meaningful runs only.
#   * The published enrichment curves are NOT fgsea::plotEnrichment() output --
#     the running enrichment score is computed directly (a +1/n_hit step
#     at each set member, -1/(n_total - n_hit) elsewhere) and drawn in base
#     graphics. That running-sum vector is not a summary statistic, so it is
#     recomputed and written out here, per gene set x axis, alongside the hit
#     positions used for the rug.
#   * set.seed(42) before the fgsea runs, for reproducibility of the
#     permutation-based null distribution.
#
# Part 3 notes:
#   * Enhancers split by sign(mean hit score) >= 0 alone (no mixture-model
#     cutpoint) -- chosen for biological interpretability.
#   * Wilcoxon (cum_abs_estimate_total ~ sign) and Fisher (triplet_pairs > 0 x
#     sign) are BH-corrected as SEPARATE families of 3, not pooled into one
#     family of 6.
#   * Bimodality diagnostics (Hartigan dip test + mclust G=1 vs G=2 BIC) are
#     reported but deliberately NOT used to gate the tests.
#
# Outputs (into Fig6/data/):
#   screen_correlation.csv        analysis, contrast, axis, n, rho, pval, fdr
#   gsea_results.csv              analysis, axis, gene_set, set_size, NES,
#                                 pval, meaningful, fdr
#   gsea_running_curves.csv       axis, gene_set, rank, es  (subsampled to
#                                 keep the file small -- see note below)
#   gsea_hit_positions.csv        axis, gene_set, rank  (rug positions)
#   hitscore_sign_tests.csv       analysis, contrast, test, n, n_pos, n_neg,
#                                 OR, CI_lo, CI_hi, pval, fdr
#   hitscore_sign_enhancers.csv   per (enhancer x contrast): mean_hit_score,
#                                 triplet_pairs, cum_abs_estimate_total,
#                                 has_triplet, hit_sign  (feeds ExtFig6f's bar)
#   bimodality_diagnostics.csv    analysis, contrast, n, dip_stat, dip_pval,
#                                 mclust_best_G, bic_G1, bic_G2

rm(list = ls())
library(tidyverse)
library(fgsea)
library(mclust)
library(diptest)

select <- dplyr::select
filter <- dplyr::filter

# ── CONTROL PANEL ─────────────────────────────────────────────────────────────
# Inputs come from Fig6/data/, placed there by
# ../enhancer_centrality_stats/copy_enhancer_stats.R (which documents their
# provenance).
data_folder <- '../../Fig6/data'
out_folder  <- '../../Fig6/data'

analysis <- 'all_edges'

contrast_order <- c('iPSvsNSC_P0', 'NSC_P0vsNSC_P2', 'NSC_P2vsNeuron')

# Running-curve output is subsampled to this many points per panel; the curves
# are ~52k steps long and plot identically at this resolution. Hit positions
# (the rug) are written in full.
curve_points <- 4000
# ──────────────────────────────────────────────────────────────────────────────

dir.create(out_folder, showWarnings = FALSE, recursive = TRUE)

stats        <- read.csv(file.path(data_folder, 'all_enhancers_stats.csv'),
                         stringsAsFactors = FALSE)
enh_overlap  <- read.csv(file.path(data_folder, 'screen_overlap_enhancers.csv'),
                         stringsAsFactors = FALSE)

scored_ids <- unique(enh_overlap$network_peak_id[!is.na(enh_overlap$hit_score)])
message('Enhancers: ', nrow(stats), '  scored (>=1 hit_score): ', length(scored_ids))

# Per-enhancer mean hit_score per contrast
scored_mean <- enh_overlap %>%
  filter(!is.na(hit_score)) %>%
  group_by(network_peak_id, compare) %>%
  summarise(mean_hit_score = mean(hit_score, na.rm = TRUE), .groups = 'drop')

scored_in_analysis <- stats$enhancer[stats$enhancer %in% scored_ids]
message('Scored enhancers present in this analysis: ', length(scored_in_analysis))

#### PART 1 — GSEA ####

message('Part 1: GSEA...')

gene_sets <- list(
  expressed_TRUE = stats$enhancer[stats$is_expressed == TRUE],
  mediating_TRUE = stats$enhancer[stats$is_mediating == TRUE]
)
axes <- list(
  triplet_pairs          = setNames(stats$triplet_pairs,          stats$enhancer),
  cum_abs_estimate_total = setNames(stats$cum_abs_estimate_total, stats$enhancer)
)

set.seed(42)
all_res <- list()
for (axis_name in names(axes)) {
  ranked <- sort(axes[[axis_name]], decreasing = TRUE)
  for (gs_name in names(gene_sets)) {
    gs <- gene_sets[[gs_name]]
    is_trivial <- (axis_name == 'triplet_pairs' && gs_name == 'mediating_TRUE')
    message('  axis=', axis_name, ' set=', gs_name, ' (n=', length(gs), ')',
            if (is_trivial) ' [SKIPPED - trivial]' else '')

    if (length(gs) == 0 || is_trivial) {
      all_res[[paste(axis_name, gs_name, sep = '__')]] <- data.frame(
        analysis = analysis, axis = axis_name, gene_set = gs_name,
        set_size = length(gs), NES = NA_real_, pval = NA_real_,
        stringsAsFactors = FALSE)
      next
    }
    res <- fgsea(pathways = setNames(list(gs), gs_name), stats = ranked,
                 minSize = 1, maxSize = Inf, nPermSimple = 10000)
    all_res[[paste(axis_name, gs_name, sep = '__')]] <- data.frame(
      analysis = analysis, axis = axis_name, gene_set = gs_name,
      set_size = length(gs), NES = res$NES, pval = res$pval,
      stringsAsFactors = FALSE)
  }
}
gsea_df <- do.call(rbind, all_res)

gsea_df$meaningful <- !(gsea_df$axis == 'triplet_pairs' &
                          gsea_df$gene_set == 'mediating_TRUE') &
                      !is.na(gsea_df$pval)
gsea_df$fdr <- NA_real_
gsea_df$fdr[gsea_df$meaningful] <- p.adjust(gsea_df$pval[gsea_df$meaningful],
                                            method = 'BH')
gsea_df <- gsea_df[order(!gsea_df$meaningful, gsea_df$fdr, gsea_df$pval), ]

write.csv(gsea_df, file.path(out_folder, 'gsea_results.csv'), row.names = FALSE)
message('  Saved: gsea_results.csv (', nrow(gsea_df), ' runs)')

# Running enrichment score per (axis, gene set): a walk down the ranked
# enhancer list that steps up at each gene-set member and down otherwise,
# vectorised via cumsum(); the published curves plot this trajectory rather
# than a single enrichment-score summary.
curve_list <- list(); hit_list <- list()
for (axis_name in names(axes)) {
  ranked <- sort(axes[[axis_name]], decreasing = TRUE)
  for (gs_name in names(gene_sets)) {
    gs      <- gene_sets[[gs_name]]
    n_total <- length(ranked)
    in_set  <- names(ranked) %in% gs
    hits    <- which(in_set)
    n_hit   <- length(hits)
    if (n_hit == 0) next

    step_up   <- 1 / n_hit
    step_down <- 1 / (n_total - n_hit)
    es_vec    <- cumsum(ifelse(in_set, step_up, -step_down))

    keep <- unique(c(1L,
                     round(seq(1, n_total, length.out = min(curve_points, n_total))),
                     n_total))
    curve_list[[paste(axis_name, gs_name)]] <- data.frame(
      axis = axis_name, gene_set = gs_name, rank = keep, es = es_vec[keep])
    hit_list[[paste(axis_name, gs_name)]] <- data.frame(
      axis = axis_name, gene_set = gs_name, rank = hits)
  }
}
write.csv(do.call(rbind, curve_list), row.names = FALSE,
          file.path(out_folder, 'gsea_running_curves.csv'))
write.csv(do.call(rbind, hit_list), row.names = FALSE,
          file.path(out_folder, 'gsea_hit_positions.csv'))
message('  Saved: gsea_running_curves.csv, gsea_hit_positions.csv')

#### PART 2a — hit score x centrality Spearman ####

message('Part 2a: hit-score x centrality correlation...')

scored_hits <- scored_mean %>%
  filter(network_peak_id %in% scored_in_analysis) %>%
  dplyr::rename(enhancer = network_peak_id)

centrality <- stats %>%
  select(enhancer, triplet_pairs, cum_abs_estimate_total) %>%
  filter(enhancer %in% scored_in_analysis)

cor_df <- scored_hits %>%
  left_join(centrality, by = 'enhancer') %>%
  filter(!is.na(mean_hit_score))

cor_results <- list()
for (ctr in contrast_order) {
  sub <- cor_df[cor_df$compare == ctr, ]
  if (nrow(sub) < 3) next
  for (axis in c('triplet_pairs', 'cum_abs_estimate_total')) {
    ct <- cor.test(sub[[axis]], sub$mean_hit_score,
                   method = 'spearman', exact = FALSE)
    cor_results[[paste(ctr, axis, sep = '__')]] <- data.frame(
      analysis = analysis, contrast = ctr, axis = axis, n = nrow(sub),
      rho = ct$estimate, pval = ct$p.value, stringsAsFactors = FALSE)
  }
}
cor_res <- do.call(rbind, cor_results)
cor_res$fdr <- p.adjust(cor_res$pval, method = 'BH')
cor_res <- cor_res[order(cor_res$fdr, cor_res$pval), ]
rownames(cor_res) <- NULL

write.csv(cor_res, file.path(out_folder, 'screen_correlation.csv'), row.names = FALSE)
message('  Saved: screen_correlation.csv (', nrow(cor_res), ' tests; FDR<0.05: ',
        sum(cor_res$fdr < 0.05, na.rm = TRUE), ')')

#### PART 3 — binarized hit-score sign vs centrality ####

message('Part 3: hit-score sign vs centrality...')

centrality3 <- stats %>%
  select(enhancer, triplet_pairs, cum_abs_estimate_total) %>%
  filter(enhancer %in% scored_in_analysis) %>%
  mutate(has_triplet = triplet_pairs > 0)

diag_list <- list(); test_list <- list(); plot_list <- list()

for (ctr in contrast_order) {
  sub <- scored_mean %>%
    filter(compare == ctr, network_peak_id %in% scored_in_analysis) %>%
    dplyr::rename(enhancer = network_peak_id) %>%
    left_join(centrality3, by = 'enhancer')

  if (nrow(sub) < 6) {
    message('  ', ctr, ': too few scored enhancers (n=', nrow(sub), ') - skipped.')
    next
  }

  sub$hit_sign <- ifelse(sub$mean_hit_score >= 0, 'positive', 'negative')
  sub$contrast <- ctr
  plot_list[[ctr]] <- sub

  # Bimodality diagnostics (reported, not a gate)
  dip_res <- dip.test(sub$mean_hit_score)
  mc_res  <- Mclust(sub$mean_hit_score, G = 1:2, verbose = FALSE)
  diag_list[[ctr]] <- data.frame(
    analysis = analysis, contrast = ctr, n = nrow(sub),
    dip_stat = dip_res$statistic, dip_pval = dip_res$p.value,
    mclust_best_G = mc_res$G, bic_G1 = mc_res$BIC[1], bic_G2 = mc_res$BIC[2],
    stringsAsFactors = FALSE)

  # Test A: Wilcoxon on the continuous axis
  wt <- wilcox.test(cum_abs_estimate_total ~ hit_sign, data = sub, exact = FALSE)
  test_list[[paste0(ctr, '__wilcox')]] <- data.frame(
    analysis = analysis, contrast = ctr,
    test = 'Wilcoxon: cum_abs_estimate_total ~ hit_score sign',
    n = nrow(sub), n_pos = sum(sub$hit_sign == 'positive'),
    n_neg = sum(sub$hit_sign == 'negative'),
    OR = NA_real_, CI_lo = NA_real_, CI_hi = NA_real_, pval = wt$p.value,
    stringsAsFactors = FALSE)

  # Test B: Fisher on the near-binary axis
  ct <- table(factor(sub$has_triplet, c(TRUE, FALSE)),
              factor(sub$hit_sign,    c('positive', 'negative')))
  ft <- fisher.test(ct)
  test_list[[paste0(ctr, '__fisher')]] <- data.frame(
    analysis = analysis, contrast = ctr,
    test = 'Fisher: triplet_pairs>0 x hit_score sign',
    n = nrow(sub), n_pos = sum(sub$hit_sign == 'positive'),
    n_neg = sum(sub$hit_sign == 'negative'),
    OR = ft$estimate, CI_lo = ft$conf.int[1], CI_hi = ft$conf.int[2],
    pval = ft$p.value, stringsAsFactors = FALSE)
}

diag_df <- do.call(rbind, diag_list); rownames(diag_df) <- NULL
write.csv(diag_df, file.path(out_folder, 'bimodality_diagnostics.csv'),
          row.names = FALSE)

test_df <- do.call(rbind, test_list); rownames(test_df) <- NULL
# BH within each test family separately (3 Wilcoxon, 3 Fisher) - not pooled.
test_df$fdr <- NA_real_
for (tst in unique(test_df$test)) {
  idx <- test_df$test == tst
  test_df$fdr[idx] <- p.adjust(test_df$pval[idx], method = 'BH')
}
test_df <- test_df[order(test_df$test, test_df$fdr), ]
write.csv(test_df, file.path(out_folder, 'hitscore_sign_tests.csv'), row.names = FALSE)
message('  Saved: bimodality_diagnostics.csv, hitscore_sign_tests.csv',
        ' (FDR<0.05: ', sum(test_df$fdr < 0.05, na.rm = TRUE), ' / ',
        nrow(test_df), ')')

plot_df <- do.call(rbind, plot_list); rownames(plot_df) <- NULL
write.csv(plot_df, file.path(out_folder, 'hitscore_sign_enhancers.csv'),
          row.names = FALSE)
message('  Saved: hitscore_sign_enhancers.csv (', nrow(plot_df), ' rows)')

message('Done.')
