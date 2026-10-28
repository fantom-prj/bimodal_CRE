#### Regeneration script: CRISPRi screen peaks overlapping network peaks ####
#
# Feeds (via copy_enhancer_stats.R -> Fig6/data/): Fig6d, Fig6e, ExtFig6e,
# ExtFig6f, and Data_and_code/screen_correlation/.
#
# Overlaps the peaks targeted by the CRISPRi neurogenesis screen (iPSC ->
# NSC -> cortical neuron, 10 pairwise contrasts) with the network's promoter
# and enhancer peaks, by genomic coordinates (GenomicRanges::findOverlaps).
# The overlap is purely positional, so network peaks are taken from the full,
# unfiltered lasso network; edge filtering is left to the downstream analyses.
# A screen peak listed in several cell types is collapsed to the union of its
# coordinates.
#
# Each output row is one (screen peak x network peak x contrast) triple, so
# every contrast's screen result is kept. Screen peaks present in the peak
# library but with no quantified result get a single row with NA contrast,
# logFC and hit_score. hit_score is the screen's signed combination of logFC
# and Wilcoxon significance (positive = promotes differentiation).
#
# Reads (see ../config.R to set these paths):
#   [raw_data_folder]/all_promoter_coordinate.bed
#     screen peak library: chr, start, end, peak_ID, cell_type (no header)
#   [raw_data_folder]/overall.promoter.base.result.osc
#     screen results, one row per (peak_ID, contrast)
#   [primary_data_folder]/5_lasso_models/lasso_network_with_randomization.rds
#
# Writes [primary_data_folder]/10_enhancer_centrality/:
#   screen_overlap_promoters.csv   (network promoter peaks)
#   screen_overlap_enhancers.csv   (network enhancer peaks)
#   columns: screen_peak_ID, peak_type, network_peak_id, compare,
#            gene_symbol, CAT_geneClass, logFC, wilcox_p, mlog_wilcox_p,
#            hit_score

rm(list = ls())
suppressPackageStartupMessages(library(GenomicRanges))

# ── CONTROL PANEL ─────────────────────────────────────────────────────────────
source('../config.R')  # defines primary_data_folder, raw_data_folder

screen_bed     <- file.path(raw_data_folder, 'all_promoter_coordinate.bed')
screen_results <- file.path(raw_data_folder, 'overall.promoter.base.result.osc')
network_file   <- file.path(primary_data_folder, '5_lasso_models/lasso_network_with_randomization.rds')
out_folder     <- file.path(primary_data_folder, '10_enhancer_centrality')
# ──────────────────────────────────────────────────────────────────────────────

dir.create(out_folder, showWarnings = FALSE, recursive = TRUE)

#### Load screen data ####

message('Loading screen data...')
bed <- read.table(screen_bed,
                  col.names        = c('chr', 'start', 'end', 'peak_ID', 'cell_type'),
                  stringsAsFactors = FALSE)

# one row per peak_ID: union of its coordinates across cell types
bed_uniq <- do.call(rbind, lapply(split(bed, bed$peak_ID), function(x)
  data.frame(chr     = x$chr[1],
             start   = min(x$start),
             end     = max(x$end),
             peak_ID = x$peak_ID[1],
             stringsAsFactors = FALSE)))

osc <- read.delim(screen_results, stringsAsFactors = FALSE)
colnames(osc) <- c('peak_ID', 'compare', 'gene_symbol', 'no_gRNA', 'geneClass',
                   'CAT_geneClass', 'Cutoff', 'logFC', 'wilcox_p',
                   'mlog_wilcox_p', 'hit_score')

message('  Screen peaks in library: ', nrow(bed_uniq))
message('  Contrasts: ', paste(unique(osc$compare), collapse = ', '))

gr_screen <- GRanges(
  seqnames = bed_uniq$chr,
  ranges   = IRanges(start = bed_uniq$start, end = bed_uniq$end)
)
mcols(gr_screen)$peak_ID <- bed_uniq$peak_ID

#### Network peaks ####

message('Loading network peak ids...')
net <- readRDS(network_file)

promoter_ids <- unique(c(
  net$source[net$source_type == 'Promoter'],
  net$target[net$target_type == 'Promoter']
))
enhancer_ids <- unique(c(
  net$source[net$source_type == 'Enhancer'],
  net$target[net$target_type == 'Enhancer']
))
rm(net)
message('  Network promoter peaks: ', length(promoter_ids))
message('  Network enhancer peaks: ', length(enhancer_ids))

# peak ids are <chr>_<start>_<end>; the chromosome name may itself contain '_'
parse_peaks <- function(peak_ids) {
  end_   <- as.integer(sub('.*_', '', peak_ids))
  tmp    <- sub('_[^_]+$', '', peak_ids)
  start_ <- as.integer(sub('.*_', '', tmp))
  chr_   <- sub('_[^_]+_[^_]+$', '', peak_ids)
  gr <- GRanges(seqnames = chr_,
                ranges   = IRanges(start = start_, end = end_))
  mcols(gr)$peak_id <- peak_ids
  gr
}

gr_promoters <- parse_peaks(promoter_ids)
gr_enhancers <- parse_peaks(enhancer_ids)

#### Overlaps ####

message('Finding overlaps...')
hits_prom <- findOverlaps(gr_screen, gr_promoters)
hits_enh  <- findOverlaps(gr_screen, gr_enhancers)
message('  Screen peaks overlapping a promoter: ', length(unique(queryHits(hits_prom))))
message('  Screen peaks overlapping an enhancer: ', length(unique(queryHits(hits_enh))))

make_overlap_pairs <- function(hits, gr_net, peak_type) {
  data.frame(
    peak_type       = peak_type,
    screen_peak_ID  = mcols(gr_screen)$peak_ID[queryHits(hits)],
    network_peak_id = mcols(gr_net)$peak_id[subjectHits(hits)],
    stringsAsFactors = FALSE
  )
}

osc_cols <- c('peak_ID', 'compare', 'gene_symbol', 'CAT_geneClass',
              'logFC', 'wilcox_p', 'mlog_wilcox_p', 'hit_score')
osc_sub  <- osc[, osc_cols]

ov_prom <- merge(make_overlap_pairs(hits_prom, gr_promoters, 'Promoter'), osc_sub,
                 by.x = 'screen_peak_ID', by.y = 'peak_ID', all.x = TRUE)
ov_enh  <- merge(make_overlap_pairs(hits_enh,  gr_enhancers, 'Enhancer'), osc_sub,
                 by.x = 'screen_peak_ID', by.y = 'peak_ID', all.x = TRUE)

#### Write outputs ####

write.csv(ov_prom, file.path(out_folder, 'screen_overlap_promoters.csv'), row.names = FALSE)
write.csv(ov_enh,  file.path(out_folder, 'screen_overlap_enhancers.csv'), row.names = FALSE)

message('Saved into ', out_folder, ':')
message('  screen_overlap_promoters.csv (', nrow(ov_prom), ' rows, ',
        length(unique(ov_prom$screen_peak_ID)), ' screen peaks)')
message('  screen_overlap_enhancers.csv (', nrow(ov_enh), ' rows, ',
        length(unique(ov_enh$screen_peak_ID)), ' screen peaks)')
message('Done.')
