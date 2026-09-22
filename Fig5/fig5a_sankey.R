#### Fig 5a — GRN edge flows as a custom Sankey / Alluvial diagram ####
#
# Traces how regulatory information flows through the three-tier GRN
# structure (TF -> enhancer -> promoter), classifying enhancers by their role
# in that flow: "sink" enhancers are regulated by a TF but do not themselves
# regulate a promoter; "mediating" enhancers are both regulated by a TF and
# regulate a promoter, i.e. genuine intermediaries; "source" enhancers
# regulate a promoter without being regulated by a TF in this network. Each
# flow is further split by whether the enhancer/promoter end is expressed
# (tCRE) or accessibility-only (aCRE). This figure shows ALL enhancer-
# promoter edges, not just positive-estimate (activating) ones, so both
# activating and repressive regulatory routes are represented.
#
# Documented EXCEPTION to the Data_and_code split used elsewhere in Fig5/Fig6/
# ExtFig7: this script stays fully self-contained (reads the raw network +
# dataset directly) because the Sankey ribbon geometry needs dozens of
# fine-grained intermediate counts that don't decompose into one clean
# summary table.
#
# Reads (see ../Data_and_code/config.R to set this path):
#   [primary_data_folder]/5_lasso_models/lasso_network_with_randomization.rds
#   [primary_data_folder]/4_preparing_dataset/dataset.rds
#   ../Data_and_code/network_loading/network_thresholds.txt
#
# Output (Fig5/out/): f5a.sankey_GRN.png / .pdf
#
# Cross-check (reported by this script's own messages): 327 TFs, 20,825 sink
# enhancers, 10,813 mediating enhancers, 20,696 source enhancers, 18,327
# promoters, TF->promoter n=31,130, TF->enhancer n=14,221 (=TF->sink +
# TF->mediating), downstream mediating->promoter n=23,695.


# Clear this environment before running, EXCEPT any control-panel override a
# caller (e.g. assemble_fig5.R) may have pre-set here before sourcing this
# script -- currently just `include_bottom_text` (see below). A plain
# `rm(list = ls())` would silently wipe that override too, since it runs
# before the override is ever read.
rm(list = setdiff(ls(), 'include_bottom_text'))
library(tidyverse)
library(ggplot2)

source('../Data_and_code/network_loading/plot_theme_and_save.R')
# Note: this panel does NOT apply bimodal_theme to the Sankey itself -- the
# plot is theme_void() by design (custom geom_rect/geom_path ribbon geometry,
# no axes) -- but it uses the shared save_panel_png_pdf() helper below for
# the standard PNG+PDF-at-same-size output convention used across this repo.

# ── CONTROL PANEL ─────────────────────────────────────────────────────────────
source('../Data_and_code/config.R')  # defines primary_data_folder

network_file <- file.path(primary_data_folder, '5_lasso_models/lasso_network_with_randomization.rds')
dataset_file <- file.path(primary_data_folder, '4_preparing_dataset/dataset.rds')

out_dir        <- 'out'
thresholds_file <- '../Data_and_code/network_loading/network_thresholds.txt'

apply_estimate_threshold        <- TRUE
positive_only_enhancer_promoter <- FALSE   # keep ALL EP edges, not positive-only
show_sink_source                 <- TRUE

# The full statistics block below the diagram duplicates what the compact
# on-node labels (built further down, e.g. "Sink enhancer\nn=20,825\n
# ATAC=17,461, TSS=3,364") already show in-panel -- it exists for
# traceability/QC when running this script standalone. TRUE unless a caller
# (e.g. an assemble_fig5.R-style figure-assembly script, which lays this
# panel out at a fraction of its native 12x7in save size, where the block's
# fixed-point-size text would visually overlap) has already set this
# variable before sourcing this script.
if (!exists('include_bottom_text')) include_bottom_text <- TRUE

node_w      <- 0.15
box_fn      <- log2
ribbon_fn   <- log2
scale_fn0   <- function(x) ifelse(x == 0, 0, box_fn(x))
ribbon_fn0  <- function(x) ifelse(x == 0, 0, ribbon_fn(x))

palette     <- 'B'
alpha_main  <- 0.85
alpha_sink  <- 0.35
alpha_src   <- 0.35

text_font      <- 'Arial'
text_size      <- 2.8
n_ribbon_pts   <- 1000

x_tf   <- 1
x_sink <- 2
x_med  <- 3
x_src  <- 4
x_prom <- 5

y_tf   <- 35
y_sink <- 82
y_med  <- 55
y_src  <- NA
y_prom <- 15
# ──────────────────────────────────────────────────────────────────────────────

palettes <- list(
  B = list(
    col_tf            = '#707070',
    col_acre          = '#5C3D99',
    col_tcre          = '#D4920A',
    col_node_tf       = '#5B7FA6',
    col_node_prom     = '#B5714A',
    col_node_enh_acre = '#C4B0E8',
    col_node_enh_tcre = '#F0D890'
  ),
  C = list(
    col_tf            = '#555555',
    col_acre          = '#2255CC',
    col_tcre          = '#CC2288',
    col_node_tf       = '#4A6A9A',
    col_node_prom     = '#3A7A55',
    col_node_enh_acre = '#B0C8F0',
    col_node_enh_tcre = '#F0B0D8'
  )
)
col_tf            <- palettes[[palette]]$col_tf
col_acre          <- palettes[[palette]]$col_acre
col_tcre          <- palettes[[palette]]$col_tcre
col_node_tf       <- palettes[[palette]]$col_node_tf
col_node_prom     <- palettes[[palette]]$col_node_prom
col_node_enh_acre <- palettes[[palette]]$col_node_enh_acre
col_node_enh_tcre <- palettes[[palette]]$col_node_enh_tcre

label_tf                 <- 'TF'
label_sink_enhancer      <- 'Sink enhancer'
label_mediating_enhancer <- 'Mediating enhancer'
label_source_enhancer    <- 'Source enhancer'
label_promoter           <- 'Promoter'

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#### Load thresholds ####

thr_lines <- readLines(thresholds_file)
thr_lines <- thr_lines[!grepl('^\\s*#', thr_lines) & grepl('=', thr_lines)]
thr_parsed <- setNames(
  as.numeric(trimws(sub('.*=', '', thr_lines))),
  trimws(sub('=.*', '', thr_lines))
)

adjpvalue_threshold <- unname(thr_parsed['adjpvalue_threshold'])
dev_ratio_threshold <- unname(thr_parsed['dev_ratio_threshold'])
abs_estimate_thresholds <- c(
  TF_Enhancer_aCRE = unname(thr_parsed['abs_estimate_TF_peak_Enhancer_aCRE']),
  TF_Enhancer_tCRE = unname(thr_parsed['abs_estimate_TF_peak_Enhancer_tCRE']),
  TF_Promoter_aCRE = unname(thr_parsed['abs_estimate_TF_peak_Promoter_aCRE']),
  TF_Promoter_tCRE = unname(thr_parsed['abs_estimate_TF_peak_Promoter_tCRE']),
  EP_aCRE_aCRE     = unname(thr_parsed['abs_estimate_Enhancer_promoter_aCRE_aCRE']),
  EP_aCRE_tCRE     = unname(thr_parsed['abs_estimate_Enhancer_promoter_aCRE_tCRE']),
  EP_tCRE_aCRE     = unname(thr_parsed['abs_estimate_Enhancer_promoter_tCRE_aCRE']),
  EP_tCRE_tCRE     = unname(thr_parsed['abs_estimate_Enhancer_promoter_tCRE_tCRE'])
)
message('Thresholds loaded from: ', thresholds_file)

#### Load and filter network ####

message('Loading network...')
network <- readRDS(network_file)

message('Loading dataset to identify expressed peaks (tCRE)...')
dataset_rownames <- rownames(readRDS(dataset_file))
expressed_peaks  <- sub('_tCRE$', '',
                        dataset_rownames[grepl('^chr', dataset_rownames) &
                                         grepl('_tCRE$', dataset_rownames)])

network <- network[network$adjpvalue <= adjpvalue_threshold, ]
network <- network[network$dev_ratio >= dev_ratio_threshold, ]

if (apply_estimate_threshold) {
  tf_key  <- paste0('TF_', network$target_type, '_', network$target_measurement_type)
  tf_thr  <- abs_estimate_thresholds[tf_key]
  ep_key  <- paste0('EP_', network$source_measurement_type, '_',
                    network$target_measurement_type)
  ep_thr  <- abs_estimate_thresholds[ep_key]
  thr <- ifelse(network$link_type == 'TF_peak', tf_thr, ep_thr)
  network <- network[!is.na(thr) & abs(network$estimate) >= thr, ]
}

if (positive_only_enhancer_promoter) {
  to_remove <- which(network$link_type == 'Enhancer_promoter' & network$estimate < 0)
  if (length(to_remove) > 0) network <- network[-to_remove, ]
}

message('  Edges after filtering: ', nrow(network))

#### Classify enhancers ####

regulated_enhancers  <- unique(network$target[network$link_type == 'TF_peak' &
                                               network$target_type == 'Enhancer'])
regulating_enhancers <- unique(network$source[network$link_type == 'Enhancer_promoter' &
                                               network$source_type == 'Enhancer'])

sink_enhancers      <- setdiff(regulated_enhancers,  regulating_enhancers)
mediating_enhancers <- intersect(regulated_enhancers, regulating_enhancers)
source_enhancers    <- setdiff(regulating_enhancers, regulated_enhancers)

message('  Sink enhancers: ',      length(sink_enhancers),
        '  Mediating: ',           length(mediating_enhancers),
        '  Source: ',              length(source_enhancers))

n_tcre_sink <- sum(sink_enhancers      %in% expressed_peaks)
n_tcre_med  <- sum(mediating_enhancers %in% expressed_peaks)
n_tcre_src  <- sum(source_enhancers    %in% expressed_peaks)
n_acre_sink <- length(sink_enhancers)      - n_tcre_sink
n_acre_med  <- length(mediating_enhancers) - n_tcre_med
n_acre_src  <- length(source_enhancers)    - n_tcre_src

vis_frac_tcre <- function(n_tcre, n_acre) {
  s <- ribbon_fn0(n_tcre) + ribbon_fn0(n_acre)
  if (s == 0) return(0)
  ribbon_fn0(n_tcre) / s
}
vis_tcre_sink <- vis_frac_tcre(n_tcre_sink, n_acre_sink)
vis_tcre_med  <- vis_frac_tcre(n_tcre_med,  n_acre_med)
vis_tcre_src  <- vis_frac_tcre(n_tcre_src,  n_acre_src)

#### Count edges per flow, split by enhancer expression status ####

expressed_sink <- intersect(sink_enhancers,      expressed_peaks)
expressed_med  <- intersect(mediating_enhancers,  expressed_peaks)
expressed_src  <- intersect(source_enhancers,     expressed_peaks)
nonexpr_sink   <- setdiff(sink_enhancers,      expressed_peaks)
nonexpr_med    <- setdiff(mediating_enhancers,  expressed_peaks)
nonexpr_src    <- setdiff(source_enhancers,     expressed_peaks)

col0 <- function(df, col) if (col %in% names(df)) df[[col]] else 0L

count_tf_to_set <- function(df, target_set) {
  df %>%
    filter(link_type == 'TF_peak', target %in% target_set) %>%
    count(target_measurement_type) %>%
    pivot_wider(names_from = target_measurement_type, values_from = n, values_fill = 0)
}

count_ep_from_set <- function(df, source_set) {
  df %>%
    filter(link_type == 'Enhancer_promoter', source %in% source_set) %>%
    count(source_measurement_type, target_measurement_type) %>%
    mutate(key = paste0(source_measurement_type, '_', target_measurement_type)) %>%
    select(key, n) %>%
    pivot_wider(names_from = key, values_from = n, values_fill = 0)
}

tf_sink_expr_c   <- count_tf_to_set(network, expressed_sink)
tf_sink_nexpr_c  <- count_tf_to_set(network, nonexpr_sink)
n_tf_sink_expr_acre  <- col0(tf_sink_expr_c,  'aCRE')
n_tf_sink_expr_tcre  <- col0(tf_sink_expr_c,  'tCRE')
n_tf_sink_expr       <- n_tf_sink_expr_acre + n_tf_sink_expr_tcre
n_tf_sink_nexpr_acre <- col0(tf_sink_nexpr_c, 'aCRE')
n_tf_sink_nexpr_tcre <- col0(tf_sink_nexpr_c, 'tCRE')
n_tf_sink_nexpr      <- n_tf_sink_nexpr_acre + n_tf_sink_nexpr_tcre
n_tf_sink            <- n_tf_sink_expr + n_tf_sink_nexpr

tf_med_expr_c   <- count_tf_to_set(network, expressed_med)
tf_med_nexpr_c  <- count_tf_to_set(network, nonexpr_med)
n_tf_med_expr_acre  <- col0(tf_med_expr_c,  'aCRE')
n_tf_med_expr_tcre  <- col0(tf_med_expr_c,  'tCRE')
n_tf_med_expr       <- n_tf_med_expr_acre + n_tf_med_expr_tcre
n_tf_med_nexpr_acre <- col0(tf_med_nexpr_c, 'aCRE')
n_tf_med_nexpr_tcre <- col0(tf_med_nexpr_c, 'tCRE')
n_tf_med_nexpr      <- n_tf_med_nexpr_acre + n_tf_med_nexpr_tcre
n_tf_mediating      <- n_tf_med_expr + n_tf_med_nexpr

tf_prom_c <- network %>%
  filter(link_type == 'TF_peak', target_type == 'Promoter') %>%
  count(target_measurement_type) %>%
  pivot_wider(names_from = target_measurement_type, values_from = n, values_fill = 0)
n_tf_prom_acre <- col0(tf_prom_c, 'aCRE')
n_tf_prom_tcre <- col0(tf_prom_c, 'tCRE')
n_tf_promoter  <- n_tf_prom_acre + n_tf_prom_tcre

med_prom_expr_c  <- count_ep_from_set(network, expressed_med)
med_prom_nexpr_c <- count_ep_from_set(network, nonexpr_med)
n_med_expr_prom_aa  <- col0(med_prom_expr_c,  'aCRE_aCRE')
n_med_expr_prom_at  <- col0(med_prom_expr_c,  'aCRE_tCRE')
n_med_expr_prom_ta  <- col0(med_prom_expr_c,  'tCRE_aCRE')
n_med_expr_prom_tt  <- col0(med_prom_expr_c,  'tCRE_tCRE')
n_med_expr_prom     <- n_med_expr_prom_aa + n_med_expr_prom_at +
                       n_med_expr_prom_ta + n_med_expr_prom_tt
n_med_nexpr_prom_aa <- col0(med_prom_nexpr_c, 'aCRE_aCRE')
n_med_nexpr_prom_at <- col0(med_prom_nexpr_c, 'aCRE_tCRE')
n_med_nexpr_prom_ta <- col0(med_prom_nexpr_c, 'tCRE_aCRE')
n_med_nexpr_prom_tt <- col0(med_prom_nexpr_c, 'tCRE_tCRE')
n_med_nexpr_prom    <- n_med_nexpr_prom_aa + n_med_nexpr_prom_at +
                       n_med_nexpr_prom_ta + n_med_nexpr_prom_tt
n_mediating_promoter <- n_med_expr_prom + n_med_nexpr_prom

src_prom_expr_c  <- count_ep_from_set(network, expressed_src)
src_prom_nexpr_c <- count_ep_from_set(network, nonexpr_src)
n_src_expr_prom_aa  <- col0(src_prom_expr_c,  'aCRE_aCRE')
n_src_expr_prom_at  <- col0(src_prom_expr_c,  'aCRE_tCRE')
n_src_expr_prom_ta  <- col0(src_prom_expr_c,  'tCRE_aCRE')
n_src_expr_prom_tt  <- col0(src_prom_expr_c,  'tCRE_tCRE')
n_src_expr_prom     <- n_src_expr_prom_aa + n_src_expr_prom_at +
                       n_src_expr_prom_ta + n_src_expr_prom_tt
n_src_nexpr_prom_aa <- col0(src_prom_nexpr_c, 'aCRE_aCRE')
n_src_nexpr_prom_at <- col0(src_prom_nexpr_c, 'aCRE_tCRE')
n_src_nexpr_prom_ta <- col0(src_prom_nexpr_c, 'tCRE_aCRE')
n_src_nexpr_prom_tt <- col0(src_prom_nexpr_c, 'tCRE_tCRE')
n_src_nexpr_prom    <- n_src_nexpr_prom_aa + n_src_nexpr_prom_at +
                       n_src_nexpr_prom_ta + n_src_nexpr_prom_tt
n_source_promoter   <- n_src_expr_prom + n_src_nexpr_prom

#### Build flow summary table ####

flows <- data.frame(
  from = c(label_tf, label_tf, label_tf,
           label_mediating_enhancer, label_source_enhancer),
  to   = c(label_sink_enhancer, label_mediating_enhancer, label_promoter,
           label_promoter, label_promoter),
  n    = c(n_tf_sink, n_tf_mediating, n_tf_promoter,
           n_mediating_promoter, n_source_promoter),
  stringsAsFactors = FALSE
)
flows <- flows[flows$n > 0, ]

n_unique_tf        <- length(unique(network$source[network$link_type == 'TF_peak' &
                                                    network$source_type == 'TF']))
n_unique_sink      <- length(sink_enhancers)
n_unique_mediating <- length(mediating_enhancers)
n_unique_promoter  <- length(unique(network$target[network$target_type == 'Promoter']))
n_unique_source    <- length(source_enhancers)

message('  Unique TFs: ', n_unique_tf,
        '  Promoters: ', n_unique_promoter)
message('  TF->promoter n=', n_tf_promoter,
        '  TF->mediating-enhancer n=', n_tf_mediating,
        ' (TF->sink n=', n_tf_sink, ', TF->all-enhancer n=', n_tf_sink + n_tf_mediating, ')',
        '  mediating->promoter n=', n_mediating_promoter)

#### Build Sankey geometry ####

plot_height <- 100
bar_max     <- plot_height * 0.25

h_tf <- h_sink <- h_mediating <- h_source <- h_promoter <- bar_max

tf_ymin   <- y_tf   * plot_height / 100 - h_tf        / 2
tf_ymax   <- y_tf   * plot_height / 100 + h_tf        / 2
sink_ymin <- y_sink * plot_height / 100 - h_sink      / 2
sink_ymax <- y_sink * plot_height / 100 + h_sink      / 2
med_ymin  <- y_med  * plot_height / 100 - h_mediating / 2
med_ymax  <- y_med  * plot_height / 100 + h_mediating / 2
prom_ymin <- y_prom * plot_height / 100 - h_promoter  / 2
prom_ymax <- y_prom * plot_height / 100 + h_promoter  / 2

if (is.na(y_src)) {
  src_ymax <- sink_ymax
  src_ymin <- src_ymax - h_source
} else {
  src_ymin <- y_src * plot_height / 100 - h_source / 2
  src_ymax <- y_src * plot_height / 100 + h_source / 2
}

sink_y_split <- sink_ymin + (sink_ymax - sink_ymin) * (1 - vis_tcre_sink)
med_y_split  <- med_ymin  + (med_ymax  - med_ymin)  * (1 - vis_tcre_med)
src_y_split  <- src_ymin  + (src_ymax  - src_ymin)  * (1 - vis_tcre_src)

sink_h_acre <- sink_y_split - sink_ymin;  sink_h_tcre <- sink_ymax - sink_y_split
med_h_acre  <- med_y_split  - med_ymin;   med_h_tcre  <- med_ymax  - med_y_split
src_h_acre  <- src_y_split  - src_ymin;   src_h_tcre  <- src_ymax  - src_y_split

scaled_bands <- function(counts, y_bot, y_top) {
  sc    <- ribbon_fn0(counts)
  fracs <- sc / sum(sc)
  bots  <- y_bot + cumsum(c(0, fracs[-length(fracs)])) * (y_top - y_bot)
  tops  <- y_bot + cumsum(fracs) * (y_top - y_bot)
  tops[length(tops)] <- y_top
  list(bots = bots, tops = tops)
}

src_prom_vis <- if (show_sink_source) n_source_promoter else 0
pb <- scaled_bands(c(n_tf_promoter, n_mediating_promoter, src_prom_vis), prom_ymin, prom_ymax)
prom_tf_bot  <- pb$bots[1];  prom_tf_top  <- pb$tops[1]
prom_med_bot <- pb$bots[2];  prom_med_top <- pb$tops[2]
prom_src_bot <- pb$bots[3];  prom_src_top <- pb$tops[3]

tf_sink_vis <- if (show_sink_source) n_tf_sink else 0
tb <- scaled_bands(c(n_tf_promoter, n_tf_mediating, tf_sink_vis), tf_ymin, tf_ymax)
tf_prom_bot <- tb$bots[1];  tf_prom_top <- tb$tops[1]
tf_med_bot  <- tb$bots[2];  tf_med_top  <- tb$tops[2]
tf_sink_bot <- tb$bots[3];  tf_sink_top <- tb$tops[3]

stack_bands <- function(counts, y_bot, y_top) {
  total <- sum(counts)
  if (total == 0) {
    bots <- rep(y_bot, length(counts))
    tops <- rep(y_bot, length(counts))
  } else {
    scale <- (y_top - y_bot) / total
    bots  <- y_bot + cumsum(c(0, counts[-length(counts)])) * scale
    tops  <- y_bot + cumsum(counts) * scale
    tops[length(tops)] <- y_top
  }
  setNames(as.list(c(bots, tops)),
           c(paste0(names(counts), '_bot'), paste0(names(counts), '_top')))
}

if (show_sink_source) {
  tf_sink_bands <- stack_bands(
    c(nexpr_acre = n_tf_sink_nexpr_acre, nexpr_tcre = n_tf_sink_nexpr_tcre,
      expr_acre  = n_tf_sink_expr_acre,  expr_tcre  = n_tf_sink_expr_tcre),
    tf_sink_bot, tf_sink_top)
}

tf_med_bands <- stack_bands(
  c(nexpr_acre = n_tf_med_nexpr_acre, nexpr_tcre = n_tf_med_nexpr_tcre,
    expr_acre  = n_tf_med_expr_acre,  expr_tcre  = n_tf_med_expr_tcre),
  tf_med_bot, tf_med_top)

if (show_sink_source) {
  sink_nexpr_in_bands <- stack_bands(
    c(acre = n_tf_sink_nexpr_acre, tcre = n_tf_sink_nexpr_tcre),
    sink_ymin, sink_y_split)
  sink_expr_in_bands <- stack_bands(
    c(acre = n_tf_sink_expr_acre, tcre = n_tf_sink_expr_tcre),
    sink_y_split, sink_ymax)
}

med_nexpr_in_bands <- stack_bands(
  c(acre = n_tf_med_nexpr_acre, tcre = n_tf_med_nexpr_tcre),
  med_ymin, med_y_split)
med_expr_in_bands <- stack_bands(
  c(acre = n_tf_med_expr_acre, tcre = n_tf_med_expr_tcre),
  med_y_split, med_ymax)

prom_med_bands <- stack_bands(
  c(nexpr_aa = n_med_nexpr_prom_aa, nexpr_at = n_med_nexpr_prom_at,
    nexpr_ta = n_med_nexpr_prom_ta, nexpr_tt = n_med_nexpr_prom_tt,
    expr_aa  = n_med_expr_prom_aa,  expr_at  = n_med_expr_prom_at,
    expr_ta  = n_med_expr_prom_ta,  expr_tt  = n_med_expr_prom_tt),
  prom_med_bot, prom_med_top)

prom_src_bands <- stack_bands(
  c(nexpr_aa = n_src_nexpr_prom_aa, nexpr_at = n_src_nexpr_prom_at,
    nexpr_ta = n_src_nexpr_prom_ta, nexpr_tt = n_src_nexpr_prom_tt,
    expr_aa  = n_src_expr_prom_aa,  expr_at  = n_src_expr_prom_at,
    expr_ta  = n_src_expr_prom_ta,  expr_tt  = n_src_expr_prom_tt),
  prom_src_bot, prom_src_top)

med_nexpr_out_bands <- stack_bands(
  c(aa = n_med_nexpr_prom_aa, at = n_med_nexpr_prom_at,
    ta = n_med_nexpr_prom_ta, tt = n_med_nexpr_prom_tt),
  med_ymin, med_y_split)
med_expr_out_bands <- stack_bands(
  c(aa = n_med_expr_prom_aa, at = n_med_expr_prom_at,
    ta = n_med_expr_prom_ta, tt = n_med_expr_prom_tt),
  med_y_split, med_ymax)

if (show_sink_source) {
  src_nexpr_out_bands <- stack_bands(
    c(aa = n_src_nexpr_prom_aa, at = n_src_nexpr_prom_at,
      ta = n_src_nexpr_prom_ta, tt = n_src_nexpr_prom_tt),
    src_ymin, src_y_split)
  src_expr_out_bands <- stack_bands(
    c(aa = n_src_expr_prom_aa, at = n_src_expr_prom_at,
      ta = n_src_expr_prom_ta, tt = n_src_expr_prom_tt),
    src_y_split, src_ymax)
}

#### Build node rectangles ####

node_rects <- data.frame(
  label    = c(label_tf, label_sink_enhancer, label_mediating_enhancer,
               label_source_enhancer, label_promoter),
  xmin     = c(x_tf, x_sink, x_med, x_src, x_prom) - node_w / 2,
  xmax     = c(x_tf, x_sink, x_med, x_src, x_prom) + node_w / 2,
  ymin     = c(tf_ymin, sink_ymin, med_ymin, src_ymin,  prom_ymin),
  ymax     = c(tf_ymax, sink_ymax, med_ymax, src_ymax,  prom_ymax),
  n_unique = c(n_unique_tf, n_unique_sink, n_unique_mediating,
               n_unique_source, n_unique_promoter)
)

#### Compact on-node labels ####
# Per co-author request: short "label\nn=..." labels sit directly above each
# node box (matching the manuscript's own Fig5a layout), IN ADDITION TO the
# full statistics block at the bottom (built further down) -- not instead of
# it. The bottom block stays for traceability/QC while this script is used
# standalone; it is expected to be cropped out at figure-assembly time, once
# the compact on-node labels carry the same information in-panel.
node_label_text <- c(
  paste0(label_tf, '\nn=', formatC(n_unique_tf, format = 'd', big.mark = ',')),
  paste0(label_sink_enhancer, '\nn=', formatC(n_unique_sink, format = 'd', big.mark = ','),
        '\nATAC=', formatC(n_acre_sink, format = 'd', big.mark = ','),
        ', TSS=',  formatC(n_tcre_sink, format = 'd', big.mark = ',')),
  paste0(label_mediating_enhancer, '\nn=', formatC(n_unique_mediating, format = 'd', big.mark = ','),
        '\nATAC=', formatC(n_acre_med, format = 'd', big.mark = ','),
        ', TSS=',  formatC(n_tcre_med, format = 'd', big.mark = ',')),
  paste0(label_source_enhancer, '\nn=', formatC(n_unique_source, format = 'd', big.mark = ','),
        '\nATAC=', formatC(n_acre_src, format = 'd', big.mark = ','),
        ', TSS=',  formatC(n_tcre_src, format = 'd', big.mark = ',')),
  paste0(label_promoter, '\nn=', formatC(n_unique_promoter, format = 'd', big.mark = ','))
)
node_label_df <- data.frame(
  x     = (node_rects$xmin + node_rects$xmax) / 2,
  y     = node_rects$ymax + 2,
  label = node_label_text,
  stringsAsFactors = FALSE
)

enh_split <- function(xmin, xmax, ymin, ymax, frac_tcre) {
  y_split <- ymin + (ymax - ymin) * (1 - frac_tcre)
  data.frame(xmin = xmin, xmax = xmax,
             ymin = c(ymin, y_split), ymax = c(y_split, ymax),
             fill = c(col_node_enh_acre, col_node_enh_tcre),
             stringsAsFactors = FALSE)
}

add_alpha <- function(df, a) { df$alpha <- a; df }

node_fill_rects <- rbind(
  add_alpha(data.frame(xmin = node_rects$xmin[1], xmax = node_rects$xmax[1],
                       ymin = node_rects$ymin[1], ymax = node_rects$ymax[1],
                       fill = col_node_tf, stringsAsFactors = FALSE), alpha_main),
  if (show_sink_source)
    add_alpha(enh_split(node_rects$xmin[2], node_rects$xmax[2],
                        node_rects$ymin[2], node_rects$ymax[2], vis_tcre_sink), alpha_sink),
  add_alpha(enh_split(node_rects$xmin[3], node_rects$xmax[3],
                      node_rects$ymin[3], node_rects$ymax[3], vis_tcre_med),  alpha_main),
  if (show_sink_source)
    add_alpha(enh_split(node_rects$xmin[4], node_rects$xmax[4],
                        node_rects$ymin[4], node_rects$ymax[4], vis_tcre_src),  alpha_src),
  add_alpha(data.frame(xmin = node_rects$xmin[5], xmax = node_rects$xmax[5],
                       ymin = node_rects$ymin[5], ymax = node_rects$ymax[5],
                       fill = col_node_prom, stringsAsFactors = FALSE), alpha_main)
)

#### Gradient ribbon helpers ####

hex_interp <- function(c1, c2, frac) {
  r1 <- as.numeric(col2rgb(c1))
  r2 <- as.numeric(col2rgb(c2))
  r  <- round(outer(r1, 1 - frac) + outer(r2, frac))
  rgb(r[1,], r[2,], r[3,], maxColorValue = 255)
}

gradient_ribbon <- function(x0, y0_lo, y0_hi, x1, y1_lo, y1_hi,
                            col_left, col_right, sub_id, n_pts = n_ribbon_pts) {
  t       <- seq(0, 1, length.out = n_pts + 1)
  sig     <- 1 / (1 + exp(-12 * (t - 0.5)))
  xs      <- x0 + (x1 - x0) * t
  top_y   <- y0_hi + (y1_hi - y0_hi) * sig
  bot_y   <- y0_lo + (y1_lo - y0_lo) * sig
  sig_mid <- (sig[-length(sig)] + sig[-1]) / 2
  data.frame(xmin = xs[-length(xs)], xmax = xs[-1],
             ymin = bot_y[-length(bot_y)], ymax = top_y[-length(top_y)],
             fill = hex_interp(col_left, col_right, sig_mid),
             sub  = sub_id, stringsAsFactors = FALSE)
}

sigmoid_path <- function(x0, y0, x1, y1, n_pts = 320) {
  t   <- seq(0, 1, length.out = n_pts)
  sig <- 1 / (1 + exp(-12 * (t - 0.5)))
  data.frame(x = x0 + (x1 - x0) * t, y = y0 + (y1 - y0) * sig)
}

stack_ribbons <- function(x0, y0_lo, y0_hi, x1, y1_lo, y1_hi,
                          counts, cols_left, cols_right, prefix,
                          label_threshold = 0.1) {
  total_left  <- y0_hi - y0_lo
  total_right <- y1_hi - y1_lo
  total_n     <- sum(counts)
  if (total_n == 0) return(list(slices = data.frame(), separators = data.frame(),
                                labels = data.frame()))
  sqrt_counts <- ribbon_fn0(counts)
  sqrt_total  <- sum(sqrt_counts)
  vis_fracs   <- sqrt_counts / sqrt_total
  true_fracs  <- counts / total_n

  cur_left  <- y0_lo
  cur_right <- y1_lo
  slices     <- vector('list', length(counts))
  separators <- vector('list', length(counts) - 1)
  labels     <- vector('list', length(counts))

  for (k in seq_along(counts)) {
    h_left  <- total_left  * vis_fracs[k]
    h_right <- total_right * vis_fracs[k]
    slices[[k]] <- gradient_ribbon(
      x0, cur_left, cur_left + h_left,
      x1, cur_right, cur_right + h_right,
      cols_left[k], cols_right[k],
      sub_id = paste0(prefix, '_', names(counts)[k]))
    if (k < length(counts)) {
      sep <- sigmoid_path(x0, cur_left + h_left, x1, cur_right + h_right)
      sep$group <- paste0(prefix, '_sep', k)
      separators[[k]] <- sep
    }
    if (true_fracs[k] < label_threshold && counts[k] > 0) {
      x_mid <- (x0 + x1) / 2
      y_bot <- (cur_left          + cur_right)           / 2
      y_top <- (cur_left + h_left + cur_right + h_right) / 2
      labels[[k]] <- data.frame(x = x_mid, y = (y_bot + y_top) / 2,
                                label = paste0(round(true_fracs[k] * 100, 1), '%'),
                                stringsAsFactors = FALSE)
    }
    cur_left  <- cur_left  + h_left
    cur_right <- cur_right + h_right
  }
  list(slices     = do.call(rbind, slices),
       separators = do.call(rbind, separators),
       labels     = do.call(rbind, labels))
}

sep_path <- function(x0, y0, x1, y1, grp) {
  d <- sigmoid_path(x0, y0, x1, y1); d$group <- grp; d
}

enh_prom_col_left  <- c(aa = col_acre, at = col_acre, ta = col_tcre, tt = col_tcre)
enh_prom_col_right <- c(aa = col_acre, at = col_tcre, ta = col_acre, tt = col_tcre)

make_tf_enh_ribbons <- function(tf_bands, nexpr_in, expr_in,
                                x_tf_r, x_enh_l, enh_y_split, prefix) {
  col_right <- c(nexpr_acre = col_acre, nexpr_tcre = col_tcre,
                 expr_acre  = col_acre, expr_tcre  = col_tcre)
  enh_in <- list(
    nexpr_acre = list(bot = nexpr_in$acre_bot, top = nexpr_in$acre_top),
    nexpr_tcre = list(bot = nexpr_in$tcre_bot, top = nexpr_in$tcre_top),
    expr_acre  = list(bot = expr_in$acre_bot,  top = expr_in$acre_top),
    expr_tcre  = list(bot = expr_in$tcre_bot,  top = expr_in$tcre_top)
  )
  keys <- names(enh_in)
  ribbons <- lapply(keys, function(k)
    gradient_ribbon(x_tf_r,  tf_bands[[paste0(k, '_bot')]], tf_bands[[paste0(k, '_top')]],
                    x_enh_l, enh_in[[k]]$bot, enh_in[[k]]$top,
                    col_left = col_tf, col_right = col_right[k],
                    sub_id = paste0(prefix, '_', k)))
  seps <- lapply(seq_len(length(keys) - 1), function(i) {
    k     <- keys[i]; knext <- keys[i + 1]
    sep_path(x_tf_r,  tf_bands[[paste0(k, '_top')]],
             x_enh_l, enh_in[[knext]]$bot, paste0(prefix, '_sep', i))
  })
  mid_sep <- sep_path(x_tf_r,  tf_bands$expr_acre_bot,
                      x_enh_l, enh_y_split, paste0(prefix, '_mid_sep'))
  list(slices = do.call(rbind, ribbons),
       seps   = rbind(do.call(rbind, seps), mid_sep))
}

#### Build all ribbon slices ####

message('Building ribbon geometry...')

if (show_sink_source) {
  tf_sink_r <- make_tf_enh_ribbons(
    tf_sink_bands, sink_nexpr_in_bands, sink_expr_in_bands,
    x_tf + node_w/2, x_sink - node_w/2, sink_y_split, 'tf_sink')
}

tf_med_r <- make_tf_enh_ribbons(
  tf_med_bands, med_nexpr_in_bands, med_expr_in_bands,
  x_tf + node_w/2, x_med - node_w/2, med_y_split, 'tf_med')

tf_prom_r <- stack_ribbons(
  x_tf + node_w/2, tf_prom_bot, tf_prom_top,
  x_prom - node_w/2, prom_tf_bot, prom_tf_top,
  counts    = c(aCRE = n_tf_prom_acre, tCRE = n_tf_prom_tcre),
  cols_left = c(aCRE = col_tf,   tCRE = col_tf),
  cols_right= c(aCRE = col_acre, tCRE = col_tcre),
  prefix = 'tf_prom')

make_enh_prom_ribbons <- function(x_enh, out_bands, x_pr, pr_bands, prefix, seg) {
  keys <- c('aa', 'at', 'ta', 'tt')
  do.call(rbind, lapply(keys, function(k)
    gradient_ribbon(x_enh, out_bands[[paste0(k,'_bot')]], out_bands[[paste0(k,'_top')]],
                    x_pr,  pr_bands[[paste0(seg,'_',k,'_bot')]], pr_bands[[paste0(seg,'_',k,'_top')]],
                    col_left  = enh_prom_col_left[k],
                    col_right = enh_prom_col_right[k],
                    sub_id = paste0(prefix, '_', k))))
}

med_prom_nexpr_slices <- make_enh_prom_ribbons(
  x_med + node_w/2, med_nexpr_out_bands, x_prom - node_w/2, prom_med_bands, 'med_prom_nexpr', 'nexpr')
med_prom_expr_slices  <- make_enh_prom_ribbons(
  x_med + node_w/2, med_expr_out_bands,  x_prom - node_w/2, prom_med_bands, 'med_prom_expr',  'expr')

if (show_sink_source) {
  src_prom_nexpr_slices <- make_enh_prom_ribbons(
    x_src + node_w/2, src_nexpr_out_bands, x_prom - node_w/2, prom_src_bands, 'src_prom_nexpr', 'nexpr')
  src_prom_expr_slices  <- make_enh_prom_ribbons(
    x_src + node_w/2, src_expr_out_bands,  x_prom - node_w/2, prom_src_bands, 'src_prom_expr',  'expr')
}

make_enh_prom_seps <- function(x_enh, out_bands, x_pr, pr_bands, prefix, seg) {
  keys <- c('aa', 'at', 'ta', 'tt')
  do.call(rbind, lapply(seq_len(length(keys) - 1), function(i)
    sep_path(x_enh, out_bands[[paste0(keys[i], '_top')]],
             x_pr,  pr_bands[[paste0(seg, '_', keys[i+1], '_bot')]],
             paste0(prefix, '_sep', i))))
}

med_prom_nexpr_seps <- make_enh_prom_seps(
  x_med + node_w/2, med_nexpr_out_bands, x_prom - node_w/2, prom_med_bands, 'med_prom_nexpr', 'nexpr')
med_prom_expr_seps  <- make_enh_prom_seps(
  x_med + node_w/2, med_expr_out_bands,  x_prom - node_w/2, prom_med_bands, 'med_prom_expr',  'expr')
if (show_sink_source) {
  src_prom_nexpr_seps <- make_enh_prom_seps(
    x_src + node_w/2, src_nexpr_out_bands, x_prom - node_w/2, prom_src_bands, 'src_prom_nexpr', 'nexpr')
  src_prom_expr_seps  <- make_enh_prom_seps(
    x_src + node_w/2, src_expr_out_bands,  x_prom - node_w/2, prom_src_bands, 'src_prom_expr',  'expr')
}

med_expr_sep <- sep_path(x_med + node_w/2, med_y_split,
                         x_prom - node_w/2, prom_med_bands$expr_aa_bot, 'med_expr_sep')
if (show_sink_source) {
  src_expr_sep <- sep_path(x_src + node_w/2, src_y_split,
                           x_prom - node_w/2, prom_src_bands$expr_aa_bot, 'src_expr_sep')
}

tag_alpha <- function(df, a) { df$alpha <- a; df }

all_slices <- rbind(
  if (show_sink_source) tag_alpha(tf_sink_r$slices,      alpha_sink),
  tag_alpha(tf_med_r$slices,         alpha_main),
  tag_alpha(tf_prom_r$slices,        alpha_main),
  tag_alpha(med_prom_nexpr_slices,   alpha_main),
  tag_alpha(med_prom_expr_slices,    alpha_main),
  if (show_sink_source) tag_alpha(src_prom_nexpr_slices, alpha_src),
  if (show_sink_source) tag_alpha(src_prom_expr_slices,  alpha_src)
)

all_separators <- rbind(
  if (show_sink_source) tf_sink_r$seps,
  tf_med_r$seps,
  tf_prom_r$separators,
  med_prom_nexpr_seps, med_prom_expr_seps, med_expr_sep,
  if (show_sink_source) src_prom_nexpr_seps,
  if (show_sink_source) src_prom_expr_seps,
  if (show_sink_source) src_expr_sep
)

#### Build bottom annotation text ####

pct_string <- function(counts, name_map, sep = ' | ') {
  true_pcts <- round(counts / sum(counts) * 100, 1)
  nms <- name_map[names(counts)]
  paste(paste0(nms, ': ', true_pcts, '%'), collapse = sep)
}
tf_names <- c(aCRE = 'ATAC', tCRE = 'TSS')
ep_names <- c(aa = 'ATAC → ATAC', at = 'ATAC → TSS',
              ta = 'TSS → ATAC',  tt = 'TSS → TSS')

text_line_h  <- 2.5
text_x       <- x_tf - node_w / 2
text_y_start <- prom_ymin - 5

acre_tcre_string <- function(n_acre, n_tcre) {
  paste0(', ATAC=', formatC(n_acre, format='d', big.mark=','),
         ', TSS=', formatC(n_tcre, format='d', big.mark=','))
}

text_lines <- c(
  paste0('TF (n=',                formatC(n_unique_tf,        format='d', big.mark=',')),
  paste0('Sink enhancer (n=',     formatC(n_unique_sink,      format='d', big.mark=','),
         acre_tcre_string(n_acre_sink, n_tcre_sink)),
  paste0('Mediating enhancer (n=',formatC(n_unique_mediating, format='d', big.mark=','),
         acre_tcre_string(n_acre_med, n_tcre_med)),
  paste0('Source enhancer (n=',   formatC(n_unique_source,    format='d', big.mark=','),
         acre_tcre_string(n_acre_src, n_tcre_src)),
  paste0('Promoter (n=',          formatC(n_unique_promoter,  format='d', big.mark=',')),
  'TF (Expression) | ATAC (Chromatin) | TSS (Transcription)',
  paste0('TF → Sink enhancer: ',
         pct_string(c(nexpr = n_tf_sink_nexpr, expr = n_tf_sink_expr),
                    c(nexpr = 'non-expr', expr = 'expr'))),
  paste0('TF → Mediating enhancer: ',
         pct_string(c(nexpr = n_tf_med_nexpr, expr = n_tf_med_expr),
                    c(nexpr = 'non-expr', expr = 'expr'))),
  paste0('TF → Promoter: ',
         pct_string(c(aCRE = n_tf_prom_acre, tCRE = n_tf_prom_tcre), tf_names)),
  paste0('Mediating (non-expr) → Promoter: ',
         pct_string(c(aa = n_med_nexpr_prom_aa, at = n_med_nexpr_prom_at,
                      ta = n_med_nexpr_prom_ta, tt = n_med_nexpr_prom_tt), ep_names)),
  paste0('Mediating (expr) → Promoter: ',
         pct_string(c(aa = n_med_expr_prom_aa, at = n_med_expr_prom_at,
                      ta = n_med_expr_prom_ta, tt = n_med_expr_prom_tt), ep_names)),
  paste0('Source (non-expr) → Promoter: ',
         pct_string(c(aa = n_src_nexpr_prom_aa, at = n_src_nexpr_prom_at,
                      ta = n_src_nexpr_prom_ta, tt = n_src_nexpr_prom_tt), ep_names)),
  paste0('Source (expr) → Promoter: ',
         pct_string(c(aa = n_src_expr_prom_aa, at = n_src_expr_prom_at,
                      ta = n_src_expr_prom_ta, tt = n_src_expr_prom_tt), ep_names))
)

bottom_text_df <- data.frame(
  x     = text_x,
  y     = text_y_start - (seq_along(text_lines) - 1) * text_line_h,
  label = text_lines,
  stringsAsFactors = FALSE
)

#### Plot ####

message('Rendering plot...')

p <- ggplot() +
  geom_rect(data = all_slices,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, alpha = alpha),
            fill = all_slices$fill, color = NA) +
  scale_alpha_identity() +
  geom_path(data = all_separators,
            aes(x = x, y = y, group = group),
            color = 'white', linewidth = 0.5, lineend = 'round') +
  geom_rect(data = node_fill_rects,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, alpha = alpha),
            fill = node_fill_rects$fill, color = NA) +
  geom_text(data = node_label_df,
            aes(x = x, y = y, label = label),
            hjust = 0.5, vjust = 0, lineheight = 0.9,
            size = text_size, family = text_font, fontface = 'bold') +
  { if (include_bottom_text)
      geom_text(data = bottom_text_df,
               aes(x = x, y = y, label = label),
               hjust = 0, vjust = 1, size = text_size, family = text_font) } +
  scale_x_continuous(expand = c(0.12, 0.18)) +
  scale_y_continuous(expand = expansion(
    mult = c(0.3, if (include_bottom_text) 0.1 else 0.22))) +
  theme_void() +
  theme(legend.position = 'none',
        # Locking to this panel's own saved aspect ratio (12 x 7 in) only
        # matters when the bottom text block is present -- its line spacing
        # is tuned for the 12x7in save size and visually overlaps itself if
        # the panel is stretched/squished to a different aspect. With that
        # block dropped (include_bottom_text = FALSE, the composite-script
        # case), an earlier version left this unconstrained (aspect.ratio =
        # NULL) so patchwork could freely stretch the panel to fill its
        # column -- avoiding a top-label clipping bug that a naive
        # aspect-ratio lock produced at the time. That free stretch is what
        # then produced a DIFFERENT visible defect: this panel's column in
        # the (a/c/d)|b layout is much taller than wide, and stretching a
        # wide, short diagram (few Sankey levels, no vertical content beyond
        # the node labels) to fill it visibly elongates it "over the
        # y-axis". Now that the title is dropped and top expansion (above)
        # is increased instead, a fixed non-NULL ratio close to this
        # panel's natural content shape (roughly 3:5, wider-than-tall but
        # not as extreme as the 12:7 standalone save) no longer clips the
        # top labels, so it is safe to lock again -- patchwork pads with
        # whitespace on the constrained axis rather than distorting.
        aspect.ratio = if (include_bottom_text) 7 / 12 else 3 / 5)

#### Save outputs ####

message('Saving PNG + PDF...')
save_panel_png_pdf(p, out_dir, 'f5a.sankey_GRN', width_in = 12, height_in = 7)

message('Done. Results written to: ', out_dir)

cat('\n--- Enhancer classification (aCRE / tCRE) ---\n')
cat('Sink enhancers (regulated, not regulating):    ', length(sink_enhancers),
    ' (aCRE=', n_acre_sink, ', tCRE=', n_tcre_sink, ')\n', sep = '')
cat('Mediating enhancers (regulated and regulating):', length(mediating_enhancers),
    ' (aCRE=', n_acre_med, ', tCRE=', n_tcre_med, ')\n', sep = '')
cat('Source enhancers (not regulated, regulating):  ', length(source_enhancers),
    ' (aCRE=', n_acre_src, ', tCRE=', n_tcre_src, ')\n', sep = '')
cat('\n--- Edge counts per flow ---\n')
print(flows[, c('from', 'to', 'n')])
cat('\nTF->promoter n=', n_tf_promoter,
    '  TF->mediating-enhancer n=', n_tf_mediating,
    '  mediating->promoter n=', n_mediating_promoter, '\n', sep='')
