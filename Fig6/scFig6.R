primary_folder <- '[primary_folder]'
path_fig6_data <- paste0(primary_folder, 'Fig6/data/')
path_fig6      <- paste0(primary_folder, 'Fig6/out/')
path_code      <- paste0(primary_folder, 'Data_and_code/')

library(ggplot2)
library(dplyr)
library(ggforce)

#=====================
theme1=theme(  panel.background = element_rect(fill = "white", color = "white", linewidth = 0.25, linetype = "solid"),
               panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid', color = "grey90"), 
               axis.text.x = element_text(color ="black"),
               axis.text.y = element_text(color ="black"),
               plot.title = element_text(hjust = 0.5,  color ="black", margin = margin(0.1,0.1,0.1,0.1, "cm")),
               plot.subtitle = element_text(color="black", hjust=0.5),
               text = element_text(size=6),
               strip.text = element_text(size=6),
               legend.background = element_rect(fill="white", color="black", size=0.25),
               legend.title = element_text(hjust=0.5),
               legend.text = element_text(lineheight = 0.6, margin = margin(l = 1, unit = "pt")),
               legend.key.size = unit(0.2, 'cm'),
               legend.margin = margin(0.02,0.02,0.02,0.02, "cm"),
               axis.line = element_line(linewidth = 0.25, colour = "black"),
               axis.ticks = element_line(linewidth = 0.25,colour = "black"),
               strip.background =element_rect(fill="grey90"),
               strip.text.x = element_text(margin = margin(0.02,0.02,0.02,0.02, "cm")),
               strip.text.y = element_text(margin = margin(0.02,0.02,0.02,0.02, "cm")),
               legend.key = element_blank(),
               legend.position = "right",
               legend.box.margin = margin(t=-10, r=0, b=0, l=-10),
               plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

#===============================================================================
# fig 6a

category_table <- read.csv(paste0(path_fig6_data, 'sweep_category_table.csv'))

category_colours <- c(fully_uncoupled = '#4477AA',
                      intermediate    = '#CCBB44',
                      fully_coupled   = '#EE6677')
category_labels  <- c(fully_uncoupled = 'Fully uncoupled (Jaccard = 0)',
                      intermediate    = 'Intermediate (0 < Jaccard < 1)',
                      fully_coupled   = 'Fully coupled (Jaccard = 1)')

category_table1=category_table[which(category_table$theme == "TF_enhancer" & category_table$sweep=="estimate_percentile"),]
category_table1$category=factor(category_table1$category, levels=c("fully_uncoupled","intermediate","fully_coupled"))

f6a=ggplot(category_table1, aes(x = as.character(sweep_value), y=pct, fill=category)) +
  geom_bar(linewidth=0.25, alpha=1,  width=0.8,  color="black", stat="identity") + 
  scale_fill_manual(values=category_colours, labels=category_labels)+
  geom_text(aes(label=n),position = position_stack(vjust = 0.5), size=1.8,color="black")+
  labs(x="|estimate| percentile per stratum", y = "% of enhancers", title="Uncoupled TF contribution to enhancer transcription and accessibility") +
  theme1
pdf(paste0(path_fig6,"f6a.stacked_TF_enhancer_estimate_pctl2.pdf"), width = 3.4, height = 2)
print(f6a)
dev.off()

#============================
# fig 6b

bin_df <- read.csv(paste0(path_fig6_data, 'distance_bin_with_tCRE.csv'))
bin_df$dist_bin <- factor(bin_df$dist_bin, levels = unique(bin_df$dist_bin))
bin_df$source_class <- factor(bin_df$source_class, levels = c("expressed enhancers","non-expressed enhancers"))

f6b <- ggplot(bin_df %>% dplyr::filter(n_enhancers >= 10),
              aes(x = dist_bin, y = y_mean,
                  colour = source_class, group = source_class)) +
  geom_line(linewidth = 0.4) +
  geom_point(aes(size = n_enhancers)) +
  scale_colour_manual(name = 'Enhancer type', labels = c("Transcribed","Non-transcribed"), values = c('darkorange','steelblue')) +
  scale_size_continuous(name = 'N enhancers', range = c(0.5, 2.5)) +
  scale_x_discrete(name = 'Enhancer-promoter distance bin (kb)') +
  scale_y_continuous(name = 'Mean elastic-net-selected edges per enhancer') +
  labs(title="Promoter connections per enhancer across EP distance") +
  theme1 + theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))

pdf(paste0(path_fig6,"f6b.edges_per_enhancer_distance_bin2.pdf"), width = 3.4, height = 2)
print(f6b)
dev.off()

#============================
# fig 6c

triplets_tCRE <- read.csv(paste0(path_fig6_data, 'triplets_tCRE.csv'))

#### Inner ring: new-route vs shared ####

n_new    <- sum(triplets_tCRE$new_route)
n_shared <- sum(!triplets_tCRE$new_route)
n_total  <- nrow(triplets_tCRE)

message('Triplets: ', n_total, '  new route ', round(100 * n_new / n_total, 1),
        '%  shared ', round(100 * n_shared / n_total, 1), '%')

inner_df <- data.frame(
  label = c('New route\n(TSS only)', 'Shared with\nATAC network'),
  count = c(n_new, n_shared),
  frac  = c(n_new / n_total, n_shared / n_total)
)

#### Outer ring: top 10 TFs in new-route triplets ####

top_n_tfs     <- 10
tf_new_counts <- sort(table(triplets_tCRE$TF[triplets_tCRE$new_route]), decreasing = TRUE)
top_tfs       <- names(tf_new_counts)[seq_len(min(top_n_tfs, length(tf_new_counts)))]
n_top         <- as.integer(tf_new_counts[top_tfs])
n_other_new   <- n_new - sum(n_top)

# 'Remaining' matches the reference figure's own wording for the collapsed
# TF wedge (was 'Other TFs' in an earlier version of this panel).
outer_new_df <- data.frame(
  label    = c(top_tfs, 'Remaining'),
  count    = c(n_top, n_other_new),
  category = 'new_route'
)
outer_shared_df <- data.frame(
  label    = 'Shared with\nATAC network',
  count    = n_shared,
  category = 'shared'
)
outer_df      <- rbind(outer_new_df, outer_shared_df)
outer_df$frac <- outer_df$count / n_total

message('  Remaining: ', round(100 * n_other_new / n_new, 1), '% of new-route')

#### Colours ####

tf_colours    <- colorRampPalette(c('#08519c', '#6baed6'))(length(top_tfs) + 1)
outer_colours <- c(tf_colours, 'grey85')
names(outer_colours) <- c(top_tfs, 'Remaining', 'Shared with\nATAC network')
inner_colours <- c('New route\n(TSS only)'     = '#2171b5',
                   'Shared with\nATAC network' = 'grey85')
all_colours   <- c(inner_colours, outer_colours)

#### Arcs ####
# Angles start at the top and proceed clockwise -- see the IMPORTANT note
# below on ggforce::geom_arc_bar()'s actual angle convention, which this
# start/end formula is written against.

make_arc_df <- function(fracs, labels, r0, r1) {
  end_ang   <- pi/2 - cumsum(fracs) * 2 * pi
  start_ang <- c(pi/2, head(end_ang, -1))
  data.frame(
    x0 = 0, y0 = 0, r0 = r0, r = r1,
    start = start_ang, end = end_ang,
    mid   = (start_ang + end_ang) / 2,
    label = labels, frac = fracs,
    stringsAsFactors = FALSE
  )
}

# Grown from the earlier 1.5/2.4/2.6/3.5 now that 'Remaining' no longer
# needs its own far-side callout (see below) -- that freed up the panel's
# left side, which was previously reserved for 'Remaining''s label, so the
# donut itself can be larger within the same saved canvas.
r_inner_in <- 1.9; r_inner_out <- 3.0
r_outer_in <- 3.2; r_outer_out <- 4.3

arc_inner <- make_arc_df(inner_df$frac, inner_df$label, r_inner_in, r_inner_out)
arc_outer <- make_arc_df(outer_df$frac, outer_df$label, r_outer_in, r_outer_out)

#### Leader-line callouts for the top-10 TFs + 'Remaining' ####
# One label per new-route wedge (the 'Shared with ATAC network' wedge is
# not labelled here -- its share is already stated by the centre label
# below, since the two inner-ring slices are complementary).
callout_df <- arc_outer[arc_outer$label != 'Shared with\nATAC network', ]
callout_df$pct_of_newroute <- round(100 * (outer_df$count[match(callout_df$label, outer_df$label)] / n_new), 1)
callout_df$callout_label   <- paste0(callout_df$label, ' ', callout_df$pct_of_newroute, '%')


callout_df$anchor_angle <- ifelse(callout_df$label == 'Remaining',
                                  callout_df$end + 0.2, callout_df$mid)

callout_df$anchor_x <- (r_outer_out + 0.05) * sin(callout_df$anchor_angle)
callout_df$anchor_y <- (r_outer_out + 0.05) * cos(callout_df$anchor_angle)


callout_df$sort_key <- ifelse(callout_df$label == 'Remaining', -Inf, callout_df$anchor_y)
callout_df <- callout_df[order(callout_df$sort_key, decreasing = TRUE), ]
n_callout <- nrow(callout_df)
row_top   <- r_outer_out + 3.0
row_step  <- 0.62
gap_before_remaining <- 5
row_offset <- (seq_len(n_callout) - 1) * row_step +
  ifelse(callout_df$label == 'Remaining', gap_before_remaining, 0)
callout_df$label_x <- r_outer_out + 1.4
callout_df$label_y <- row_top - row_offset
callout_df$label_hjust <- 0


callout_df$text_colour <- ifelse(callout_df$label == 'Remaining', '#2171b5',
                                 unname(outer_colours[callout_df$label]))

#### Centre label: new-route share of ALL triplets ####

centre_label <- paste0('TSS only\n', round(100 * n_new / n_total, 1), '%')

f6c <- ggplot() +
  ggforce::geom_arc_bar(data = arc_inner,
                        aes(x0 = x0, y0 = y0, r0 = r0, r = r,
                            start = start, end = end, fill = label),
                        colour = 'white', linewidth = 0.5) +
  ggforce::geom_arc_bar(data = arc_outer,
                        aes(x0 = x0, y0 = y0, r0 = r0, r = r,
                            start = start, end = end, fill = label),
                        colour = 'white', linewidth = 0.3) +
  annotate('text', x = 0, y = 0, label = centre_label,
           size = 2.4, fontface = 'bold',
           colour = '#2171b5', lineheight = 0.9) +
  geom_segment(data = callout_df,
               aes(x = anchor_x, y = anchor_y, xend = label_x, yend = label_y),
               colour = 'grey50', linewidth = 0.25) +
  geom_text(data = callout_df,
            aes(x = label_x, y = label_y, label = callout_label, hjust = label_hjust),
            colour = callout_df$text_colour, size = 1.9,
            vjust = 0.5) +
  scale_fill_manual(values = all_colours) +
  coord_fixed(clip = 'off') +
  xlim(-(r_outer_out + 0.3), r_outer_out + 5.8) +
  ylim(-(r_outer_out + 0.5), r_outer_out + 3.2) +
  theme_void() +
  theme(legend.position = 'none',
        plot.margin = margin(5, 5, 5, 5))

pdf(paste0(path_fig6,"f6c.new_route_donut.pdf"), width = 3.66, height = 3.4)
print(f6c)
dev.off()

#============================
# fig 6d

analysis <- 'all_edges'
stats       <- read.csv(paste0(path_fig6_data, 'all_enhancers_stats.csv'),
                        stringsAsFactors = FALSE)
enh_overlap <- read.csv(paste0(path_fig6_data, 'screen_overlap_enhancers.csv'),
                        stringsAsFactors = FALSE)
scored_ids <- unique(enh_overlap$network_peak_id[!is.na(enh_overlap$hit_score)])

stats$screen_coverage <- factor(
  dplyr::case_when(stats$enhancer %in% scored_ids ~ 'scored',
                   TRUE                           ~ 'not_covered'),
  levels = c('scored', 'not_covered'))

n_s  <- sum(stats$screen_coverage == 'scored')
n_nc <- sum(stats$screen_coverage == 'not_covered')
message('Enhancers: ', nrow(stats), '  scored: ', n_s, '  not covered: ', n_nc)

colour_vals <- c(scored = '#D7191C', not_covered = 'grey75')
size_vals   <- c(scored = 0.8,       not_covered = 0.25)   # source: 2.5 / 0.8
alpha_vals  <- c(scored = 0.95,      not_covered = 0.35)
label_vals  <- c(scored      = paste0('Scored (', n_s, ')'),
                 not_covered = paste0('Not scored (', n_nc, ')'))

stats0=stats[which(stats$screen_coverage=="not_covered"),]
stats1=stats[which(stats$screen_coverage!="not_covered"),]
f6d <- ggplot() +
  geom_point_rast(data=stats0, aes(x=log10(triplet_pairs + 1), y=cum_abs_estimate_total, color = screen_coverage, size  = screen_coverage, alpha = screen_coverage)) +
  geom_point_rast(data=stats1, aes(x=log10(triplet_pairs + 1), y=cum_abs_estimate_total, color = screen_coverage, size  = screen_coverage, alpha = screen_coverage)) +
  scale_color_manual(values = colour_vals, labels = label_vals,
                     name = 'Screen coverage') +
  scale_size_manual(values = size_vals, labels = label_vals,
                    name = 'Screen coverage') +
  scale_alpha_manual(values = alpha_vals, labels = label_vals,
                     name = 'Screen coverage') +
  scale_x_continuous(name = 'Triplet-pair count (log10 scale)') +
  scale_y_continuous(name = 'Cumulative |estimate| (in + out edges)') +
  labs(title="Network centrality of 52,334 enhancers") +
  theme1+theme(legend.position = "bottom")

pdf(paste0(path_fig6,"f6d.screen_coverage_scatter2.pdf"), width = 2.5, height = 3)
print(f6d)
dev.off()

#============================
# fig 6e

analysis <- 'all_edges'

contrast_order  <- c('NSC_P2vsNeuron', 'NSC_P0vsNSC_P2', 'iPSvsNSC_P0')
contrast_labels <- c(iPSvsNSC_P0    = 'iPS vs NSC P0',
                     NSC_P0vsNSC_P2 = 'NSC P0 vs NSC P2',
                     NSC_P2vsNeuron = 'NSC P2 vs Neuron')

cor_res <- read.csv(paste0(path_fig6_data, 'screen_correlation.csv'),
                    stringsAsFactors = FALSE)

message('Tests: ', nrow(cor_res), '  FDR<0.05: ',
        sum(cor_res$fdr < 0.05, na.rm = TRUE), '  n per test: ', cor_res$n[1])

cor_res$contrast_f <- factor(cor_res$contrast, levels = contrast_order)
cor_res$axis_label <- ifelse(cor_res$axis == 'triplet_pairs',
                             'Triplet-pair count', 'Cumulative |estimate|')

f6e <- ggplot(cor_res, aes(y = contrast_f, x = axis_label, fill = rho)) +
  geom_tile(color = 'grey90', linewidth = 0.4) +
  geom_tile(data = dplyr::filter(cor_res, fdr < 0.05),
            color = 'black', linewidth = 0.6, fill = NA) +
  scale_x_discrete(labels = contrast_labels, position = "top") +
  scale_fill_gradient2(low = '#2166AC', mid = 'white', high = '#B2182B',
                       midpoint = 0, name = 'Spearman rho') +
  geom_text(aes(label = paste0('rho=', round(rho, 2), '\n',
                               'FDR=', signif(fdr, 2))),
            size = 1.8) +
  labs(title="Correlations between CRISPRi hit\nscores and enhancer centrality")+
  labs(x = NULL, y = NULL) +
  theme1 +
  theme(axis.text.x     = element_text(angle = 25, hjust = 0),
        legend.key.size = unit(0.4, 'cm'),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom")

pdf(paste0(path_fig6,"f6e.screen_correlation_heatmap2.pdf"), width = 2, height = 3)
print(f6e)
dev.off()

#===============================================================================


