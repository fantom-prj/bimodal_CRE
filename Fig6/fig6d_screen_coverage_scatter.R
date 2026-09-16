#### Figure 6d — enhancer centrality scatter coloured by screen coverage ####
#
# Panel:   Fig 6d. Two-axis enhancer centrality scatter (log10 triplet-pair
#          count -- number of distinct TF/promoter pairs bridged by the
#          enhancer -- vs cumulative |estimate| -- total regulatory strength
#          over all edges touching it), every network enhancer plotted,
#          coloured by whether it is covered by the CRISPRi neurogenesis
#          screen with a quantified hit score ('scored') or not. Screen peaks
#          that overlap an enhancer but carry no quantified score are merged
#          into 'not_covered', since they are not biologically distinguishable
#          from a true screen miss.
#          Analysis variant: all_edges.
# Reads:   Fig6/data/all_enhancers_stats.csv
#          Fig6/data/screen_overlap_enhancers.csv
#          (both placed there by
#           ../Data_and_code/enhancer_centrality_stats/copy_enhancer_stats.R)
# Writes:  Fig6/out/f6d.screen_coverage_scatter.{png,pdf}

library(dplyr)
library(ggplot2)

primary_folder <- '[primary_folder]'
path_fig6_data <- paste0(primary_folder, 'Fig6/data/')
path_fig6      <- paste0(primary_folder, 'Fig6/out/')
path_code      <- paste0(primary_folder, 'Data_and_code/')

source(paste0(path_code, 'network_loading/plot_theme_and_save.R'))

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

f6d <- ggplot(stats,
              aes(x     = log10(triplet_pairs + 1),
                  y     = cum_abs_estimate_total,
                  color = screen_coverage,
                  size  = screen_coverage,
                  alpha = screen_coverage)) +
  geom_point() +
  scale_color_manual(values = colour_vals, labels = label_vals,
                     name = 'Screen coverage') +
  scale_size_manual(values = size_vals, labels = label_vals,
                    name = 'Screen coverage') +
  scale_alpha_manual(values = alpha_vals, labels = label_vals,
                     name = 'Screen coverage') +
  scale_x_continuous(name = 'Triplet-pair count (log10 scale)') +
  scale_y_continuous(name = 'Cumulative |estimate| (in + out edges)') +
  ggtitle(paste0('Enhancer centrality - screen coverage (', analysis, ')')) +
  bimodal_theme +
  theme(legend.position  = 'bottom',
        legend.key.size  = unit(0.25, 'cm'),
        panel.grid.major = element_line(colour = 'grey92'))

save_panel_png_pdf(f6d, path_fig6, 'f6d.screen_coverage_scatter',
                   width_in = 3.4, height_in = 3.2)
