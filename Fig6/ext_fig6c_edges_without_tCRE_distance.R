#### Extended Figure 6c — edges per enhancer by distance bin, without-tCRE ####
#
# Panel:   Ext Fig 6c. Same quantity as Fig 6b, but computed on the
#          elastic-net network fitted WITHOUT enhancer tCRE measurements --
#          the equal-footing control: if expressed enhancers still look
#          enriched for elastic-net-selected edges when they can only
#          contribute aCRE edges
#          (i.e. cannot benefit from having a second, transcription-based
#          measurement), the Fig 6b effect is not merely "expressed
#          enhancers get a second chance at an edge".
# Reads:   Fig6/data/distance_bin_without_tCRE.csv
#          (written by ../Data_and_code/tCRE_topology/build_tCRE_topology_data.R)
# Writes:  Fig6/out/ext_f6c.edges_without_tCRE_distance_bin.{png,pdf}
#
# Bins with fewer than 10 enhancers are dropped to avoid displaying noisy
# per-enhancer averages from very small bins.

library(ggplot2)
library(dplyr)

primary_folder <- '[primary_folder]'
path_fig6_data <- paste0(primary_folder, 'Fig6/data/')
path_fig6      <- paste0(primary_folder, 'Fig6/out/')
path_code      <- paste0(primary_folder, 'Data_and_code/')

source(paste0(path_code, 'network_loading/plot_theme_and_save.R'))

bin_df <- read.csv(paste0(path_fig6_data, 'distance_bin_without_tCRE.csv'))
bin_df$dist_bin <- factor(bin_df$dist_bin, levels = unique(bin_df$dist_bin))

ext_f6c <- ggplot(bin_df %>% filter(n_enhancers >= 10),
                  aes(x = dist_bin, y = y_mean,
                      colour = source_class, group = source_class)) +
  geom_line(linewidth = 0.4) +
  geom_point(aes(size = n_enhancers)) +
  scale_colour_manual(name = '',
                      values = c('expressed enhancers'     = 'darkorange',
                                 'non-expressed enhancers' = 'steelblue')) +
  scale_size_continuous(name = 'N enhancers', range = c(0.5, 2.5)) +
  scale_x_discrete(name = 'Enhancer-promoter distance bin (kb)') +
  scale_y_continuous(name = 'Mean elastic-net-selected edges per enhancer') +
  # No title -- the panel letter + figure caption identify it.
  bimodal_theme +
  theme(axis.text.x      = element_text(angle = 45, hjust = 1),
        legend.position  = 'bottom',
        legend.box       = 'vertical',
        legend.key.size  = unit(0.25, 'cm'),
        panel.grid.major = element_line(colour = 'grey92'))

# Sized for the Ext Fig 6 layout (A|B)/C/(D|(E/F)): C has its own row but
# stays half-width (~3.66in), matching the A-D main-figure panel convention,
# per explicit user instruction, rather than stretching to the full A4 width
# that (A|B) uses.
save_panel_png_pdf(ext_f6c, path_fig6, 'ext_f6c.edges_without_tCRE_distance_bin',
                   width_in = 3.66, height_in = 2.6)
