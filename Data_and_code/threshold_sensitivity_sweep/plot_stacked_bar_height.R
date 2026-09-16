#### Shared plotting helper for the six threshold-sensitivity stacked bars ####
#
# Sourced by:
#   Fig5/fig5d_stacked_promoter_pctl.R          (TF_promoter x estimate_percentile)
#   Fig5/ext_fig5d_stacked_promoter_adjp.R      (TF_promoter x adjpvalue_threshold)
#   Fig5/ext_fig5e_stacked_promoter_devratio.R  (TF_promoter x dev_ratio_threshold)
#   Fig6/fig6a_stacked_enhancer_pctl.R          (TF_enhancer x estimate_percentile)
#   Fig6/ext_fig6a_stacked_enhancer_adjp.R      (TF_enhancer x adjpvalue_threshold)
#   Fig6/ext_fig6b_stacked_enhancer_devratio.R  (TF_enhancer x dev_ratio_threshold)
#
# The plotting logic is identical across all six panels -- only theme_name and
# sweep_name differ -- so it lives here rather than being copied six times.
#
# Draws the linear-height variant: total bar height is the raw n_nodes of
# that grid point, so the shrinking sample size at stringent thresholds is
# visible directly in the figure (a log10-height variant would flatten out
# that shrinkage and is not used here).
#
# A grid point with n_nodes = 0 renders as an empty bar labelled "n = 0" so
# every tested threshold value stays on the x-axis. No grid point in the
# current sweep hits this, but the code path is kept for robustness.
#
# Requires bimodal_theme + save_panel_png_pdf from
# Data_and_code/network_loading/plot_theme_and_save.R (sourced by the caller).

library(ggplot2)
library(dplyr)

# Category colours / labels used consistently across all six panels.
category_colours <- c(fully_uncoupled = '#4477AA',
                      intermediate    = '#CCBB44',
                      fully_coupled   = '#EE6677')
category_labels  <- c(fully_uncoupled = 'Fully uncoupled (Jaccard = 0)',
                      intermediate    = 'Intermediate (0 < Jaccard < 1)',
                      fully_coupled   = 'Fully coupled (Jaccard = 1)')

sweep_x_labels <- c(adjpvalue_threshold = 'adjpvalue threshold (<=)',
                    dev_ratio_threshold = 'dev_ratio threshold (>=)',
                    estimate_percentile = '|estimate| percentile per stratum')

# category_table: as written by build_sweep_data.R (sweep, sweep_value, theme,
#                 category, pct, n, n_nodes)
# theme_name    : 'TF_enhancer' or 'TF_promoter'
# sweep_name    : one of names(sweep_x_labels)
make_stacked_bar_linheight <- function(category_table, theme_name, sweep_name) {

  all_sweep_values <- sort(unique(category_table$sweep_value[category_table$sweep == sweep_name]))

  n_lookup <- category_table %>%
    filter(sweep == sweep_name, theme == theme_name) %>%
    distinct(sweep_value, n_nodes) %>%
    mutate(height_n = n_nodes)          # 'linear' height transform

  zero_n_values <- n_lookup$sweep_value[n_lookup$n_nodes == 0]
  max_height_n  <- max(n_lookup$height_n, na.rm = TRUE)

  df <- category_table %>%
    filter(theme == theme_name, sweep == sweep_name) %>%
    left_join(n_lookup, by = c('sweep_value', 'n_nodes')) %>%
    mutate(category    = factor(category, levels = names(category_colours)),
           sweep_value = factor(sweep_value, levels = all_sweep_values),
           bar_height  = pct / 100 * height_n)

  label_df <- if (length(zero_n_values) > 0) {
    data.frame(sweep_value = factor(zero_n_values, levels = all_sweep_values),
               bar_height  = max_height_n * 0.015)
  } else {
    data.frame(sweep_value = factor(character(0), levels = all_sweep_values),
               bar_height  = numeric(0))
  }

  ggplot(df, aes(x = sweep_value, y = bar_height, fill = category)) +
    geom_col(width = 0.7) +
    { if (nrow(label_df) > 0)
        geom_text(data = label_df, aes(x = sweep_value, y = bar_height, label = 'n = 0'),
                  inherit.aes = FALSE, size = 1.8, colour = 'grey40', vjust = -0.5,
                  family = 'Arial') } +
    scale_fill_manual(name = '', values = category_colours, labels = category_labels) +
    scale_y_continuous(name = 'N targets',
                       limits = c(0, max_height_n * 1.08), expand = c(0, 0)) +
    scale_x_discrete(name = sweep_x_labels[[sweep_name]], drop = FALSE) +
    ggtitle(paste0(theme_name, ': TF-set coupling vs ', sweep_x_labels[[sweep_name]]),
            subtitle = 'bar height = N targets') +
    bimodal_theme +
    theme(legend.position     = 'bottom',
          legend.key.size     = unit(0.25, 'cm'),
          legend.direction    = 'vertical',
          panel.grid.major.x  = element_blank())
}
