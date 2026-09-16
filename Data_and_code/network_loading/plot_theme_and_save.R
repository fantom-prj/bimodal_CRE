#### Shared plotting helpers: Arial 5-7pt theme + PNG+PDF save ####
#
# Sourced by every panel script in Fig5/, Fig6/, ExtFig7/. Per co-author
# request: all text Arial, sized 5-7pt; every panel saved as both PNG and PDF
# at the same physical size.

library(ggplot2)

# Arial is resolved directly from the system font registry by ggsave's PNG
# device (via ragg/systemfonts) and by cairo_pdf() for the PDF device — no
# extrafont font-database step needed on a machine where Arial is installed
# (confirmed present as /System/Library/Fonts/Supplemental/Arial.ttf here).

bimodal_theme <- theme_bw(base_family = 'Arial', base_size = 6) +
  theme(
    text             = element_text(family = 'Arial', size = 6, colour = 'black'),
    axis.text        = element_text(family = 'Arial', size = 5, colour = 'black'),
    axis.title       = element_text(family = 'Arial', size = 6, colour = 'black'),
    plot.title       = element_text(family = 'Arial', size = 7, colour = 'black', hjust = 0.5),
    plot.subtitle    = element_text(family = 'Arial', size = 6, colour = 'black', hjust = 0.5),
    legend.text      = element_text(family = 'Arial', size = 5, colour = 'black'),
    legend.title     = element_text(family = 'Arial', size = 6, colour = 'black'),
    strip.text       = element_text(family = 'Arial', size = 6, colour = 'black'),
    panel.grid.minor = element_blank()
  )

# Save a ggplot object as both PNG and PDF at the same physical size.
# width_in / height_in are inches; PNG is rendered at 300 dpi.
save_panel_png_pdf <- function(plot, out_dir, base_name, width_in, height_in, dpi = 300) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(out_dir, paste0(base_name, '.png')), plot,
         width = width_in, height = height_in, units = 'in', dpi = dpi)
  ggsave(file.path(out_dir, paste0(base_name, '.pdf')), plot,
         width = width_in, height = height_in, units = 'in', device = cairo_pdf)
}

# Same, but for a pre-rendered grid/ComplexHeatmap object drawn via draw()/
# grid.draw() inside a graphics device (heatmaps, custom grid layouts) rather
# than a ggplot object ggsave() can serialize directly.
save_grid_png_pdf <- function(draw_fn, out_dir, base_name, width_in, height_in, dpi = 300) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  png(file.path(out_dir, paste0(base_name, '.png')),
      width = width_in, height = height_in, units = 'in', res = dpi, family = 'Arial')
  draw_fn()
  dev.off()
  cairo_pdf(file.path(out_dir, paste0(base_name, '.pdf')),
            width = width_in, height = height_in, family = 'Arial')
  draw_fn()
  dev.off()
}
