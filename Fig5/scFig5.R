library(ComplexHeatmap)
library(circlize)
library(grid)
library(tidyverse)

primary_folder <- '[primary_folder]'
primary_folder <- "/analysisdata/fantom6/Interactome/single_cell_wallace/Figure/"
path_fig5_data <- paste0(primary_folder, 'Fig5/data/')
path_fig5      <- paste0(primary_folder, 'Fig5/out/')
path_code      <- paste0(primary_folder, 'Data_and_code/')

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
source(paste0(path_code, 'network_loading/plot_theme_and_save.R'))

# The right-side logRP (log rank-product) bar restates, as a bar chart, the
# same overall-importance ranking that the row order (top = highest rank
# within each cluster) already conveys -- it was found to read as visual
# clutter next to the heatmap body without adding information not already
# implicit in row order, so it is now OFF by default. Set to TRUE (e.g.
# `include_logrp_annotation <- TRUE` before sourcing) to bring it back.
if (!exists('include_logrp_annotation')) include_logrp_annotation <- FALSE

#===============================
#f5b
d <- readRDS(paste0(path_fig5_data, 'tf_effects_clustering.rds'))

sh             <- d$sh
sh_sqrt        <- d$sh_sqrt
rk             <- d$rk
row_group_fine <- d$row_group_fine
col_order      <- d$col_order
cdend          <- d$cdend
stat_label     <- d$stat_label

sqrt_signed <- function(x) sqrt(abs(x)) * sign(x)

# colour scale on the scaled values, centred at 0
col_fun <- colorRamp2(
  c(quantile(sh, 0.025, na.rm = TRUE), 0, quantile(sh, 0.975, na.rm = TRUE)),
  c('blue', 'white', 'red')
)

# column annotation: measurement type (ATAC/TSS) nested under target class
col_anno_df <- data.frame(
  target      = c('Mediating Enhancer', 'Mediating Enhancer', 'Promoter', 'Promoter'),
  measurement = c('ATAC', 'TSS', 'ATAC', 'TSS'),
  row.names   = col_order
)
col_anno_colors <- list(
  target      = c('Promoter' = '#4DAF4A', 'Mediating Enhancer' = '#984EA3'),
  measurement = c(ATAC = '#377EB8', TSS = '#E41A1C')
)
ha_col <- HeatmapAnnotation(
  measurement = col_anno_df$measurement,
  target      = col_anno_df$target,
  col         = col_anno_colors,
  annotation_name_side  = 'left',
  annotation_name_gp    = gpar(fontfamily = 'Arial', fontsize = 6),
  annotation_legend_param = list(
    measurement = list(title_gp  = gpar(fontfamily = 'Arial', fontsize = 6),
                       labels_gp = gpar(fontfamily = 'Arial', fontsize = 5)),
    target      = list(title_gp  = gpar(fontfamily = 'Arial', fontsize = 6),
                       labels_gp = gpar(fontfamily = 'Arial', fontsize = 5))
  ),
  simple_anno_size = unit(2.5, 'mm')
)

# row annotation: log rank-product bar (optional, see include_logrp_annotation)
ha_row <- if (include_logrp_annotation) {
  rowAnnotation(
    logRP = anno_barplot(rk, border = FALSE,
                         gp = gpar(col = 'grey50', fill = 'grey50'),
                         axis_param = list(gp = gpar(fontfamily = 'Arial', fontsize = 5))),
    annotation_name_gp = gpar(fontfamily = 'Arial', fontsize = 6),
    width = unit(8, 'mm')
  )
} else {
  NULL
}

# Display sh (scaled); cluster on sh_sqrt (compressed) within each split group.
ht <- Heatmap(sh,
              name                     = paste0('Scaled TF effect\n[', stat_label, ']'),
              col                      = col_fun,
              row_split                = row_group_fine,
              cluster_columns          = cdend,
              cluster_rows             = TRUE,
              clustering_distance_rows = function(m) dist(sqrt_signed(m)),
              clustering_method_rows   = 'ward.D2',
              top_annotation           = ha_col,
              right_annotation         = ha_row,
              show_row_names           = FALSE,
              show_column_names        = FALSE,
              row_title_gp             = gpar(fontfamily = 'Arial', fontsize = 6),
              column_title_gp          = gpar(fontfamily = 'Arial', fontsize = 7),
              heatmap_legend_param     = list(
                title_gp     = gpar(fontfamily = 'Arial', fontsize = 6),
                labels_gp    = gpar(fontfamily = 'Arial', fontsize = 5),
                direction    = 'horizontal',
                grid_height  = unit(2.5, 'mm'),
                legend_width = unit(25, 'mm'),
                title_position = 'topcenter'
              ))

save_grid_png_pdf(function() draw(ht, heatmap_legend_side = 'bottom',
                                  annotation_legend_side = 'bottom',
                                  legend_grouping = 'original'),
                  path_fig5, 'f5b.heatmap_TF_effects2',
                  width_in = 3.4, height_in = 7)

#=====================
# f5c
significance_threshold <- 0.1   # FDR cutoff for the black significance outline

d <- readRDS(paste0(path_fig5_data, 'tf_effects_clustering.rds'))
h_wide      <- d$h_wide
rk          <- d$rk
row_group   <- d$row_group
results_bar <- d$results_bar

#### TF characteristics table ####

tf_characteristics <- read.table(
  paste0(path_fig3_data, 'TF_grouping_by_chromvar_correlation_and_literature.n258.tsv.gz'),
  header = TRUE, check.names = FALSE)

tf_characteristics <- tf_characteristics[
  !grepl('::', tf_characteristics$motifName, fixed = TRUE), ]
tf_characteristics$motifName <- toupper(tf_characteristics$motifName)
tf_characteristics <- unique(tf_characteristics[,c(36:39,41:44)])
colnames(tf_characteristics) <- gsub('+', 'plus',  colnames(tf_characteristics), fixed = TRUE)
colnames(tf_characteristics) <- gsub('-', 'minus', colnames(tf_characteristics), fixed = TRUE)

to_remove_motifs <- tf_characteristics$motifName[duplicated(tf_characteristics$motifName)]
tf_characteristics <- tf_characteristics[
  !(tf_characteristics$motifName %in% to_remove_motifs), ]

#### TF-level statistics table ####

all_peaks_wide <- results_bar[['All']] %>%
  dplyr::select(source, target_measurement_type, log_total_effect) %>%
  mutate(col = paste0('all_peaks_', target_measurement_type)) %>%
  dplyr::select(source, col, log_total_effect) %>%
  pivot_wider(names_from = col, values_from = log_total_effect) %>%
  column_to_rownames('source') %>%
  as.matrix()

common_tfs <- intersect(rownames(h_wide), rownames(all_peaks_wide))
tf_statistics <- data.frame(
  h_wide[common_tfs, ],
  all_peaks_wide[common_tfs, ],
  logRP                = rk[common_tfs],
  dominant_effect_sign = row_group[common_tfs],
  motifName            = toupper(common_tfs)
)

tf_df <- merge(tf_statistics, tf_characteristics, by = 'motifName', all = FALSE)
message('TFs matched to characteristics table: ', nrow(tf_df))

#### Welch t-test per (TF effect x TF characteristic) ####

cont_vars <- c('mediating_enhancer_aCRE', 'mediating_enhancer_tCRE',
               'promoter_aCRE', 'promoter_tCRE',
               'all_peaks_aCRE', 'all_peaks_tCRE')
tf_char_names <- setdiff(colnames(tf_characteristics), 'motifName')

char_results <- data.frame(
  expand.grid(TF_effect = cont_vars, TF_characteristic = tf_char_names),
  diff   = NA_real_,
  pvalue = NA_real_
)
for (iter in seq_len(nrow(char_results))) {
  ef  <- as.character(char_results$TF_effect[iter])
  ch  <- as.character(char_results$TF_characteristic[iter])
  tmp <- t.test(as.formula(paste0('`', ef, '` ~ `', ch, '`')), data = tf_df)
  char_results$diff[iter]   <- tmp$estimate[2] - tmp$estimate[1]
  char_results$pvalue[iter] <- tmp$p.value
}
char_results$qvalue      <- p.adjust(char_results$pvalue, method = 'fdr')
char_results$adjpval     <- -log10(char_results$qvalue)
char_results$significant <- char_results$qvalue <= significance_threshold

message('Significant (FDR <= ', significance_threshold, '): ',
        sum(char_results$significant), ' of ', nrow(char_results))

#### Dot plot ####

x_labels <- c(
  mediating_enhancer_aCRE = 'Mediating enh. (ATAC)',
  mediating_enhancer_tCRE = 'Mediating enh. (TSS)',
  promoter_aCRE           = 'Promoter (ATAC)',
  promoter_tCRE           = 'Promoter (TSS)',
  all_peaks_aCRE          = 'All peaks (ATAC)',
  all_peaks_tCRE          = 'All peaks (TSS)'
)
y_labels <- c(
  TF_activator = 'Activator',
  TF_repressor = 'Repressor',
  TF_pioneer   = 'Pioneer',
  #TF_closer    = 'Closer',
  ATACplus     = 'ATAC+',
  ATACminus    = 'ATAC-',
  TSSplus      = 'TSS+',
  TSSminus     = 'TSS-'
)

char_results$x_lab <- factor(x_labels[as.character(char_results$TF_effect)],
                             levels = x_labels)
char_results$y_lab <- factor(y_labels[as.character(char_results$TF_characteristic)],
                             levels = rev(y_labels))

f5c <- ggplot(char_results,
              aes(x = x_lab, y = y_lab,
                  fill = diff, colour = significant, size = adjpval)) +
  geom_point(shape = 21, stroke = 0.4) +
  scale_x_discrete(name = 'TF cumulative effect') +
  scale_y_discrete(name = 'TF annotation from ChromVar / literature') +
  scale_colour_manual(name = 'Significant',
                      values = c('FALSE' = 'grey80', 'TRUE' = 'black')) +
  scale_fill_gradient2(name = 'Average\nDifference',
                       low = 'blue', mid = 'white', high = 'red', midpoint = 0) +
  labs(title = " Associations of TF regulatory modes") +
  scale_size_continuous(name = '-log10(adj. p)', range = c(1, 4.5)) +
  theme1 +
  theme(axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.major = element_line(colour = 'grey92'),
        legend.key.size  = unit(0.3, 'cm'))

pdf(paste0(path_fig5,"f5c.TF_characteristics_dotplot2.pdf"), width = 3.4, height = 2.2)
print(f5c)
dev.off()

#=====================
# f5d
category_table <- read.csv(paste0(path_fig5_data, 'sweep_category_table.csv'))

category_colours <- c(fully_uncoupled = '#4477AA',
                      intermediate    = '#CCBB44',
                      fully_coupled   = '#EE6677')
category_labels  <- c(fully_uncoupled = 'Fully uncoupled (Jaccard = 0)',
                      intermediate    = 'Intermediate (0 < Jaccard < 1)',
                      fully_coupled   = 'Fully coupled (Jaccard = 1)')

category_table1=category_table[which(category_table$theme == "TF_promoter" & category_table$sweep=="estimate_percentile"),]
category_table1$category=factor(category_table1$category, levels=c("fully_uncoupled","intermediate","fully_coupled"))

f5d=ggplot(category_table1, aes(x = as.character(sweep_value), y=pct, fill=category)) +
  geom_bar(linewidth=0.25, alpha=1,  width=0.8,  color="black", stat="identity") + 
  scale_fill_manual(values=category_colours, labels=category_labels)+
  geom_text(aes(label=n),position = position_stack(vjust = 0.5), size=1.8,color="black")+
  labs(x="|estimate| percentile per stratum", y = "% of promoters", title="Uncoupled TF contribution to CRE transcription and accessibility") +
  theme1
pdf(paste0(path_fig5,"f5d.stacked_TF_promoter_estimate_pctl2.pdf"), width = 3.4, height = 2)
print(f5d)
dev.off()

#===============================================================================


