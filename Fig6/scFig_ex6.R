primary_folder <- '[primary_folder]'
path_fig6_data <- paste0(primary_folder, 'Fig6/data/')
path_fig6      <- paste0(primary_folder, 'Fig6/out/')
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
# fig ex6a

category_table <- read.csv(paste0(path_fig6_data, 'sweep_category_table.csv'))
category_colours <- c(fully_uncoupled = '#4477AA',
                      intermediate    = '#CCBB44',
                      fully_coupled   = '#EE6677')
category_labels  <- c(fully_uncoupled = 'Fully uncoupled (Jaccard = 0)',
                      intermediate    = 'Intermediate (0 < Jaccard < 1)',
                      fully_coupled   = 'Fully coupled (Jaccard = 1)')

category_table1=category_table[which(category_table$theme == "TF_enhancer" & category_table$sweep %in% c("adjpvalue_threshold","dev_ratio_threshold")),]
category_table1$category=factor(category_table1$category, levels=c("fully_uncoupled","intermediate","fully_coupled"))

ext_f6a=ggplot(category_table1, aes(x = as.character(sweep_value), y=pct, fill=category)) +
  facet_wrap(~sweep, ncol=2, scales="free")+
  geom_bar(linewidth=0.25, alpha=1,  width=0.8,  color="black", stat="identity") + 
  scale_fill_manual(values=category_colours, labels=category_labels)+
  geom_text(aes(label=n),position = position_stack(vjust = 0.5), size=1.8,color="black")+
  labs(x=NULL, y = "% of promoters", title="Uncoupled TF contribution to enhancer transcription and accessibility") +
  theme1
pdf(paste0(path_fig6,"ext_f6a.stacked_TF_enhancer_coupling.pdf"), width = 3.5, height = 2)
print(ext_f6a)
dev.off()

#============================
# fig ex6b

bin_df <- read.csv(paste0(path_fig6_data, 'distance_bin_without_tCRE.csv'))
bin_df$dist_bin <- factor(bin_df$dist_bin, levels = unique(bin_df$dist_bin))
bin_df$source_class <- factor(bin_df$source_class, levels = c("expressed enhancers","non-expressed enhancers"))

ext_f6b <- ggplot(bin_df %>% dplyr::filter(n_enhancers >= 10),
                  aes(x = dist_bin, y = y_mean,
                      colour = source_class, group = source_class)) +
  geom_line(linewidth = 0.4) +
  geom_point(aes(size = n_enhancers)) +
  scale_colour_manual(name = 'Enhancer type', labels = c("Transcribed","Non-transcribed"), values = c('darkorange','steelblue')) +
  scale_size_continuous(name = 'N enhancers', range = c(0.5, 2.5)) +
  scale_x_discrete(name = 'Enhancer-promoter distance bin (kb)') +
  scale_y_continuous(name = 'Mean elastic-net-selected edges per enhancer') +
  labs(title="Promoter connections per enhancer across EP distance (no TSS signal)") +
  theme1 + theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))

pdf(paste0(path_fig6,"ext_f6b.edges_without_TSS_distance_bin.pdf"), width = 3.6, height = 2)
print(ext_f6b)
dev.off()

#============================
# fig ex6c

triplets_tCRE      <- read.csv(paste0(path_fig6_data, 'triplets_tCRE.csv'))
new_route_triplets <- triplets_tCRE[triplets_tCRE$new_route, ]

message('New-route triplets: ', nrow(new_route_triplets))

#### Top 20 TFs by new-route triplet count ####

tf_counts <- sort(table(new_route_triplets$TF), decreasing = TRUE)
n_tf      <- min(20, length(tf_counts))
tf_bar_df <- data.frame(
  feature    = factor(names(tf_counts)[seq_len(n_tf)],
                      levels = rev(names(tf_counts)[seq_len(n_tf)])),
  count = as.integer(tf_counts)[seq_len(n_tf)],
  group = "Top 20 TFs"
)

gene_counts <- sort(table(new_route_triplets$gene), decreasing = TRUE)
n_gene      <- min(20, length(gene_counts))
gene_bar_df <- data.frame(
  feature  = factor(names(gene_counts)[seq_len(n_gene)],
                    levels = rev(names(gene_counts)[seq_len(n_gene)])),
  count = as.integer(gene_counts)[seq_len(n_gene)],
  group = "Top 20 target genes"
)

all_df=rbind(tf_bar_df,gene_bar_df)
all_df$group=factor(all_df$group, levels=c("Top 20 TFs","Top 20 target genes"))

ext_f6c <- ggplot(all_df, aes(x = feature, y = count, fill=group)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~group, ncol=2, scales="free") +
  scale_fill_manual(values=c("steelblue","darkorange"), guide = NULL) +
  coord_flip() +
  scale_y_continuous(name = 'Number of new-route TSS triplets',
                     expand = expansion(mult = c(0, 0.08))) +
  scale_x_discrete(name = '') +
  ggtitle('TFs & target genes from the new enhancer TSS routes') +
  theme1 

pdf(paste0(path_fig6,"ext_f6c.new_route_triplet_bars_TF_and_target.pdf"), width = 3.9, height = 2.7)
print(ext_f6c)
dev.off()

#============================
# fig ex6d

gsea_df <- read.csv(paste0(path_fig6_data, 'gsea_results.csv'), stringsAsFactors = FALSE)
curves  <- read.csv(paste0(path_fig6_data, 'gsea_running_curves.csv'), stringsAsFactors = FALSE)
hits_df <- read.csv(paste0(path_fig6_data, 'gsea_hit_positions.csv'), stringsAsFactors = FALSE)

axis_labels  <- c(triplet_pairs = 'Rank (triplet-pair count)',
                  cum_abs_estimate_total = 'Rank (cumulative |estimate|)')
panel_labels <- c(expressed_TRUE = 'Transcribed enhancers',
                  mediating_TRUE = 'Mediating enhancers')
plot_sets <- c('expressed_TRUE', 'mediating_TRUE')
axes_names <- c('triplet_pairs', 'cum_abs_estimate_total')

make_gsea_cell <- function(axis_name, gs_name) {
  
  row <- gsea_df[gsea_df$axis == axis_name & gsea_df$gene_set == gs_name, ]
  fdr_label <- if (nrow(row) == 0) {
    'FDR=NA'
  } else if (!row$meaningful) {
    'FDR=N/A (trivial)'
  } else if (!is.na(row$fdr)) {
    paste0('FDR=', scientific(signif(row$fdr, 2)))
  } else 'FDR=NA'
  nes_label <- if (nrow(row) > 0 && !is.na(row$NES))
    paste0('NES=', round(row$NES, 2)) else ''
  annotation_colour <- if (nrow(row) > 0 && !is.na(row$fdr) && row$fdr < 0.05)
    'red3' else 'grey30'
  
  cv <- curves[curves$axis == axis_name & curves$gene_set == gs_name, ]
  hp <- hits_df[hits_df$axis == axis_name & hits_df$gene_set == gs_name, ]
  
  if (nrow(cv) == 0) return(patchwork::plot_spacer())
  
  col_curve <- if (nrow(row) > 0 && !is.na(row$NES) && row$NES > 0) '#B2182B' else '#2166AC'
  ylim_range <- range(cv$es)
  if (diff(ylim_range) < 0.01) ylim_range <- c(-0.05, 0.05)
  rug_y <- ylim_range[1] - 0.05 * diff(ylim_range)
  
  p <- ggplot(cv, aes(x = rank, y = es)) +
    geom_hline(yintercept = 0, colour = 'grey60', linetype = 'dashed', linewidth = 0.2) +
    geom_line(colour = col_curve, linewidth = 0.3) +
    { if (nrow(hp) > 0)
      geom_segment(data = hp, aes(x = rank, xend = rank, y = rug_y, yend = ylim_range[1]),
                   colour = col_curve, linewidth = 0.2, inherit.aes = FALSE) } +
    annotate('text', x = -Inf, y = Inf, label = paste(nes_label, fdr_label, sep = '  '),
             hjust = -0.05, vjust = 1.4, size = 1.8, colour = annotation_colour) +
    scale_y_continuous(name = 'Enrichment score',
                       limits = c(rug_y, ylim_range[2] * 1.05)) +
    scale_x_continuous(name = axis_labels[axis_name]) +
    ggtitle(panel_labels[gs_name]) +
    theme1 + theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
}

cells <- lapply(axes_names, function(axis_name)
  lapply(plot_sets, function(gs_name) make_gsea_cell(axis_name, gs_name)))

ext_f6d <- (cells[[1]][[1]] | cells[[2]][[1]] | cells[[2]][[2]])

pdf(paste0(path_fig6,"ext_f6d.gsea_enrichment_plots.pdf"), width = 3.3, height = 1.7)
print(ext_f6d)
dev.off()

#============================
# fig ex6e

plot_df  <- read.csv(paste0(path_fig6_data, 'hitscore_sign_enhancers.csv'), stringsAsFactors = FALSE)
tests_df <- read.csv(paste0(path_fig6_data, 'hitscore_sign_tests.csv'), stringsAsFactors = FALSE)

contrast_order <- c('iPSvsNSC_P0', 'NSC_P0vsNSC_P2', 'NSC_P2vsNeuron')
contrast_labels <- c(
  iPSvsNSC_P0    = 'iPS vs NSC P0',
  NSC_P0vsNSC_P2 = 'NSC P0 vs NSC P2',
  NSC_P2vsNeuron = 'NSC P2 vs Neuron'
)

plot_df$contrast_f <- factor(dplyr::recode(plot_df$contrast, !!!contrast_labels),
                             levels = unname(contrast_labels[contrast_order]))

bar_df <- plot_df %>%
  group_by(contrast_f, hit_sign, has_triplet) %>%
  summarise(n = n(), .groups = 'drop')

# Fisher-test FDR annotation per contrast (from precomputed hitscore_sign_tests.csv)
fisher_tests <- tests_df[tests_df$test == 'Fisher: triplet_pairs>0 x hit_score sign', ]
fisher_tests$contrast_f <- factor(dplyr::recode(fisher_tests$contrast, !!!contrast_labels),
                                  levels = unname(contrast_labels[contrast_order]))
fdr_label_df <- fisher_tests %>%
  mutate(label = paste0('FDR=', signif(fdr, 2)))

ext_f6e <- ggplot(bar_df, aes(x = hit_sign, y = n, fill = has_triplet)) +
  geom_col(position = 'fill', width=0.6) +
  facet_wrap(~contrast_f, nrow = 1) +
  geom_text(data = fdr_label_df, aes(x = 1.5, y = 1.1, label = label),
            inherit.aes = FALSE, size = 1.8) +
  scale_fill_manual(values = c(`TRUE` = '#B2182B', `FALSE` = 'grey80'),
                    name = 'Triplet participation',
                    labels = c(`TRUE` = 'triplet_pairs > 0', `FALSE` = 'triplet_pairs = 0')) +
  scale_y_continuous(name = 'Proportion', expand = expansion(mult = c(0, 0.14))) +
  coord_cartesian(ylim = c(0, 1.12)) +
  labs(title="Association between enhancer involvement and CRISPRi score") +
  labs(x = 'Hit-score sign') +
  theme1 +
  theme(legend.position = 'bottom',
        legend.key.size = unit(0.25, 'cm'))

pdf(paste0(path_fig6,"ext_f6e.hitscore_sign_triplet_bar.pdf"), width = 3.3, height = 1.25)
print(ext_f6e)
dev.off()

#===============================================================================
