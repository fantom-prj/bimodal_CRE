primary_folder <- '[primary_folder]'
path_extfig7_data <- paste0(primary_folder, 'ExtFig7/data/')
path_extfig7       <- paste0(primary_folder, 'ExtFig7/out/')
path_code          <- paste0(primary_folder, 'Data_and_code/')

library(ggplot2)
library(dplyr)
library(ggraph)
library(igraph)
library(patchwork)


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
# fig ex7a
significance_threshold <- 0.05
res_collapsed <- read.csv(paste0(path_extfig7_data, 'ep_enrichment_collapsed.csv'),
                          stringsAsFactors = FALSE)

feature_display <- c(
  any.27ac_enhancer  = 'H3K27ac',   any.27ac_promoter  = 'H3K27ac',
  any.27m3_enhancer  = 'H3K27me3',  any.27m3_promoter  = 'H3K27me3',
  iPSC_SE_enhancer   = 'iPSC SE',   iPSC_SE_promoter   = 'iPSC SE',
  NSC_SE_enhancer    = 'NSC SE',    NSC_SE_promoter    = 'NSC SE',
  Neuron_SE_enhancer = 'Neuron SE', Neuron_SE_promoter = 'Neuron SE',
  CGI_enhancer       = 'CGI',       CGI_promoter       = 'CGI',
  TATA_box_enhancer  = 'TATA box',  TATA_box_promoter  = 'TATA box'
)
feature_levels  <- c('H3K27ac', 'H3K27me3', 'iPSC SE', 'NSC SE', 'Neuron SE', 'CGI', 'TATA box')
stratum_display <- c(all = 'All enhancers', tCRE_all = 'Transcribed enhancers',
                     mediating = 'Mediating enhancers')

min_pval <- min(res_collapsed$adjpvalue[res_collapsed$adjpvalue > 0], na.rm = TRUE)

res_plot <- res_collapsed %>%
  mutate(
    adjpvalue   = ifelse(adjpvalue == 0, min_pval / 2, adjpvalue),
    log10p      = -log10(adjpvalue),
    significant = adjpvalue <= significance_threshold & !is.na(adjpvalue),
    odds_ratio  = ifelse(!significant, NA_real_, odds_ratio),
    log10p      = ifelse(!significant, NA_real_, log10p),
    enhancer    = factor(dplyr::recode(enhancer, !!!feature_display), levels = feature_levels),
    promoter    = factor(dplyr::recode(promoter, !!!feature_display), levels = rev(feature_levels)),
    stratum     = factor(dplyr::recode(stratum, !!!stratum_display),
                         levels = unname(stratum_display))
  )

ext_f7a <- ggplot(res_plot, aes(x = enhancer, y = promoter)) +
  geom_point(aes(size = log10p, fill = odds_ratio, colour = significant),
             shape = 21, stroke = 0.5) +
  geom_text(aes(label = ifelse(significant, round(odds_ratio, 2), '')),
            size = 1.6, colour = 'black', fontface = 'bold',
            nudge_y = 0.35, vjust = 0) +
  scale_fill_gradient2(name = 'Odds ratio',
                       low = 'blue', mid = 'white', high = 'red',
                       midpoint = 1, na.value = 'grey92') +
  scale_colour_manual(name = 'Significant',
                      values = c('FALSE' = 'grey80', 'TRUE' = 'black'),
                      guide  = 'none') +
  scale_size_continuous(name = '-log10(adj. p)', range = c(1, 5)) +
  facet_wrap(~stratum) +
  theme1 +
  theme(
    axis.text.x      = element_text(angle = 25, hjust = 1, vjust = 1),
    panel.grid.major = element_line(colour = 'grey92'),
    legend.key.size  = unit(0.25, 'cm')
  ) +
  labs(x = 'Enhancer feature', y = 'Promoter feature', title="Enrichment of chromatin feature between connected enhancer-promoter")

pdf(paste0(path_extfig7,"ext_f7a.ep_enrichment_collapsed2.pdf"), width = 7, height = 2.5)
print(ext_f7a)
dev.off()

#==========================
# fig ex7b

d <- readRDS(paste0(path_extfig7_data, 'tf_umap_dominance.rds'))

layout_all        <- d$layout_all
edge_direct       <- d$edge_direct
edge_indirect     <- d$edge_indirect
edge_combined     <- d$edge_combined
node_names        <- d$node_names
umap_fixed        <- d$umap_fixed
num_label_per_cat <- d$num_label_per_cat
dom_colours       <- d$dom_colours

# Labels: top num_label_per_cat TFs per dominance category by combined PageRank
tfs_to_label_all <- unlist(lapply(
  c('Direct-dominant', 'Indirect-dominant', 'Mixed'),
  function(cat) {
    idx <- which(layout_all$dominance == cat)
    if (length(idx) == 0) return(character(0))
    idx_top <- idx[order(layout_all$pr[idx], decreasing = TRUE)[
      seq_len(min(num_label_per_cat, length(idx)))]]
    node_names[idx_top]
  }
))


make_dom_panel <- function(layout, edge_df_panel, tfs_label, title_str) {
  edge_lims <- if (nrow(edge_df_panel) > 0)
    c(quantile(edge_df_panel$value, 0.1), quantile(edge_df_panel$value, 0.9))
  else c(0, 1)
  g_panel <- igraph::graph_from_data_frame(edge_df_panel,
                                           vertices = data.frame(name = node_names))
  lay <- ggraph::create_layout(g_panel, layout = umap_fixed[node_names, ])
  lay$pr        <- layout$pr
  lay$dominance <- layout$dominance
  ggraph(graph = lay) +
    geom_edge_diagonal(aes(width = value), color = 'grey75', alpha = 0.4) +
    scale_edge_width_continuous(range = c(0.01, 0.3), limits = edge_lims, guide = 'none') +
    geom_node_point(aes(size = pr, colour = dominance), alpha = 0.85) +
    geom_node_text(aes(label = name, filter = name %in% tfs_label),
                   size = 1.8, repel = TRUE) +
    scale_size_continuous(name = 'Reverse\nPageRank', range = c(0.5, 3)) +
    scale_colour_manual(name = 'Dominance category', values = dom_colours, drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 2.5))) +
    theme_graph(background = 'white',  base_size = 6, base_family="Helvetica") +
    theme(legend.text  = element_text(size = 5),
          legend.title = element_text(size = 6),
          plot.title   = element_text(size = 7, face="plain"),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    labs(title = title_str)
}

p_panel_direct   <- make_dom_panel(layout_all, edge_direct,   tfs_to_label_all,
                                   'Direct edges only\n(TF -> promoter(TF))')
p_panel_indirect <- make_dom_panel(layout_all, edge_indirect, tfs_to_label_all,
                                   'Indirect edges only\n(TF -> enhancer -> promoter(TF))')
p_panel_combined <- make_dom_panel(layout_all, edge_combined, tfs_to_label_all,
                                   'Combined edges\n(stronger of direct/indirect per pair)')

ext_f7b <- (p_panel_direct | p_panel_indirect | p_panel_combined) +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom', legend.box = 'vertical')

pdf(paste0(path_extfig7,"ext_f7b.tf_umap_dominance2.pdf"), width = 7.2, height = 3.5)
print(ext_f7b)
dev.off()

#===============================================================================


