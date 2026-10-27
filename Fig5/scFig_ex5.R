primary_folder <- '[primary_folder]'
path_fig5_data <- paste0(primary_folder, 'Fig5/data/')
path_fig5      <- paste0(primary_folder, 'Fig5/out/')
path_code      <- paste0(primary_folder, 'Data_and_code/')

library(ggplot2)
library(dplyr)

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
#source(paste0(path_code, 'network_loading/plot_theme_and_save.R'))

res <- read.csv(paste0(path_fig5_data, 'chipseq_f1_either_all.csv'))

network_labels <- c(
  base                   = 'Base, final',
  lasso                  = 'Elastic net, final',
  base_no_tobias         = 'Base, no TOBIAS',
  lasso_no_tobias        = 'Elastic net, no TOBIAS',
  base_no_hic            = 'Base, no Hi-C',
  lasso_no_hic           = 'Elastic net, no Hi-C',
  base_no_hic_no_tobias  = 'Base, no Hi-C, no TOBIAS',
  lasso_no_hic_no_tobias = 'Elastic net, no Hi-C, no TOBIAS'
)
network_levels <- unname(network_labels)
networks_colors <- c(
  base                   = '#6BAED6',
  base_no_tobias         = '#FC8D59',
  base_no_hic            = '#74C476',
  base_no_hic_no_tobias  = '#BCBDDC',
  lasso                  = '#08519C',
  lasso_no_tobias        = '#B30000',
  lasso_no_hic           = '#238B45',
  lasso_no_hic_no_tobias = '#6A51A3'
)
display_colors <- setNames(networks_colors, network_labels[names(networks_colors)])

#==================
# ex5a
res <- res %>%
  mutate(network = factor(dplyr::recode(network, !!!network_labels), levels = network_levels))

res$TF=factor(res$TF, levels=c("TEAD4","ONECUT2"))
ext_f5a <- ggplot(res, aes(x = dev_ratio, y = value, color = network)) +
  geom_line(linewidth = 0.3) +
  labs(title="TF-CRE interaction recovery from ChIP-seq data") +
  facet_wrap(~TF, scales = 'free') +
  scale_x_continuous(name = 'Deviance ratio threshold',
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)) +
  scale_y_continuous(name = 'F1 score') +
  scale_color_manual(name = 'Network', values = display_colors) +
  theme1 +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))
pdf(paste0(path_fig5,"ext_f5a.chipseq_F1_either_all2.pdf"), width = 4, height = 2)
print(ext_f5a)
dev.off()

#===============
# ex5b
results <- read.csv(paste0(path_fig5_data, 'eqtl_f1_all_interactions.csv'))

network_linetypes <- setNames(
  ifelse(grepl('no_hic', names(network_labels)), 'dashed', 'solid'),
  network_labels)

results <- results %>%
  mutate(network = factor(dplyr::recode(network, !!!network_labels), levels = network_levels),
         hic = factor(ifelse(grepl('no Hi-C', network), 'No Hi-C', 'With Hi-C'),
                      levels = c('With Hi-C', 'No Hi-C')))

results=results[-grep("TOBIAS",results$network ),] # remove redundancy
ext_f5b <- ggplot(results, aes(x = dev_ratio, y = value, color = network, linetype = hic)) +
  geom_line(linewidth = 0.3) +
  scale_x_continuous(name = 'Deviance ratio threshold',
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)) +
  scale_y_continuous(name = 'F1 score') +
  scale_color_manual(values = display_colors, guide = 'none') +
  scale_linetype_manual(name = 'Hi-C', values = c('With Hi-C' = 'solid', 'No Hi-C' = 'dashed')) +
  labs(title="Enhancer-promoter interaction recovery from eQTL data") +
  theme1 +
  theme(legend.position = 'bottom',
        legend.key.size = unit(0.25, 'cm'),
        legend.key.width = unit(0.5, 'cm'),
        legend.direction = 'horizontal',
        axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))
pdf(paste0(path_fig5,"ext_f5b.eqtl_F1_all_interactions2.pdf"), width = 3, height = 2)
print(ext_f5b)
dev.off()

#===============
# ex5c
d <- read.csv(paste0(path_fig5_data, 'network_node_counts.csv'))
d2 <- read.csv(paste0(path_fig5_data, 'network_edge_counts.csv'))
d$group="Network node counts"
d2$group="Network edge counts"
d=rbind(d, d2)
d$group=factor(d$group, levels=c("Network edge counts","Network node counts"))
d$network[which(d$network=="Base")]="Base, final"
d$network[which(d$network=="Elastic net")]="Elastic net, final"
d$metric=gsub("Number of ","",d$metric)
d$metric[which(d$metric=="edges")]="All edges"
d$metric=gsub(" edges","",d$metric)

#node_levels <- c('Number of TFs', 'Number of enhancers', 'Number of promoters')
raw_levels <- c('TFs', 'enhancers', 'promoters',
                'All', 'TF-to-peak', 'TF-to-promoter', 'TF-to-enhancer', 'enhancer-to-promoter')
d$metric  <- factor(d$metric,  levels = raw_levels)
d$network <- factor(d$network, levels = network_levels)

display_colors <- setNames(d$colour, d$network)[network_levels]
display_colors <- display_colors[!is.na(display_colors)]

ext_f5c <- ggplot(d, aes(x = network, y = value, fill = network)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~group+metric, nrow=2, ncol=5, scales = 'free') +
  scale_fill_manual(values = display_colors) +
  ggtitle('Content of 8 GRN variants') +
  theme1 +
  theme(legend.position = c(0.8,0.25),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

pdf(paste0(path_fig5,"ext_f5c.network_all_counts.pdf"), width = 6.8, height = 2.5)
print(ext_f5c)
dev.off()

#===============
# ex5d

tier_comparison <- read.csv(paste0(path_fig5_data, 'tier_comparison.csv'), stringsAsFactors = FALSE)

tier_levels <- c('base', 'elastic_net', 'filtered')
tier_labels <- c(base = 'Base', elastic_net = 'Elastic\nnet', filtered = 'High-\nconfidence')

variant_labels <- c(iPSC = 'iPSC', NSC = 'NSC', Neuron = 'Neuron',
                    max_across_celltypes = 'Max across cell types')
variant_colors <- c(iPSC = '#6BAED6', NSC = '#74C476', Neuron = '#FC8D59',
                    `max_across_celltypes` = '#252525')

tier_comparison <- tier_comparison %>%
  mutate(tier = factor(dplyr::recode(tier, !!!tier_labels), levels = unname(tier_labels)),
         abc_variant = factor(dplyr::recode(abc_variant, !!!variant_labels),
                              levels = unname(variant_labels)))
display_colors <- setNames(variant_colors, variant_labels[names(variant_colors)])

ext_f5d <- ggplot(tier_comparison, aes(x = tier, y = f1, colour = abc_variant, group = abc_variant)) +
  geom_line(linewidth = 0.3) +
  geom_point(size = 1.2) +
  scale_x_discrete(name = 'Network tier') +
  scale_y_continuous(name = 'F1 (ABC positivity vs network edge)') +
  labs(title="Comparison with ABC model") +
  scale_colour_manual(name = 'ABC variant', values = display_colors) +
  theme1
pdf(paste0(path_fig5,"ext_f5d.abc_tier_comparison_f12.pdf"), width = 2.6, height = 1.8)
print(ext_f5d)
dev.off()

#===============
# ex5e

category_table <- read.csv(paste0(path_fig5_data, 'sweep_category_table.csv'))
category_colours <- c(fully_uncoupled = '#4477AA',
                      intermediate    = '#CCBB44',
                      fully_coupled   = '#EE6677')
category_labels  <- c(fully_uncoupled = 'Fully uncoupled (Jaccard = 0)',
                      intermediate    = 'Intermediate (0 < Jaccard < 1)',
                      fully_coupled   = 'Fully coupled (Jaccard = 1)')

category_table1=category_table[which(category_table$theme == "TF_promoter" & category_table$sweep %in% c("adjpvalue_threshold","dev_ratio_threshold")),]
category_table1$category=factor(category_table1$category, levels=c("fully_uncoupled","intermediate","fully_coupled"))

ext_f5e=ggplot(category_table1, aes(x = as.character(sweep_value), y=pct, fill=category)) +
  facet_wrap(~sweep, ncol=2, scales="free")+
  geom_bar(linewidth=0.25, alpha=1,  width=0.8,  color="black", stat="identity") + 
  scale_fill_manual(values=category_colours, labels=category_labels)+
  geom_text(aes(label=n),position = position_stack(vjust = 0.5), size=1.8,color="black")+
  labs(x=NULL, y = "% of promoters", title="Uncoupled TF contribution to CRE transcription and accessibility") +
  theme1
pdf(paste0(path_fig5,"ext_f5e.stacked_TF_promoter_coupling.pdf"), width = 4.5, height = 1.8)
print(ext_f5e)
dev.off()

#===============================================================================



