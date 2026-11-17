# ANALYSIS: Marker genes, orthology and phylogenetic trees
# MANUSCRIPT SECTION: Multilocus sequence comparison; Phylogenetic reconstructions
# USED FOR: Fig. S2.8, Figs. S4.1 to S4.6
# FOLDER: 04_markers_orthology_trees (see README.md in the parent folder)
suppressPackageStartupMessages({library(ggplot2); library(ggtree); library(ape)})
tr <- read.tree("tuf_tree.treefile"); tr <- root(tr, outgroup = "A_laidlawii", resolve.root = TRUE); tr <- ladderize(tr)
p <- ggtree(tr, size = .3); d <- p$data; d$support <- suppressWarnings(as.numeric(d$label))
p <- p + geom_tiplab(size = 2.2, colour = "black") +
  geom_nodepoint(data = d[!d$isTip & !is.na(d$support) & d$support >= 95, ], size = 0.9, colour = "#2166ac") +
  geom_treescale(width = 0.1, fontsize = 2.5, offset = -1) + xlim(0, max(d$x) * 1.35)
ggsave("tuf_tree_figure.png", p, width = 7.2, height = 12, dpi = 600, bg = "white")
ggsave("tuf_tree_figure.pdf", p, width = 7.2, height = 12, bg = "white")
ggsave("tuf_tree_figure_preview.png", p, width = 7.2, height = 12, dpi = 110, bg = "white")
