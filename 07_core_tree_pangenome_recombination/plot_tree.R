# ANALYSIS: Core genome tree, pangenome and recombination
# MANUSCRIPT SECTION: Phylogenetic reconstructions; Comparative Genomics Analyses
# USED FOR: Fig. 5, Fig. 6, Figs. S5.1 to S5.6
# FOLDER: 07_core_tree_pangenome_recombination (see README.md in the parent folder)
suppressPackageStartupMessages({library(ape); library(ggplot2); library(ggtree)})
a <- commandArgs(TRUE); f <- a[1]; out <- a[2]; ttl <- a[3]; og <- if (length(a) >= 4 && a[4] != "none") a[4] else NA
tr <- read.tree(f)
if (!is.na(og) && og %in% tr$tip.label) tr <- root(tr, outgroup = og, resolve.root = TRUE) else if (requireNamespace("phangorn", quietly = TRUE)) tr <- phangorn::midpoint(tr)
fix <- c(CicWB_2022 = "CicWB-2022", ZT3_1 = "ZT3-1", OVSP_L01 = "OVSP-L01", OVSP_S01 = "OVSP-S01", A_laidlawii = "A. laidlawii")
tr <- ladderize(tr); tr$tip.label <- ifelse(tr$tip.label %in% names(fix), fix[tr$tip.label], tr$tip.label); n <- Ntip(tr); mx <- max(node.depth.edgelength(tr))
p <- ggtree(tr, size = 0.5, color = "black"); p$data$sup <- suppressWarnings(as.numeric(p$data$label))
p <- p + geom_tiplab(size = ifelse(n > 40, 2.6, 4), color = "black") +
  geom_text2(aes(subset = !isTip & !is.na(sup) & sup >= 70, label = sup), size = ifelse(n > 40, 2.2, 3), color = "black", hjust = 1.15, vjust = -0.4) +
  geom_treescale(fontsize = 4, linesize = 0.6, color = "black") + xlim_tree(mx * 1.5) +
  labs(title = ttl, caption = "Node labels: ultrafast bootstrap support of 70 or more. Scale bar: substitutions per site.") +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5, color = "black"), plot.caption = element_text(size = 10, color = "black", hjust = 0),
        plot.margin = margin(10, 15, 10, 10), plot.background = element_rect(fill = "white", color = NA))
h <- max(4, 0.26 * n + 1.5)
ggsave(paste0(out, ".png"), p, width = 9, height = h, dpi = 600, bg = "white")
ggsave(paste0(out, ".pdf"), p, width = 9, height = h, bg = "white")
ggsave(paste0(out, "_preview.png"), p, width = 9, height = h, dpi = 80, bg = "white")
