# ANALYSIS: Marker genes, orthology and phylogenetic trees
# MANUSCRIPT SECTION: Multilocus sequence comparison; Phylogenetic reconstructions
# USED FOR: Fig. S2.8, Figs. S4.1 to S4.6
# FOLDER: 04_markers_orthology_trees (see README.md in the parent folder)
suppressMessages({library(ape); library(ggtree); library(ggplot2)})
norm <- function(x) gsub("-", "_", x)
PPS <- c("PPS1A","PPS2B","PPS6A","PPS6B","PPS9B","PPS10A")
cl <- list("#4C78A8" = c(PPS, "CanS4","MW1","VAC","MA1"), "#F58518" = c("CX","6A1","ZT3_1","PR2021","JR1","PDA15"),
           "#54A24B" = c("CicWB_2022","ChTDIII","OVSP_L01","OVSP_S01","Vc33"), "#B279A2" = "LWB")
colmap <- unlist(lapply(names(cl), function(k) setNames(rep(k, length(cl[[k]])), cl[[k]])))
getsup <- function(t) { if (is.null(t$node.label)) return(rep(NA_real_, t$Nnode))
  vapply(strsplit(t$node.label, "/"), function(v) if (length(v) == 0) NA_real_ else suppressWarnings(as.numeric(tail(v, 1))), numeric(1)) }
draw <- function(f, out, bl = FALSE) {
  t <- read.tree(f); t$tip.label <- norm(t$tip.label)
  og <- grep("laidlawii", t$tip.label, value = TRUE)
  if (length(og) == 1) t <- root(t, outgroup = og, resolve.root = TRUE)
  else if (requireNamespace("phangorn", quietly = TRUE)) t <- phangorn::midpoint(t)
  else cat("phangorn missing: not midpoint rooted:", f, "\n")
  sup <- getsup(t); n <- Ntip(t)
  p0 <- ggtree(ladderize(t), size = 0.45); dd <- p0$data
  dd$sup <- NA_real_; ni <- !dd$isTip; dd$sup[ni] <- sup[dd$node[ni] - n]
  tp <- dd[dd$isTip, ]; tp$fill <- unname(colmap[tp$label]); tp$own <- tp$label %in% PPS; tp$label <- gsub("_", "-", tp$label)
  xm <- max(dd$x); off <- 0.012 * xm
  p <- p0 + geom_point(data = tp[!is.na(tp$fill), ], aes(x = x, y = y, fill = fill), shape = 21, size = 2.4, colour = "black") + scale_fill_identity() +
    geom_text(data = tp[tp$own, ], aes(x = x + off, y = y, label = label), hjust = 0, size = 2.7, fontface = "bold", colour = "black") +
    geom_text(data = tp[!tp$own, ], aes(x = x + off, y = y, label = label), hjust = 0, size = 2.7, colour = "black") +
    geom_text(data = dd[ni & !is.na(dd$sup) & dd$sup >= 70, ], aes(x = x, y = y, label = round(sup)), hjust = 1.15, vjust = -0.5, size = 2.2, colour = "black") +
    geom_treescale(x = 0, y = -1, fontsize = 2.6, linesize = 0.4) + scale_x_continuous(limits = c(NA, xm * 1.5), expand = c(0, 0))
  if (bl) p <- p + geom_text(data = dd[!is.na(dd$branch.length) & ni, ], aes(x = branch, y = y, label = sprintf("%.3f", branch.length)), size = 1.9, vjust = 1.4, colour = "black")
  h <- max(6, 0.15 * n + 1)
  ggsave(paste0(out, ".png"), p, width = 9, height = h, dpi = 600, limitsize = FALSE, bg = "white")
  ggsave(paste0(out, ".pdf"), p, width = 9, height = h, limitsize = FALSE, bg = "white"); cat("done", out, "\n")
}
draw("09_core_tree/broad84_reused/core_tree/core84_tree.contree", "12_marker_gene_trees/figures_pub/tree_core84_panaroo")
