# ANALYSIS: Potential mobile unit (PMU) genes in the 22 genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Fig. 7, Tables S6.1 to S6.5
# FOLDER: 08_PMU (see README.md in the parent folder)
suppressMessages(library(ggplot2))
w <- read.delim("results/pmu_window_genes.tsv"); w$Genome <- gsub("OVSP_", "OVSP-", w$Genome)
GN <- c("smc", "ftsH", "hupB", "ssb", "sigF")
cols <- c(smc = "#56B4E9", ftsH = "#E69F00", hupB = "#009E73", ssb = "#CC79A7", sigF = "#D55E00")
info <- do.call(rbind, lapply(unique(w$Genome), function(g) { d <- w[w$Genome == g, ]
  data.frame(Genome = g, x0 = min(d$Start), W = max(d$End) - min(d$Start)) }))
info$panel <- as.integer(cut(info$W, c(0, 7000, 20000, Inf)))
w <- merge(w, info, by = "Genome")
w$xs <- w$Start - w$x0; w$xe <- w$End - w$x0
w$dir <- ifelse(w$Strand == "+", 1, -1); w$Gene <- factor(w$Gene, levels = GN)
draw <- function(k, first, last) {
  inf <- info[info$panel == k, ]; inf$yi <- (nrow(inf):1) * 2.0
  d <- merge(w[w$panel == k, ], inf[, c("Genome", "yi")], by = "Genome")
  span <- c(7000, 20000, 38000)[k]
  poly <- do.call(rbind, lapply(split(d, seq_len(nrow(d))), function(r) {
    len <- r$xe - r$xs; hl <- min(0.5 * len, 0.022 * span)
    u <- c(0, len - hl, len - hl, len, len - hl, len - hl, 0); v <- c(-.2, -.2, -.33, 0, .33, .2, .2)
    tl <- if (r$dir == 1) r$xs else r$xe
    data.frame(x = tl + r$dir * u, y = r$yi + v, id = paste(r$Genome, r$Gene), Gene = r$Gene) }))
  d <- d[order(d$Genome, d$xs), ]
  d$ly <- d$yi + ave(d$xs, d$Genome, FUN = function(z) rep(c(0.5, 0.85, 1.2), length.out = length(z)))
  ggplot() + geom_segment(data = inf, aes(x = 0, xend = W, y = yi, yend = yi), linewidth = 0.3, colour = "black") +
    geom_polygon(data = poly, aes(x, y, group = id, fill = Gene), colour = "black", linewidth = 0.25) +
    geom_text(data = d, aes(x = (xs + xe) / 2, y = ly, label = Gene), size = 2.8, fontface = "italic", colour = "black", vjust = 0) +
    scale_fill_manual(values = cols, name = NULL, drop = FALSE) +
    scale_y_continuous(breaks = inf$yi, labels = inf$Genome, expand = expansion(add = c(0.6, 1.9))) +
    scale_x_continuous(labels = function(v) v / 1000, limits = c(0, NA), expand = expansion(mult = c(0.035, 0.03))) +
    labs(x = if (last) "Position from the start of the window (kb); contig coordinates as annotated" else NULL, y = NULL) +
    theme_classic(base_size = 9) + theme(text = element_text(colour = "black"), axis.text = element_text(colour = "black"),
      axis.text.y = element_text(face = "bold"), legend.position = if (first) "top" else "none", legend.key.size = unit(0.35, "cm")) + coord_cartesian(clip = "off")
}
pl <- lapply(1:3, function(k) draw(k, k == 1, k == 3)); nr <- as.numeric(table(factor(info$panel, 1:3))); hts <- nr * 1.25 + c(1.8, 1.2, 1.6)
comb <- if (requireNamespace("patchwork", quietly = TRUE)) patchwork::wrap_plots(pl, ncol = 1, heights = hts) else if (requireNamespace("cowplot", quietly = TRUE)) cowplot::plot_grid(plotlist = pl, ncol = 1, align = "v", axis = "lr", rel_heights = hts) else stop("neither patchwork nor cowplot is installed")
ggsave("figures/pmu_22genomes.png", comb, width = 7.2, height = 8, dpi = 600)
ggsave("figures/pmu_22genomes.pdf", comb, width = 7.2, height = 8)
ggsave("figures/pmu_22genomes_preview.png", comb, width = 7.2, height = 8, dpi = 80)
