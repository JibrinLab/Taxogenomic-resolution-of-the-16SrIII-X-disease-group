# ANALYSIS: Core genome tree, pangenome and recombination
# MANUSCRIPT SECTION: Phylogenetic reconstructions; Comparative Genomics Analyses
# USED FOR: Fig. 5, Fig. 6, Figs. S5.1 to S5.6
# FOLDER: 07_core_tree_pangenome_recombination (see README.md in the parent folder)
suppressPackageStartupMessages({library(ape); library(ggplot2); library(ggtree)})
a <- commandArgs(TRUE); dir <- a[1]; pre <- file.path(dir, a[2]); out <- a[3]; ttl <- a[4]
tr <- read.tree(paste0(pre, ".labelled_tree.newick"))
imp <- read.table(paste0(pre, ".importation_status.txt"), header = TRUE, stringsAsFactors = FALSE)
xr <- paste0(pre, ".position_cross_reference.txt")
L <- max(imp$End, if (file.exists(xr)) length(strsplit(readLines(xr, n = 1), ",")[[1]]) else 0)
n <- Ntip(tr); p <- ggtree(tr, size = 0.5, color = "black"); d <- p$data; mx <- max(d$x)
TW <- 3.2; LW <- 1.4; PW <- 7; s <- TW / mx; x0 <- mx + LW / s; pw <- PW / s
br <- d[d$parent != d$node, ]
rows <- data.frame(y = br$y, xmin = x0, xmax = x0 + pw)
imp$y <- d$y[match(imp$Node, d$label)]; imp <- imp[!is.na(imp$y), ]
fix <- c(CicWB_2022 = "CicWB-2022", ZT3_1 = "ZT3-1", OVSP_L01 = "OVSP-L01", OVSP_S01 = "OVSP-S01", A_laidlawii = "A. laidlawii")
lab <- p$data$label; hit <- p$data$isTip & lab %in% names(fix); p$data$label[hit] <- fix[lab[hit]]
imp$xmin <- x0 + imp$Beg / L * pw; imp$xmax <- x0 + imp$End / L * pw
tk <- pretty(c(0, L), n = 6); tk <- tk[tk <= L]
ax <- data.frame(x = x0 + tk / L * pw, lab = format(tk, big.mark = ",", scientific = FALSE, trim = TRUE))
p2 <- p + geom_tiplab(align = TRUE, linetype = "dotted", linesize = 0.3, size = 4, color = "black", offset = 0.01 * mx) +
  geom_rect(data = rows, aes(xmin = xmin, xmax = xmax, ymin = y - 0.38, ymax = y + 0.38), fill = "#87CEEB", inherit.aes = FALSE) +
  geom_rect(data = imp, aes(xmin = xmin, xmax = pmax(xmax, xmin + 0.002 * pw), ymin = y - 0.38, ymax = y + 0.38), fill = "#B22222", inherit.aes = FALSE) +
  geom_segment(data = ax, aes(x = x, xend = x, y = 0.3, yend = 0.45), inherit.aes = FALSE) +
  geom_text(data = ax, aes(x = x, y = -0.1, label = lab), size = 3.5, color = "black", inherit.aes = FALSE) +
  annotate("text", x = x0 + pw / 2, y = -0.9, label = "Position in alignment (bp)", size = 4.2, color = "black") +
  geom_treescale(x = 0, y = 0.1, fontsize = 4, linesize = 0.6, color = "black") +
  labs(title = ttl, caption = "Blue: clonal frame. Red: inferred recombination (ClonalFrameML). Scale bar: substitutions per site.") +
  coord_cartesian(xlim = c(-0.02 * mx, x0 + pw), ylim = c(-1.5, n + 1), expand = FALSE, clip = "off") +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5, color = "black"), plot.caption = element_text(size = 10, color = "black", hjust = 0),
        plot.margin = margin(10, 25, 10, 10), plot.background = element_rect(fill = "white", color = NA))
w <- TW + LW + PW + 1.5; h <- max(5, 0.30 * 2 * n + 2.5)
ggsave(paste0(out, ".png"), p2, width = w, height = h, dpi = 600, bg = "white")
ggsave(paste0(out, ".pdf"), p2, width = w, height = h, bg = "white")
ggsave(paste0(out, "_preview.png"), p2, width = w, height = h, dpi = 80, bg = "white")
imp$bp <- imp$End - imp$Beg + 1; imp$n <- 1
sm <- aggregate(cbind(n, bp) ~ Node, imp, sum)
write.table(sm, paste0(out, "_branch_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(a[2], "events", nrow(imp), "total bp", sum(imp$bp), "alignment length used", L, "\n")
