# ANALYSIS: Whole-genome synteny of the PPS strains against four reference genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Figs. S1.4 to S1.7
# FOLDER: 11_synteny (see README.md in the parent folder)
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
lens <- read.delim("genome_contig_lengths.tsv", header = FALSE, col.names = c("g", "contig", "len"), stringsAsFactors = FALSE)
cn <- c("S1","E1","S2","E2","L1","L2","IDY","LR","LQ","CR","CQ","TR","TQ")
prep <- function(rep, q) {
  rl <- lens[lens$g == rep, ]; rl$off <- cumsum(c(0, head(rl$len, -1))); ql <- lens[lens$g == q, ]
  co <- read.delim(sprintf("coords/%s_vs_%s.coords", rep, q), header = FALSE, col.names = cn, stringsAsFactors = FALSE)
  co$rs <- co$S1 + rl$off[match(co$TR, rl$contig)] - 1; co$re <- co$E1 + rl$off[match(co$TR, rl$contig)] - 1; co$rv <- co$S2 > co$E2
  a <- co %>% group_by(TQ) %>% summarise(score = weighted.mean((rs + re) / 2, L1), flip = sum(L2[rv]) > sum(L2[!rv]))
  ql$score <- a$score[match(ql$contig, a$TQ)]; ql$flip <- a$flip[match(ql$contig, a$TQ)]; ql$flip[is.na(ql$flip)] <- FALSE
  ql <- ql[order(is.na(ql$score), ql$score, -ql$len), ]; ql$qoff <- cumsum(c(0, head(ql$len, -1)))
  m <- match(co$TQ, ql$contig); fl <- ql$flip[m]; ln <- ql$len[m]
  tr <- function(p) ql$qoff[m] + ifelse(fl, ln - p + 1, p)
  data.frame(Strain = q, x1 = co$rs, x2 = co$re, y1 = tr(co$S2), y2 = tr(co$E2), rev = ifelse(xor(co$rv, fl), "Reverse", "Forward"), stringsAsFactors = FALSE)
}
blk <- element_text(colour = "black")
plot_rep <- function(rep, qs, tag, nc, w, h) {
  d <- bind_rows(lapply(qs, function(q) prep(rep, q))); d$Strain <- factor(d$Strain, levels = qs); rlen <- sum(lens$len[lens$g == rep])
  p <- ggplot(d) + geom_segment(aes(x = x1 / 1e6, xend = x2 / 1e6, y = y1 / 1e6, yend = y2 / 1e6, colour = rev), linewidth = .35, lineend = "round") +
    scale_colour_manual(values = c(Forward = "#2166ac", Reverse = "#b2182b"), name = NULL) +
    facet_wrap(~Strain, ncol = nc, scales = "free_y") + coord_cartesian(xlim = c(0, rlen / 1e6)) +
    labs(x = paste0(rep, " (Mb)"), y = "Genome (Mb)") + theme_bw(base_size = 9) +
    theme(text = blk, axis.text = blk, axis.title = blk, strip.text = blk, legend.text = blk, strip.background = element_blank(), panel.grid = element_blank(), panel.border = element_rect(colour = "black"), legend.position = "top")
  b <- sprintf("figures/synteny_%s_vs_%s", rep, tag)
  ggsave(paste0(b, ".png"), p, width = w, height = h, dpi = 600, bg = "white"); ggsave(paste0(b, ".pdf"), p, width = w, height = h, bg = "white")
  ggsave(paste0(b, "_preview.png"), p, width = w, height = h, dpi = 110, bg = "white")
}
P6 <- c("PPS1A", "PPS2B", "PPS6A", "PPS6B", "PPS9B", "PPS10A")
for (r in c("CanS4", "CX", "CicWB-2022", "LWB")) { plot_rep(r, P6, "PPS6", 3, 7.2, 5); plot_rep(r, P6[1:4], "PPS4", 2, 5.2, 5.2) }
