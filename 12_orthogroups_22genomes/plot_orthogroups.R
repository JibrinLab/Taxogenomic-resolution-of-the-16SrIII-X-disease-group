# ANALYSIS: OrthoFinder orthogroups of the 22 genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Figs. S6.3 to S6.5
# FOLDER: 12_orthogroups_22genomes (see README.md in the parent folder)
suppressMessages({library(ggplot2); library(dplyr)})
cats <- c("present in all genomes", "shared (PPS and reference)", "exclusive to own set")
cols <- c("#4C6A92", "#B9C7D9", "#E39B3D")
draw <- function(tsv, out, w) {
  d <- read.delim(tsv, check.names = FALSE, stringsAsFactors = FALSE)
  d <- d[d$Orthogroups > 0 | d$Category %in% cats[1:2], ]
  d$Category <- ifelse(d$Category %in% c("PPS only", "reference only"), cats[3], d$Category)
  d <- d %>% group_by(Genome, Category) %>% summarise(n = sum(Orthogroups), .groups = "drop")
  d$Genome <- factor(d$Genome, levels = unique(read.delim(tsv, check.names = FALSE)$Genome))
  d$Category <- factor(d$Category, levels = rev(cats))
  p <- ggplot(d, aes(Genome, n, fill = Category)) + geom_col(width = 0.8, colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = setNames(cols, cats), breaks = cats, name = NULL) +
    labs(x = NULL, y = "Orthogroups containing the genome") + theme_classic(base_size = 11) +
    theme(text = element_text(colour = "black"), axis.text = element_text(colour = "black"),
          axis.text.x = element_text(angle = 60, hjust = 1), axis.title = element_text(colour = "black"),
          legend.text = element_text(colour = "black"), legend.position = "top", axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"))
  ggsave(paste0("figures/", out, ".png"), p, width = w, height = 4.2, dpi = 600, bg = "white")
  ggsave(paste0("figures/", out, ".pdf"), p, width = w, height = 4.2, bg = "white")
}
draw("pg_22.tsv", "orthogroups_22genomes", 7.5)
draw("pg_4pps.tsv", "orthogroups_4PPS_16ref", 7.0)
draw("pg_cur84.tsv", "orthogroups_curated84_21genomes", 7.2)
