# ANALYSIS: Secretome (SignalP, DIAMOND) and effector search in the 22 genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Fig. S6.2, Tables S8.x, S9.x, S10.x
# FOLDER: 10_secretome_effectors (see README.md in the parent folder)
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
cn <- read.delim("results/Table_secretome_Counts.tsv", check.names = FALSE)
pg <- read.delim("results/Table_secretome_conservation_Per_genome.tsv", check.names = FALSE)
ord <- cn$Genome; cn$Genome <- factor(cn$Genome, levels = ord)
blk <- element_text(colour = "black")
th <- theme_classic(base_size = 9) + theme(text = blk, axis.text = blk, axis.title = blk, legend.text = blk, legend.position = "top", axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black"))
A <- ggplot(cn, aes(Genome, Secreted_SP_Sec_SPI)) + geom_col(fill = "#9ecae1", colour = "black", linewidth = .2, width = .75) +
  geom_text(aes(label = paste0(Secreted_pct_of_proteins, "%")), vjust = -0.4, size = 2.4, colour = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, .12))) + labs(x = NULL, y = "SP(Sec/SPI) proteins") + th +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
d <- pg %>% transmute(Genome = factor(Genome, levels = ord), `Hit in a reference` = With_hit_in_reference, `No reference hit` = No_hit_in_any_reference) %>%
  pivot_longer(-Genome, names_to = "Class", values_to = "n")
d$Class <- factor(d$Class, levels = c("No reference hit", "Hit in a reference"))
B <- ggplot(d, aes(Genome, n, fill = Class)) + geom_col(colour = "black", linewidth = .2, width = .75) +
  scale_fill_manual(values = c("Hit in a reference" = "#9ecae1", "No reference hit" = "#fdae6b"), breaks = c("Hit in a reference", "No reference hit"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) + labs(x = NULL, y = "Secreted candidates") + th +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
p <- if (requireNamespace("patchwork", quietly = TRUE)) patchwork::wrap_plots(A, B, ncol = 1) else cowplot::plot_grid(A, B, ncol = 1, align = "v", axis = "lr")
ggsave("figures/secretome_22genomes.png", p, width = 7.2, height = 6, dpi = 600, bg = "white")
ggsave("figures/secretome_22genomes.pdf", p, width = 7.2, height = 6, bg = "white")
ggsave("figures/secretome_22genomes_preview.png", p, width = 7.2, height = 6, dpi = 120, bg = "white")
