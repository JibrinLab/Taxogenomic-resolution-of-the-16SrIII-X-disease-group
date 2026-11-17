# Taxogenomic-resolution-of-the-16SrIII-X-disease-group
# Scripts for the PPS taxonomy manuscript: Taxogenomic resolution of the 16SrIII (X-disease) group: elevation of ‘Candidatus Phytoplasma pruni’-related strains to species rank and proposal of three novel Candidatus species ‘Candidatus Phytoplasma caryae’, ‘Candidatus Phytoplasma conii’ and ‘Candidatus Phytoplasma lactucae’.

Folders are numbered in the order of the Methods. Scripts are copies of the versions that produced the results in the manuscript and supplementary files. File paths inside the scripts are the absolute paths on Pete where they were run.

| Folder | Methods section | Results produced |
|---|---|---|
| 00_genome_tables | Genome sequencing and additional genomes | Genome statistics table of the six PPS genomes (make_table1.py) and the master table of the 22 genomes (build_master_table.py) |
| 01_genome_quality | Additional genomes used | CheckM v1 (lineage_wf) of the 299 genomes and of A. laidlawii (Table S1.6) |
| 02_ANI | Average nucleotide identity | FastANI of the 83 curated and 299 genomes (make_curated_ani.py, build_curated_workbook.py, make_pruni_ani_table.py); ANI of the 22 genomes plus PG-8A with FastANI, ANIm and ANIb (22genomes_PG8A: run_ani23.slurm, check_ani23.py, ani_cluster_summary*.py); heatmaps and tables (make_ani_figures.py, make_ani_figures_pub.py); linkage sensitivity at 94, 95 and 96% (run_ani_tests.py) |
| 03_16S | Multilocus sequence comparison; 16S classification | 16S rRNA gene sequences, group-specific positions and identity tables (Tables S2.1 to S2.5, Figs. S2.1 to S2.7: evidence_16S.py, group_snps16S.py with group_unique16S.py and phyto_unique16S.py, which it imports); comparison with the 60 described species (Table S2.6: fetch16S.py, lib16s.py, step2.py, stepA.py to stepG.py); 16S sequence sets (build_16S_sets.py, extract_16S.py) |
| 04_markers_orthology_trees | Multilocus comparison; phylogenetic reconstructions | Marker gene extraction and identity (a_long.py, b_sum.py, add_pbt_markers.py, add_cx_rpoB.py, run_secY_masked.slurm, fig_diff.py, prep_gene_aln.py; Fig. S2.8); OrthoFinder and Panaroo on the 84 genomes, the 92-gene tree and the core tree (find_markers.py, make_*.py, markers_pipeline.py, build_concat.py, run_orthofinder.slurm, run_concat_tree.slurm, run_panaroo84.slurm, run_core_tree.slurm); tuf trees (extract_tuf.py, run_tuf.slurm, plot_tuf.R); tree drawing (draw_trees.R, draw_core84.R, draw_fn.R) |
| 05_AAI_annotation | Phylogenetic reconstructions | Reciprocal best BLASTP identity between PPS1A and the other groups (rbh_aai.py); Prokka with genetic code 11 (run_prokka_t11.slurm) |
| 06_GTDB-Tk | ANI, dDDH and GTDB clustering | GTDB-Tk v2.6.1 on 53 genomes (Tables S5.1 to S5.3) |
| 07_core_tree_pangenome_recombination | Phylogenetic reconstructions; comparative genomics | Panaroo and core genome tree of the 22 genomes (Fig. 5); pangenome, Heaps' law and ClonalFrameML for Groups 1, 2 and 3 (Fig. 6, Figs. S5.1 to S5.6) |
| 08_PMU | Comparative genomics | Potential mobile unit genes (Fig. 7, Tables S6.1 to S6.5) |
| 09_metabolism | Comparative genomics | Metabolic gene panel, Prokka and BLASTp calls (Fig. 8, Tables S7.x) |
| 10_secretome_effectors | Comparative genomics | SignalP, DIAMOND comparison and effector search (Tables S8.x, S9.x, S10.x, Fig. S6.2) |
| 11_synteny | Comparative genomics | Whole-genome synteny of the PPS strains against four reference genomes (Figs. S1.4 to S1.7) |
| 12_orthogroups_22genomes | Comparative genomics | OrthoFinder on the 22 genomes (Figs. S6.3 to S6.5) |
| 13_trees_taxa | Figures | Trees with the four species marked (tree_taxa.py, tree_styles.py, build_taxon_map.py, draw_taxa_trees_v3.sh) |

## Steps without a script in this folder
- Read trimming, host-read removal and assembly (fastp, Trim Galore, Bowtie2, SPAdes) and assembly statistics.
- dDDH: computed on the GGDC web server (formulas 1, 2 and 3, BLAST+ option).
- PCR, qPCR and Sanger sequencing: laboratory work.
- Tree visualisation in iTOL.

