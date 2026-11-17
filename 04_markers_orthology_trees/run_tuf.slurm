#!/bin/bash
# ANALYSIS: Marker genes, orthology and phylogenetic trees
# MANUSCRIPT SECTION: Multilocus sequence comparison; Phylogenetic reconstructions
# USED FOR: Fig. S2.8, Figs. S4.1 to S4.6
# FOLDER: 04_markers_orthology_trees (see README.md in the parent folder)
#SBATCH -J tuf_tree
#SBATCH -p batch
#SBATCH -N 1 -n 8
#SBATCH --mem=16G
#SBATCH -t 4:00:00
#SBATCH -o /projects/dsn011/ppss_publication/2026-10-27_22genome_reanalysis/12_marker_gene_trees/tuf_%j.out
export PATH=/projects/dsn011/urbello_conda/envs/mafft_env/bin:/opt/blast+/2.16.0/prebuilt/bin:$PATH
cd /projects/dsn011/ppss_publication/2026-10-27_22genome_reanalysis/12_marker_gene_trees
G22=/projects/dsn011/ppss_publication/2026-10-27_22genome_reanalysis/06_secretome_effector_22genomes/genomes22
cp /projects/dsn011/ppss_publication/2026-10-20_new_species_taxonomy/14_pbt_genbank_refs/curated84_contigs.fna combined.fna
for f in $G22/*.fna; do l=$(basename $f .fna); grep -qx "${l//-/_}" curated84_labels.txt || grep -qx "$l" curated84_labels.txt || { echo "added $l"; sed "s/^>/>${l}|/" $f >> combined.fna; }; done
OUT=$(grep -i -m1 laidlawii curated84_labels.txt); echo "outgroup $OUT"
makeblastdb -in combined.fna -dbtype nucl -out tuf_db > /dev/null
blastn -version | head -n 1; mafft --version 2>&1 | head -n 1; iqtree --version 2>&1 | head -n 1
blastn -task blastn -query tuf_query_PPS1A.fna -db tuf_db -evalue 1e-20 -max_target_seqs 20000 -num_threads 8 -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" -out tuf_blast.tsv
python3 extract_tuf.py combined.fna tuf_blast.tsv
for f in $G22/*.fna; do l=$(basename $f .fna); u=${l//-/_}; [ "$u" != "$l" ] && sed -i "s/^>${u}\$/>${l}/" tuf_sequences.fna; done
mafft --auto --thread 8 tuf_sequences.fna > tuf.aln.fasta
iqtree -s tuf.aln.fasta -m MFP -B 1000 -T 8 -o $OUT -redo --prefix tuf_tree
echo finished
