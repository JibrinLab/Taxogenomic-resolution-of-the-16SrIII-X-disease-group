# ANALYSIS: Whole-genome synteny of the PPS strains against four reference genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Figs. S1.4 to S1.7
# FOLDER: 11_synteny (see README.md in the parent folder)
module load mummer/4.0.0rc1
G=../../06_secretome_effector_22genomes/genomes22
N=/projects/dsn011/ppss_publication/2026-10-20_new_species_taxonomy/16_ddh_upload_PPS1A
mkdir -p old_assembly_backup; cp -n synteny_all_pairs.tsv genome_contig_lengths.tsv old_assembly_backup/
for s in PPS6A PPS10A; do mv dnadiff/*_vs_$s.* old_assembly_backup/ 2>/dev/null; done
for r in CanS4 CX CicWB-2022 LWB; do for s in PPS6A PPS10A; do dnadiff -p dnadiff/${r}_vs_$s $G/$r.fna $N/$s.fna > dnadiff/${r}_vs_$s.log 2>&1 & done; done; wait
grep -v -P "^(PPS6A|PPS10A)\t" genome_contig_lengths.tsv > t.tsv
for s in PPS6A PPS10A; do awk -v g=$s '/^>/{if(n!="")print g"\t"n"\t"l; n=substr($1,2); l=0; next}{l+=length($0)}END{print g"\t"n"\t"l}' $N/$s.fna; done >> t.tsv; mv t.tsv genome_contig_lengths.tsv
ls dnadiff/*.report | wc -l
