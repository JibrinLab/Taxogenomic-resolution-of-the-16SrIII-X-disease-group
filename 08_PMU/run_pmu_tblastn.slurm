#!/bin/bash
# ANALYSIS: Potential mobile unit (PMU) genes in the 22 genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Fig. 7, Tables S6.1 to S6.5
# FOLDER: 08_PMU (see README.md in the parent folder)
#SBATCH -J pmu_blast
#SBATCH -p batch
#SBATCH -c 4
#SBATCH --mem=4G
#SBATCH -t 01:00:00
#SBATCH -o pmu_blast_%j.out
module load blast+/2.16.0
cd /projects/dsn011/ppss_publication/2026-10-27_22genome_reanalysis/04_PMU_22genomes
mkdir -p db tblastn
for g in PPS1A PPS2B PPS6A PPS6B PPS9B PPS10A CanS4 MW1 VAC MA1 CX 6A1 ZT3-1 PR2021 JR1 PDA15 CicWB-2022 Vc33 ChTDIII OVSP_L01 OVSP_S01 LWB; do
  f=/projects/dsn011/ppss_publication/ani_23_with_PG8A/genomes/$g.fna
  makeblastdb -in $f -dbtype nucl -out db/$g -logfile /dev/null
  tblastn -query pmu_queries.faa -db db/$g -db_gencode 4 -evalue 1e-3 -num_threads 4 -outfmt "6 qseqid sseqid pident length qcovs evalue bitscore sstart send sframe" > tblastn/$g.tsv
done
echo done
