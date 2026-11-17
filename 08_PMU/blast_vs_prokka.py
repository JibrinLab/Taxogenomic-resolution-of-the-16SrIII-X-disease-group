# ANALYSIS: Potential mobile unit (PMU) genes in the 22 genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Fig. 7, Tables S6.1 to S6.5
# FOLDER: 08_PMU (see README.md in the parent folder)
import csv, re, os
W = "/projects/dsn011/ppss_publication/2026-10-27_22genome_reanalysis"
T = W + "/05_metabolism_22genomes/gff_tree"
FA = "/projects/dsn011/ppss_publication/ani_23_with_PG8A/genomes/"
AR = re.compile(r'(\w+)=([^;]*)')
def fa_len(p):
    d = {}; k = None
    for l in open(p):
        if l[0] == ">": k = l[1:].split()[0]; d[k] = 0
        else: d[k] += len(l.strip())
    return d
def gff(p):
    reg = {}; cds = []
    for l in open(p):
        if l.startswith("##FASTA"): break
        if l.startswith("##sequence-region"): x = l.split(); reg[x[1]] = int(x[3])
        elif l[0] != "#":
            c = l.rstrip("\n").split("\t")
            if len(c) > 8 and c[2] == "CDS":
                a = dict(AR.findall(c[8])); cds.append((c[0], int(c[3]), int(c[4]), a.get("gene", ""), a.get("product", "")))
    return reg, cds
rows = [r.rstrip("\n").split("\t") for r in open("pmu_blast_best.tsv")]
out = csv.writer(open("results/pmu_blast_vs_prokka.tsv", "w"), delimiter="\t")
out.writerow(["Genome", "Query", "pident", "qcov", "evalue", "FastaContig", "Start", "End", "ProkkaGene", "ProkkaProduct"])
cache = {}
for g, q, pid, cov, ev, ctg, s, e in rows:
    lab = g.replace("_", "-"); p = T + "/" + lab + "/" + lab + ".gff"
    if not os.path.exists(p): continue
    if g not in cache: cache[g] = (fa_len(FA + g + ".fna"),) + gff(p)
    fl, reg, cds = cache[g]
    L = fl.get(ctg)
    names = [n for n, ln in reg.items() if n == ctg or (ln == L and sum(1 for v in reg.values() if v == L) == 1)]
    lo, hi = sorted((int(s), int(e)))
    hit = [x for x in cds if x[0] in names and min(hi, x[2]) - max(lo, x[1]) > 0.5 * (hi - lo)]
    if not names: out.writerow([g, q, pid, cov, ev, ctg, s, e, "NO_MATCHING_CONTIG", ""])
    elif not hit: out.writerow([g, q, pid, cov, ev, ctg, s, e, "NO_CDS_OVERLAP", ""])
    for x in hit: out.writerow([g, q, pid, cov, ev, ctg, s, e, x[3], x[4]])
