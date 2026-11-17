# ANALYSIS: Metabolic gene panel (66 genes, 7 pathways) in the 22 genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Fig. 8, Tables S7.1 to S7.10
# FOLDER: 09_metabolism (see README.md in the parent folder)
import sys, os, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/scripts")
import metab_lib as M
W = "/projects/dsn011/ppss_publication/2026-10-27_22genome_reanalysis/05_metabolism_22genomes"
GT = W + "/gff_tree"
OLD = "/projects/dsn011/ppss_publication/2026-10-16_metabolism_20genomes/tables/metabolic_all20_matrix.tsv"
GROUPS = {"Group 1": "PPS1A PPS2B PPS6A PPS6B PPS9B PPS10A CanS4 MW1 VAC MA1",
          "Group 2": "CX 6A1 ZT3-1 PR2021 JR1 PDA15",
          "Group 3": "CicWB-2022 ChTDIII OVSP-L01 OVSP-S01 Vc33", "Group 4": "LWB"}
norm = lambda s: s.replace("_", "-").replace("-Ca.P.pruni-16SrIII", "")
ordr = [x for g in GROUPS.values() for x in g.split()]
labs = sorted(os.listdir(GT), key=lambda l: ordr.index(norm(l)) if norm(l) in ordr else 99)
print("labels:", [norm(l) for l in labs], "n =", len(labs))
print("not in groups:", [l for l in labs if norm(l) not in ordr])
mat = M.run(labs, GT, "all22")
nl = {norm(l): l for l in labs}
O = M.OUT
w = csv.writer(open(O + "/metabolic_pps4_matrix.tsv", "w"), delimiter="\t")
p4 = ["PPS1A", "PPS2B", "PPS6A", "PPS6B"]
w.writerow(["Pathway", "Gene"] + p4)
for g in M.order: w.writerow([M.G[g][0], g] + [mat[g][nl[x]] for x in p4])
w = csv.writer(open(O + "/gene_by_group.tsv", "w"), delimiter="\t")
w.writerow(["Pathway", "Gene"] + ["%s (n=%d)" % (k, len(v.split())) for k, v in GROUPS.items()])
for g in M.order:
    w.writerow([M.G[g][0], g] + ["%d" % sum(1 for x in v.split() if mat[g][nl[x]] > 0) for v in GROUPS.values()])
w = csv.writer(open(O + "/genes_detected_per_genome.tsv", "w"), delimiter="\t")
w.writerow(["Genome", "Group", "GenesDetected", "of"])
for x in ordr:
    gr = [k for k, v in GROUPS.items() if x in v.split()][0]
    w.writerow([x, gr, sum(1 for g in M.order if mat[g][nl[x]] > 0), len(M.order)])
try:
    rd = csv.reader(open(OLD), delimiter="\t"); hd = next(rd); n = 0
    for r in rd:
        for j in range(2, len(hd)):
            x = norm(hd[j])
            if x in nl and r[1] in mat and (float(r[j] or 0) > 0) != (mat[r[1]][nl[x]] > 0):
                print("DIFF", r[1], x, "old", r[j], "new", mat[r[1]][nl[x]]); n += 1
    print("presence differences vs old matrix:", n)
except Exception as e: print("old-matrix compare skipped:", e)
