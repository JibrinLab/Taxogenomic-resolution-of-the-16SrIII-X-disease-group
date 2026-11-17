# ANALYSIS: Potential mobile unit (PMU) genes in the 22 genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Fig. 7, Tables S6.1 to S6.5
# FOLDER: 08_PMU (see README.md in the parent folder)
import csv, collections
rows = list(csv.DictReader(open("results/pmu_hits_long.tsv"), delimiter="\t"))
print("Gene values:", sorted({r["Gene"] for r in rows}))
def cls(g):
    g = g.lower()
    for k, p in (("sigF","sigf"),("ssb","ssb"),("hupB","hup"),("ftsH","ftsh"),("smc","smc")):
        if g.startswith(p): return k
by = collections.defaultdict(set)
for r in rows:
    c = cls(r["Gene"])
    if c: by[(r["Genome"], r["Contig"])].add((int(r["Start"]), int(r["End"]), r["Strand"], c, r["GffName"]))
win = {}
for (g, ct), S in by.items():
    L = sorted(S)
    for i in range(len(L)):
        for j in range(i, len(L)):
            sub = L[i:j+1]
            if {x[3] for x in sub} == {"sigF","ssb","hupB","ftsH","smc"}:
                span = max(x[1] for x in sub) - sub[0][0] + 1
                if g not in win or span < win[g][0]: win[g] = (span, ct, sub)
                break
out = open("results/pmu_windows.tsv", "w")
out.write("Genome\tContig\tSpan_bp\tGene\tStart\tEnd\tStrand\tGffName\tCopiesInWindow\n")
hyb = {r["Genome"]: r for r in csv.DictReader(open("results/pmu_hybrid_cluster.tsv"), delimiter="\t")}
for g in hyb:
    if g in win:
        span, ct, sub = win[g]; first = {}
        for x in sub: first.setdefault(x[3], x)
        for k, x in first.items():
            out.write(f"{g}\t{ct}\t{span}\t{k}\t{x[0]}\t{x[1]}\t{x[2]}\t{x[4]}\t{sum(1 for y in sub if y[3]==k)}\n")
    print(g, "| mine:", win[g][0] if g in win else "none", "| hybrid:", hyb[g]["AllFiveOnOneContig"], hyb[g]["MinSpan_bp"])
