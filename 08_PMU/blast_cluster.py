# ANALYSIS: Potential mobile unit (PMU) genes in the 22 genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Fig. 7, Tables S6.1 to S6.5
# FOLDER: 08_PMU (see README.md in the parent folder)
import glob, os, itertools, csv
Q = ["sigF", "ssb", "hupB", "ftsH2_1", "smc_1"]
out = csv.writer(open("results/pmu_blast_cluster.tsv", "w"), delimiter="\t")
out.writerow(["Genome", "AllFiveOnOneContig", "MinSpan_bp", "Contig", "Order", "Strands"])
for f in sorted(glob.glob("tblastn/*.tsv")):
    g = os.path.basename(f)[:-4]; H = {q: {} for q in Q}
    for l in open(f):
        c = l.rstrip("\n").split("\t")
        q = c[0].replace("PPS1A_", "")
        if q in H and float(c[5]) <= 1e-20 and float(c[2]) >= 60 and float(c[4]) >= 70:
            a, b = sorted((int(c[7]), int(c[8])))
            H[q][(c[1], a // 300)] = (c[1], a, b, "+" if int(c[9]) > 0 else "-")
    H = {q: list(v.values()) for q, v in H.items()}
    best = None
    if all(H[q] for q in Q):
        for ctg in set.intersection(*[set(h[0] for h in H[q]) for q in Q]):
            L = [[h for h in H[q] if h[0] == ctg] for q in Q]
            for combo in itertools.product(*L):
                s = min(h[1] for h in combo); e = max(h[2] for h in combo)
                if best is None or e - s < best[0]:
                    o = sorted(zip(combo, Q), key=lambda x: x[0][1])
                    best = (e - s, ctg, "-".join(q for _, q in o), "".join(h[3] for h, _ in o))
    out.writerow([g, "yes" if best else "no"] + (list(best) if best else ["", "", "", ""]))
