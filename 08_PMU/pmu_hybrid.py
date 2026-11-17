# ANALYSIS: Potential mobile unit (PMU) genes in the 22 genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Fig. 7, Tables S6.1 to S6.5
# FOLDER: 08_PMU (see README.md in the parent folder)
import re, csv, itertools, os
W = "/projects/dsn011/ppss_publication/2026-10-27_22genome_reanalysis"
T = W + "/05_metabolism_22genomes/gff_tree"
FA = "/projects/dsn011/ppss_publication/ani_23_with_PG8A/genomes/"
AR = re.compile(r'(\w+)=([^;]*)')
GENES = ["sigF", "ssb", "hupB", "ftsH", "smc"]
GN = "PPS1A PPS2B PPS6A PPS6B PPS9B PPS10A CanS4 MW1 VAC MA1 CX 6A1 ZT3-1 PR2021 JR1 PDA15 CicWB-2022 Vc33 ChTDIII OVSP_L01 OVSP_S01 LWB".split()
def b(x): return re.sub(r'_\d+$', '', x)
NAME = {"ssb": lambda x: x in ("ssb", "ssbA"), "ftsH": lambda x: x.startswith("ftsH"), "smc": lambda x: x == "smc",
        "sigF": lambda x: x == "sigF", "hupB": lambda x: x == "hupB"}
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
                a = dict(AR.findall(c[8])); cds.append((c[0], int(c[3]), int(c[4]), c[6], (a.get("gene") or a.get("Name") or "")))
    return reg, cds
def window(H):
    best = None
    if not all(H[g] for g in GENES): return None
    for ctg in set.intersection(*[set(h[0] for h in H[g]) for g in GENES]):
        L = [[h for h in H[g] if h[0] == ctg] for g in GENES]
        for combo in itertools.product(*L):
            s = min(h[1] for h in combo); e = max(h[2] for h in combo)
            if best is None or e - s < best[0]:
                o = sorted(zip(combo, GENES), key=lambda x: x[0][1])
                best = (e - s, ctg, "-".join(g for _, g in o), "".join(h[3] for h, _ in o), [(g, h[1], h[2], h[3]) for h, g in o])
    return best
out = csv.writer(open("results/pmu_hybrid_cluster.tsv", "w"), delimiter="\t")
out.writerow(["Genome", "AllFiveOnOneContig", "MinSpan_bp", "Contig", "Order", "Strands"] + [g + "_loci" for g in GENES])
WG = csv.writer(open("results/pmu_window_genes.tsv", "w"), delimiter="\t")
WG.writerow(["Genome", "Contig", "Span_bp", "Gene", "Start", "End", "Strand"])
for g in GN:
    lab = g.replace("_", "-"); p = T + "/" + lab + "/" + lab + ".gff"
    if not os.path.exists(p): print("missing GFF", g); continue
    fl = fa_len(FA + g + ".fna"); reg, cds = gff(p)
    def mp(n):
        if n in fl: return n
        m = [k for k, v in fl.items() if v == reg.get(n)]
        return m[0] if len(m) == 1 else None
    H = {k: {} for k in GENES}
    for ctg, s, e, st, nm in cds:
        for k in GENES:
            if nm and NAME[k](b(nm)) and mp(ctg): H[k][(mp(ctg), s // 300)] = (mp(ctg), s, e, st)
    for l in open("tblastn/" + g + ".tsv"):
        c = l.rstrip("\n").split("\t"); q = c[0].replace("PPS1A_", "")
        if q in ("sigF", "hupB") and float(c[5]) <= 1e-20 and float(c[2]) >= 60 and float(c[4]) >= 80:
            a, z = sorted((int(c[7]), int(c[8]))); H[q][(c[1], a // 300)] = (c[1], a, z, "+" if int(c[9]) > 0 else "-")
    H = {k: list(v.values()) for k, v in H.items()}
    w = window(H)
    out.writerow([g, "yes" if w else "no"] + (list(w[:4]) if w else ["", "", "", ""]) + [len(H[k]) for k in GENES])
    if w: [WG.writerow([g, w[1], w[0], n, a, z, st]) for n, a, z, st in w[4]]
