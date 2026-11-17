# ANALYSIS: Potential mobile unit (PMU) genes in the 22 genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Fig. 7, Tables S6.1 to S6.5
# FOLDER: 08_PMU (see README.md in the parent folder)
import re, csv, itertools, os
W = "/projects/dsn011/ppss_publication/2026-10-27_22genome_reanalysis"
T = W + "/05_metabolism_22genomes/gff_tree"
OUT = W + "/04_PMU_22genomes/results"
os.makedirs(OUT, exist_ok=True)
G = {"Group 1": "PPS1A PPS2B PPS6A PPS6B PPS9B PPS10A CanS4 MW1 VAC MA1".split(),
     "Group 2": "CX 6A1 ZT3-1 PR2021 JR1 PDA15".split(),
     "Group 3": "CicWB-2022 Vc33 ChTDIII OVSP-L01 OVSP-S01".split(),
     "Group 4": ["LWB"]}
labs = [l for v in G.values() for l in v]
grp = {l: n for n, m in G.items() for l in m}
AR = re.compile(r'(\w+)=([^;]*)')
GENES = ["sigF", "ssb", "hupB", "ftsH2", "smc"]
def b(x): return re.sub(r'_\d+$', '', x)
STRICT = {"sigF": lambda x: x == "sigF", "ssb": lambda x: x in ("ssb", "ssbA"),
          "hupB": lambda x: x == "hupB", "ftsH2": lambda x: x == "ftsH2", "smc": lambda x: x == "smc"}
RELAX = dict(STRICT); RELAX["ftsH2"] = lambda x: x.startswith("ftsH")
def parse(p):
    out = []
    for line in open(p):
        if line.startswith("##FASTA"): break
        if line.startswith("#"): continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 9 or c[2] != "CDS": continue
        a = dict(AR.findall(c[8]))
        out.append((c[0], int(c[3]), int(c[4]), c[6], (a.get("gene") or a.get("Name") or "").strip(), a.get("product", "")))
    return out
def hits(feats, rules):
    h = {g: [] for g in GENES}
    for f in feats:
        n = b(f[4])
        for g in GENES:
            if n and rules[g](n): h[g].append(f)
    return h
def window(h):
    best = None
    for c in set.intersection(*[set(f[0] for f in h[g]) for g in GENES]):
        L = [[f for f in h[g] if f[0] == c] for g in GENES]
        for combo in itertools.product(*L):
            s = min(f[1] for f in combo); e = max(f[2] for f in combo)
            if best is None or e - s < best[0]:
                o = sorted(zip(combo, GENES), key=lambda x: x[0][1])
                best = (e - s, c, "-".join(g for _, g in o), "".join(f[3] for f, _ in o))
    return best
T3 = lambda p: csv.writer(open(OUT + p, "w"), delimiter="\t")
lg = T3("/pmu_hits_long.tsv"); lg.writerow(["Genome", "Group", "Gene", "Contig", "Start", "End", "Strand", "GffName", "Product"])
cn = T3("/pmu_copy_numbers.tsv"); cn.writerow(["Genome", "Group"] + [g + "_copies" for g in GENES] + [g + "_contigs" for g in GENES])
cs = T3("/pmu_cluster_summary.tsv"); cs.writerow(["Genome", "Group", "Rule", "AllFiveOnOneContig", "MinSpan_bp", "Contig", "Order", "Strands"])
for l in labs:
    p = T + "/" + l + "/" + l + ".gff"
    if not os.path.exists(p): print("WARNING missing GFF", l); continue
    feats = parse(p)
    h = hits(feats, STRICT)
    for g in GENES:
        for f in h[g]: lg.writerow([l, grp[l], g] + list(f))
    cn.writerow([l, grp[l]] + [len(h[g]) for g in GENES] + [len(set(f[0] for f in h[g])) for g in GENES])
    for name, rules in (("strict-ftsH2-name", STRICT), ("any-ftsH", RELAX)):
        hh = hits(feats, rules)
        w = window(hh) if all(hh[g] for g in GENES) else None
        cs.writerow([l, grp[l], name, "yes" if w else "no"] + (list(w) if w else ["", "", "", ""]))
    print(l, [len(h[g]) for g in GENES], len(feats), "CDS")
