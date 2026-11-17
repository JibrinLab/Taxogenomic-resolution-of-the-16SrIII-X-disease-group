# ANALYSIS: Marker genes, orthology and phylogenetic trees
# MANUSCRIPT SECTION: Multilocus sequence comparison; Phylogenetic reconstructions
# USED FOR: Fig. S2.8, Figs. S4.1 to S4.6
# FOLDER: 04_markers_orthology_trees (see README.md in the parent folder)
import sys, collections
rc = lambda s: s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]
fa, cur = {}, None
for l in open(sys.argv[1]):
    l = l.rstrip()
    if l.startswith(">"): cur = l[1:].split()[0]; fa[cur] = []
    else: fa[cur].append(l)
fa = {k: "".join(v) for k, v in fa.items()}
H = collections.defaultdict(list); qlen = 1
for l in open(sys.argv[2]):
    t = l.split("\t"); qlen = int(t[4]); a, b = int(t[8]), int(t[9])
    H[(t[1], "+" if a < b else "-")].append((min(a, b), max(a, b), int(t[3]), float(t[2]), int(t[6]), int(t[7])))
L = collections.defaultdict(list)
def flush(c, k):
    lo, hi = min(x[0] for x in c), max(x[1] for x in c); aln = sum(x[2] for x in c); qc = (max(x[5] for x in c) - min(x[4] for x in c) + 1) / qlen
    L[k[0].split("|")[0]].append((aln, k[0], lo, hi, k[1], sum(x[3] * x[2] for x in c) / aln, qc))
for k, hs in H.items():
    hs.sort(); c = [hs[0]]
    for h in hs[1:]:
        if h[0] - max(x[1] for x in c) <= 300: c.append(h)
        else: flush(c, k); c = [h]
    flush(c, k)
out, tab = open("tuf_sequences.fna", "w"), open("tuf_extraction_table.tsv", "w")
tab.write("Genome\tContig\tStart\tEnd\tStrand\tMean_identity_pct\tQuery_coverage_pct\tRegion_bp\tOther_tuf_like_loci\n")
for g in sorted(L):
    b = sorted(L[g], reverse=True)[0]; s = fa[b[1]][b[2] - 1:b[3]]; s = rc(s) if b[4] == "-" else s
    out.write(f">{g}\n{s}\n"); tab.write("\t".join(str(x) for x in (g, b[1], b[2], b[3], b[4], f"{b[5]:.1f}", f"{100 * b[6]:.0f}", len(s), len(L[g]) - 1)) + "\n")
print(len(L), "genomes with a tuf hit")
