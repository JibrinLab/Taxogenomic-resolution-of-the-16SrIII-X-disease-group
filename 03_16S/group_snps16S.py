#!/usr/bin/env python3
# ANALYSIS: 16S rRNA gene analysis
# MANUSCRIPT SECTION: Multilocus sequence comparison and 16S classification
# USED FOR: Figs. S2.1 to S2.7, Tables S2.1 to S2.6
# FOLDER: 03_16S (see README.md in the parent folder)
"""Group-diagnostic nucleotide positions in the 16S rRNA gene alignment of Groups 1-4.
usage: python3 group_snps16S.py OUTDIR ALIGNED_16S_FASTA [--bg NCBI_FASTA]
A position is diagnostic for a group when every sequence of the group has the same base (no gap) and no
sequence of any other group has that base. Each is checked for the same 21-nt window (10 nt either side)
in the NCBI sequences of described 'Candidatus Phytoplasma' species."""
import sys, os, argparse
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import phyto_unique16S as P
from group_unique16S import group_of

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("out"); ap.add_argument("fasta"); ap.add_argument("--bg"); ap.add_argument("--cap", type=int, default=30); ap.add_argument("--frac", type=float, default=1.0)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    recs = [(h, s.upper().replace("U", "T")) for h, s in P.read_fasta(open(a.fasta).read())]
    L = {len(s) for h, s in recs}
    if len(L) != 1: sys.exit("sequences are not aligned (different lengths): %s" % sorted(L)[:5])
    L = L.pop(); G = defaultdict(list)
    for h, s in recs:
        g = group_of(h)
        if g: G[g].append((h.split()[0], s))
        else: print("not assigned:", h)
    bg = defaultdict(list)
    if a.bg:
        for h, s in P.read_fasta(open(a.bg).read()):
            sp = P.species_of(h)
            if sp and len(bg[sp]) < a.cap: bg[sp].append((h.split()[0], s))
    # identity between groups (alignment columns without gaps, mean over sequence pairs)
    def ident(x, y):
        n = d = 0
        for i in range(L):
            if x[i] in "ACGT" and y[i] in "ACGT": n += 1; d += x[i] != y[i]
        return 100.0 * (1 - d / n) if n else float("nan")
    gs = sorted(G); lines = ["group_a\tgroup_b\tmean_identity_percent\tmin_identity_percent"]
    for i, g1 in enumerate(gs):
        for g2 in gs[i:]:
            v = [ident(x[1], y[1]) for x in G[g1] for y in G[g2] if x is not y]
            if v: lines.append("%s\t%s\t%.2f\t%.2f" % (g1, g2, sum(v) / len(v), min(v)))
    open(os.path.join(a.out, "16S_identity_between_groups.tsv"), "w").write("\n".join(lines) + "\n")
    ref0 = next((s for h, s in recs if h.startswith("PPS1A")), recs[0][1])
    rows = []
    for g in gs:
        seqs = [s for h, s in G[g]]; others = {o: [s for h, s in G[o]] for o in gs if o != g}
        gref = seqs[0]
        for c in range(L):
            col = [s[c] for s in seqs]
            base = max("ACGT", key=col.count)
            if col.count(base) < a.frac * len(col): continue
            ob = {o: sorted({s[c] for s in ss}) for o, ss in others.items()}
            if any(base in v for v in ob.values()): continue
            rp = sum(1 for x in gref[:c + 1] if x != "-"); pp = sum(1 for x in ref0[:c + 1] if x != "-")
            ungapped = gref.replace("-", ""); w0 = max(0, rp - 1 - 10); win = ungapped[w0:rp + 10]
            hit = sorted({sp for sp, lst in bg.items() for h, s in lst if win in s.replace("-", "")}) if bg else []
            rows.append((g, c + 1, rp, pp, "%s (%d of %d)" % (base, col.count(base), len(col)), ";".join("%s:%s" % (o, "/".join(v)) for o, v in ob.items()), win, ",".join(hit) if hit else "none"))
    with open(os.path.join(a.out, "group_diagnostic_16S_positions.tsv"), "w") as f:
        f.write("group\talignment_column\tposition_in_group_sequence\tposition_in_PPS1A\tbase_in_group\tbases_in_other_groups\twindow_21nt\twindow_found_in_described_species\n")
        for r in rows: f.write("\t".join(str(x) for x in r) + "\n")
    cnt = defaultdict(lambda: [0, 0])
    for r in rows: cnt[r[0]][0] += 1; cnt[r[0]][1] += r[7] == "none"
    for g in gs: print("%s: %d sequences, %d diagnostic positions, %d with a window absent from described species" % (g, len(G[g]), cnt[g][0], cnt[g][1]))
    print(open(os.path.join(a.out, "16S_identity_between_groups.tsv")).read())
if __name__ == "__main__": main()
