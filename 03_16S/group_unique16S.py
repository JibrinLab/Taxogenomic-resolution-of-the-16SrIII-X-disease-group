#!/usr/bin/env python3
# ANALYSIS: 16S rRNA gene analysis
# MANUSCRIPT SECTION: Multilocus sequence comparison and 16S classification
# USED FOR: Figs. S2.1 to S2.7, Tables S2.1 to S2.6
# FOLDER: 03_16S (see README.md in the parent folder)
"""Unique 16S rRNA gene regions of Groups 1-4 (our 22 genomes).
usage: python3 group_unique16S.py OUTDIR 16S_FASTA [--bg NCBI_FASTA] [--k 10] [--frac 1.0] [--same G2=pruni]
16S_FASTA : 16S sequences of the 22 genomes (aligned or not; gaps are removed), one or more copies per genome
--bg      : NCBI 16S sequences of described Ca. Phytoplasma species (e.g. seqs/_all_candidatus.fasta)
--same    : background species that belongs to a group (removed from that group's background), e.g. G2=pruni"""
import sys, os, re, csv, argparse
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import phyto_unique16S as P

GROUPS = {"Group 1": ["PPS1A","PPS2B","PPS6A","PPS6B","PPS9B","PPS10A","CanS4","MW1","VAC","MA1"],
          "Group 2": ["CX","6A1","ZT3-1","PR2021","JR1","PDA15"],
          "Group 3": ["CicWB-2022","Vc33","ChTDIII","OVSP-L01","OVSP-S01"],
          "Group 4": ["LWB"]}
key = lambda s: re.sub(r"[^A-Za-z0-9]", "", s).upper()
LK = {key(n): g for g, m in GROUPS.items() for n in m}
def group_of(h):
    k = key(h.split()[0])
    if k in LK: return LK[k]
    for kk in sorted(LK, key=len, reverse=True):
        if kk in key(h): return LK[kk]
    return None

def best_window(seq, pool, L=20):
    if len(seq) <= L: return None
    best = None
    for i in range(len(seq) - L + 1):
        mm, who = P.min_mismatch(seq[i:i + L], pool)
        if best is None or mm > best[0]: best = (mm, i, who)
    return best

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("out"); ap.add_argument("fasta")
    ap.add_argument("--bg"); ap.add_argument("--k", type=int, default=10); ap.add_argument("--frac", type=float, default=1.0)
    ap.add_argument("--cap", type=int, default=30); ap.add_argument("--same", action="append", default=[])
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    same = dict(x.split("=") for x in a.same)                       # e.g. {"G2": "pruni"}
    G = defaultdict(list)
    for h, s in P.read_fasta(open(a.fasta).read()):
        s = s.replace("-", "").replace(".", ""); g = group_of(h)
        if g: G[g].append((h, s))
        else: print("not assigned to a group:", h)
    bg = defaultdict(list)
    if a.bg:
        for h, s in P.read_fasta(open(a.bg).read()):
            sp = P.species_of(h)
            if sp and len(s) >= 1400 and sum(c not in "ACGT" for c in s) <= 5 and len(bg[sp]) < a.cap: bg[sp].append((h, s))
    with open(os.path.join(a.out, "group_16S_sequences.tsv"), "w") as f:
        f.write("group\tgenome\tlength\n")
        for g in sorted(G):
            for h, s in G[g]: f.write("%s\t%s\tlen %d\n" % (g, h.split()[0], len(s)))
    refs = {}; rows = []
    for g in sorted(G):
        S = {g: G[g]}
        pool_groups = {o: G[o] for o in G if o != g}
        pool_all = dict(pool_groups)
        for sp, lst in bg.items():
            if same.get(g.replace("roup ", "")) == sp: continue
            pool_all["ncbi_" + sp] = lst
        print("%s: %d sequences, lengths %s" % (g, len(G[g]), sorted({len(s) for h, s in G[g]})))
        for tag, pool in (("vs_other_groups", pool_groups), ("vs_other_groups_and_described_species", pool_all)):
            ref, regs = P.regions_for(g, S, a.k, a.frac, pool)
            if ref is None: continue
            refs[g] = ref
            cand = []
            for (st, en, L, seq, ns, n) in regs:
                mm, who = P.min_mismatch(seq, pool)
                cand.append(dict(group=g, comparison=tag, reference=ref[0].split()[0], start=st, end=en, length=L, sequence=seq,
                                 strains_with_region="%d of %d" % (ns, n), min_mismatches_to_others=mm, closest_other=who, core20=""))
            cand.sort(key=lambda r: (-r["min_mismatches_to_others"], -r["length"]))
            for r in cand[:12]:                                     # best 20-nt window inside long regions
                if r["length"] > 25 and r["min_mismatches_to_others"] >= 2:
                    bw = best_window(r["sequence"], pool)
                    if bw: r["core20"] = "%d-%d %s (%d mm)" % (r["start"] + bw[1], r["start"] + bw[1] + 19, r["sequence"][bw[1]:bw[1] + 20], bw[0])
            rows += cand
            print("  %s: %d regions, %d with >=2 mismatches" % (tag, len(cand), sum(r["min_mismatches_to_others"] >= 2 for r in cand)))
    cols = ["group", "comparison", "reference", "start", "end", "length", "sequence", "strains_with_region", "min_mismatches_to_others", "closest_other", "core20"]
    with open(os.path.join(a.out, "group_unique_16S_regions.tsv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(rows)
    with open(os.path.join(a.out, "group_reference_16S.fasta"), "w") as f:
        for g, (h, s) in refs.items(): f.write(">%s | %s\n%s\n" % (g, h, s))
    print("done ->", a.out)
if __name__ == "__main__": main()
