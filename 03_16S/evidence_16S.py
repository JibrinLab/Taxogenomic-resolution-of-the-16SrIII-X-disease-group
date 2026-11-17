#!/usr/bin/env python3
# ANALYSIS: 16S rRNA gene analysis
# MANUSCRIPT SECTION: Multilocus sequence comparison and 16S classification
# USED FOR: Figs. S2.1 to S2.7, Tables S2.1 to S2.6
# FOLDER: 03_16S (see README.md in the parent folder)
"""16S rRNA gene evidence package for the four groups (tables, figures, proof).
usage: python3 evidence_16S.py OUTDIR ALIGNED_16S_FASTA [--bg NCBI_FASTA] [--frac 0.9] [--cap 30]
Run with the conda base env (numpy, pandas, matplotlib, openpyxl). No network needed.
Outputs (OUTDIR): Table1..Table5 tsv, evidence_16S_tables.xlsx, Fig1..Fig4 png/pdf, group_reference_16S.fasta, README.md"""
import sys, os, argparse
from collections import defaultdict
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import phyto_unique16S as P
from group_unique16S import group_of

REFS = {"Group 1": "PPS1A", "Group 2": "CX", "Group 3": "CicWB2022", "Group 4": "LWB"}
OUTLIERS = ("ChTDIII", "OVSP_L01_c2")
GCOL = {"Group 1": "#1b9e77", "Group 2": "#d95f02", "Group 3": "#7570b3", "Group 4": "#e7298a"}
BCOL = {"A": "#7fc97f", "C": "#80b1d3", "G": "#fdc086", "T": "#fb8072", "-": "#e0e0e0", "N": "#ffffff"}
k = lambda s: "".join(ch for ch in s if ch.isalnum()).upper()

def save(fig, base):
    fig.savefig(base + ".png", dpi=300, bbox_inches="tight"); fig.savefig(base + ".pdf", bbox_inches="tight"); plt.close(fig)

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("out"); ap.add_argument("fasta")
    ap.add_argument("--bg"); ap.add_argument("--cap", type=int, default=30); ap.add_argument("--frac", type=float, default=0.9)
    a = ap.parse_args(); out = a.out; os.makedirs(out, exist_ok=True)
    recs = [(h.split()[0], s.upper().replace("U", "T")) for h, s in P.read_fasta(open(a.fasta).read())]
    Ls = {len(s) for h, s in recs}
    if len(Ls) != 1: sys.exit("sequences are not aligned: %s" % sorted(Ls)[:5])
    L = Ls.pop(); G = defaultdict(list); seqs = {}
    for h, s in recs:
        g = group_of(h)
        if g: G[g].append(h); seqs[h] = s
    gs = sorted(G); order = [h for g in gs for h in sorted(G[g])]
    ref = {g: next((h for h in G[g] if k(h).startswith(k(REFS[g]))), G[g][0]) for g in gs}
    pps = seqs[ref["Group 1"]]
    cum = lambda s: np.cumsum([c != "-" for c in s])
    cumref = {h: cum(seqs[h]) for h in order}
    # ---- Table 1: reference 16S per group
    t1 = []
    with open(os.path.join(out, "group_reference_16S.fasta"), "w") as f:
        for g in gs:
            h = ref[g]; u = seqs[h].replace("-", "")
            gc = 100.0 * (u.count("G") + u.count("C")) / len(u)
            t1.append((g, h, len(u), "%.1f" % gc, len(G[g]), ", ".join(sorted(G[g]))))
            f.write(">%s %s %d bp\n%s\n" % (g.replace(" ", ""), h, len(u), u))
    T1 = pd.DataFrame(t1, columns=["group", "reference_sequence", "length_bp", "GC_percent", "n_sequences", "sequences"])
    # ---- Table 2: identity between sequences (alignment columns with no gap in either) and between groups
    A = np.array([list(seqs[h]) for h in order]); ok = np.isin(A, list("ACGT"))
    n = len(order); ID = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(n):
            m = ok[i] & ok[j]
            if m.sum(): ID[i, j] = 100.0 * (A[i][m] == A[j][m]).mean()
    gi = {g: [i for i, h in enumerate(order) if group_of(h) == g] for g in gs}
    rows = []
    for x in gs:
        for y in gs:
            v = [ID[i, j] for i in gi[x] for j in gi[y] if i != j]
            rows.append((x, y, "%.2f" % np.mean(v) if v else "", "%.2f" % min(v) if v else "", "%.2f" % max(v) if v else ""))
    T2 = pd.DataFrame(rows, columns=["group_a", "group_b", "mean_identity_percent", "min_identity_percent", "max_identity_percent"])
    T2s = pd.DataFrame(ID, index=order, columns=order).round(2)
    # ---- diagnostic positions
    bg = defaultdict(list)
    if a.bg:
        for h, s in P.read_fasta(open(a.bg).read()):
            sp = P.species_of(h)
            if sp and len(bg[sp]) < a.cap: bg[sp].append(s.upper().replace("-", ""))
    def diag(frac, exclude=()):
        res = []
        for g in gs:
            mem = [h for h in G[g] if not any(k(h).startswith(k(e)) for e in exclude)]
            oth = {o: [h for h in G[o] if not any(k(h).startswith(k(e)) for e in exclude)] for o in gs if o != g}
            for c in range(L):
                col = [seqs[h][c] for h in mem]
                if not col: continue
                b = max("ACGT", key=col.count)
                if col.count(b) < frac * len(col): continue
                if any(b in {seqs[h][c] for h in hs} for hs in oth.values()): continue
                res.append((g, c))
        return res
    D = diag(a.frac); D_all = set(diag(1.0)); D_ex = set(diag(1.0, OUTLIERS))
    t3 = []
    for g, c in D:
        h0 = ref[g]; u = seqs[h0].replace("-", "")
        rp = int(cumref[h0][c]); pp = int(cumref[ref["Group 1"]][c])
        col = [seqs[h][c] for h in G[g]]; b = max("ACGT", key=col.count)
        win = u[max(0, rp - 11):rp + 10]
        hit = sorted({sp for sp, lst in bg.items() if any(win in s for s in lst)}) if bg else []
        oth = "; ".join("%s: %s" % (o, "/".join(sorted({seqs[h][c] for h in G[o]}))) for o in gs if o != g)
        t3.append((g, pp, rp, c + 1, b, "%d of %d" % (col.count(b), len(col)), "yes" if (g, c) in D_all else "no",
                   "yes" if (g, c) in D_ex else "no", oth, win, ", ".join(hit) if hit else ("none" if bg else "not tested")))
    T3 = pd.DataFrame(t3, columns=["group", "position_PPS1A", "position_group_reference", "alignment_column", "group_base", "sequences_with_base",
                                   "in_all_sequences", "still_diagnostic_without_outlier_sequences", "bases_in_other_groups", "reference_21nt_window",
                                   "window_found_in_described_species"])
    # per-sequence proof matrix
    cols = [(g, c) for g, c in D]
    T4 = pd.DataFrame({("%s|%s" % (g, int(cumref[ref["Group 1"]][c]))): [seqs[h][c] for h in order] for g, c in cols}, index=order) if cols else pd.DataFrame(index=order)
    T4.insert(0, "group", [group_of(h) for h in order])
    for name, df, idx in (("Table1_reference_16S", T1, False), ("Table2_identity_between_groups", T2, False), ("Table3_diagnostic_positions", T3, False),
                          ("Table4_base_per_sequence", T4, True), ("Table5_identity_per_sequence", T2s, True)):
        df.to_csv(os.path.join(out, name + ".tsv"), sep="\t", index=idx)
    with pd.ExcelWriter(os.path.join(out, "evidence_16S_tables.xlsx")) as w:
        for name, df, idx in (("1_references", T1, False), ("2_group_identity", T2, False), ("3_diagnostic_positions", T3, False),
                              ("4_base_per_sequence", T4, True), ("5_sequence_identity", T2s, True)):
            df.to_excel(w, sheet_name=name, index=idx)
    # ---- Fig 1: identity heatmap per sequence
    fig, ax = plt.subplots(figsize=(7.5, 6.8))
    vmin = np.nanmin(ID); im = ax.imshow(ID, cmap="viridis", vmin=max(90, np.floor(vmin)), vmax=100)
    ax.set_xticks(range(n)); ax.set_yticks(range(n)); ax.set_xticklabels(order, rotation=90, fontsize=6); ax.set_yticklabels(order, fontsize=6)
    for lab in ax.get_xticklabels() + ax.get_yticklabels(): lab.set_color(GCOL[group_of(lab.get_text())])
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.set_label("16S identity (%)")
    save(fig, os.path.join(out, "Fig1_16S_identity_heatmap"))
    # ---- Fig 2: base per sequence at every diagnostic position
    if cols:
        fig, ax = plt.subplots(figsize=(max(5, 0.45 * len(cols) + 3), 0.22 * n + 1.8))
        for i, h in enumerate(order):
            for j, (g, c) in enumerate(cols):
                b = seqs[h][c]; ax.add_patch(Rectangle((j, n - 1 - i), 1, 1, fc=BCOL.get(b, "w"), ec="white", lw=0.8))
                ax.text(j + .5, n - 1 - i + .5, b, ha="center", va="center", fontsize=6)
                if group_of(h) == g: ax.add_patch(Rectangle((j, n - 1 - i), 1, 1, fc="none", ec="black", lw=1.0))
        ax.set_xlim(0, len(cols)); ax.set_ylim(0, n); ax.set_yticks([n - 1 - i + .5 for i in range(n)]); ax.set_yticklabels(order, fontsize=6)
        ax.set_xticks([j + .5 for j in range(len(cols))]); ax.set_xticklabels([str(int(cumref[ref["Group 1"]][c])) for g, c in cols], rotation=90, fontsize=6)
        for lab in ax.get_yticklabels(): lab.set_color(GCOL[group_of(lab.get_text())])
        for j, (g, c) in enumerate(cols): ax.add_patch(Rectangle((j, n + .1), 1, .5, fc=GCOL[g], ec="none", clip_on=False))
        ax.tick_params(length=0); [s.set_visible(False) for s in ax.spines.values()]
        save(fig, os.path.join(out, "Fig2_diagnostic_positions_per_sequence"))
    # ---- Fig 3 (one per group): sequence context around each diagnostic position
    for g in gs:
        pc = [c for gg, c in cols if gg == g]
        if not pc: continue
        W = 12; npn = len(pc); fig, axs = plt.subplots(1, npn, figsize=(2.6 * npn + 1.4, 0.2 * n + 1.2), squeeze=False)
        for ax, c in zip(axs[0], pc):
            lo, hi = max(0, c - W), min(L, c + W + 1)
            for i, h in enumerate(order):
                for j in range(lo, hi):
                    b = seqs[h][j]; ax.add_patch(Rectangle((j - lo, n - 1 - i), 1, 1, fc=BCOL.get(b, "w"), ec="white", lw=.5))
                    ax.text(j - lo + .5, n - 1 - i + .5, b, ha="center", va="center", fontsize=4.5, fontweight="bold" if j == c else "normal")
            ax.add_patch(Rectangle((c - lo, 0), 1, n, fc="none", ec="black", lw=1.2))
            ax.set_xlim(0, hi - lo); ax.set_ylim(0, n); ax.set_xticks([c - lo + .5]); ax.set_xticklabels([str(int(cumref[ref["Group 1"]][c]))], fontsize=7)
            ax.set_yticks([n - 1 - i + .5 for i in range(n)] if ax is axs[0][0] else []);
            if ax is axs[0][0]:
                ax.set_yticklabels(order, fontsize=5)
                for lab in ax.get_yticklabels(): lab.set_color(GCOL[group_of(lab.get_text())])
            ax.tick_params(length=0); [s.set_visible(False) for s in ax.spines.values()]
        save(fig, os.path.join(out, "Fig3_%s_diagnostic_context" % g.replace(" ", "")))
    # ---- Fig 4: unique windows summary: counts
    fig, ax = plt.subplots(figsize=(3.6, 2.8))
    cnt = [sum(1 for gg, c in cols if gg == g) for g in gs]
    ax.bar(range(len(gs)), cnt, color=[GCOL[g] for g in gs]); ax.set_xticks(range(len(gs))); ax.set_xticklabels(gs, fontsize=8)
    ax.set_ylabel("Diagnostic positions"); [ax.spines[s].set_visible(False) for s in ("top", "right")]
    save(fig, os.path.join(out, "Fig4_diagnostic_position_counts"))
    # README
    open(os.path.join(out, "README.md"), "w").write(
        "# 16S rRNA gene evidence for the four groups\n\n## Methods\n16S rRNA gene sequences of the 22 genomes were aligned with MAFFT. Identity between sequences was calculated over alignment columns without gaps. "
        "A position is group-diagnostic when at least %.0f%% of the sequences of the group carry the same base and no sequence of the other three groups carries that base. "
        "Each diagnostic position was checked in a 21-nt window against 16S sequences of described Candidatus Phytoplasma species from NCBI (up to %d per species).\n\n"
        "## Files\n- Table1_reference_16S: reference sequence per group (group_reference_16S.fasta)\n- Table2_identity_between_groups and Table5_identity_per_sequence\n"
        "- Table3_diagnostic_positions: positions, bases, 21-nt windows, result of the check against described species\n- Table4_base_per_sequence: base of every sequence at every diagnostic position\n"
        "- Fig1: identity heatmap; Fig2: base per sequence at the diagnostic positions (black outline = member of that group); Fig3_Group*: sequence context of each position; Fig4: counts\n"
        % (a.frac * 100, a.cap))
    print(T1.to_string(index=False)); print(T2.to_string(index=False)); print(T3.drop(columns=["bases_in_other_groups", "reference_21nt_window"]).to_string(index=False))
if __name__ == "__main__": main()
