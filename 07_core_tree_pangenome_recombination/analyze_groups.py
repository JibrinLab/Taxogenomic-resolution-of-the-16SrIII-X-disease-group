#!/usr/bin/env python3
# ANALYSIS: Core genome tree, pangenome and recombination
# MANUSCRIPT SECTION: Phylogenetic reconstructions; Comparative Genomics Analyses
# USED FOR: Fig. 5, Fig. 6, Figs. S5.1 to S5.6
# FOLDER: 07_core_tree_pangenome_recombination (see README.md in the parent folder)
"""Pangenome openness (Heaps' law) and ClonalFrameML summary for Groups 1-3.
Run in the 15_group_pangenome_recombination folder:  python3 analyze_groups.py
Needs numpy, pandas, scipy-free, matplotlib, openpyxl (conda env tree_viz)."""
import os, re, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

W = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()
COL = {1: "#C1440E", 2: "#2E6DB4", 3: "#2A8F6F"}
NPERM = 1000
rng = np.random.default_rng(1)
plt.rcParams.update({"font.family": "DejaVu Sans", "text.color": "black", "axes.labelcolor": "black",
                     "xtick.color": "black", "ytick.color": "black", "axes.edgecolor": "black",
                     "pdf.fonttype": 42, "svg.fonttype": "none"})

def curves(M):
    n = M.shape[0]; pan = np.zeros((NPERM, n)); core = np.zeros((NPERM, n))
    k = np.arange(1, n + 1)[:, None]
    for p in range(NPERM):
        cum = np.cumsum(M[rng.permutation(n)], axis=0)
        pan[p] = (cum > 0).sum(1); core[p] = (cum == k).sum(1)
    return pan, core

def heaps(new, Ns):
    """new = genes added at genome N (N>=2); fit new = kappa * N^-alpha on log scale."""
    m = new > 0
    if m.sum() < 2: return np.nan, np.nan
    b, a = np.polyfit(np.log(Ns[m]), np.log(new[m]), 1)
    return -b, np.exp(a)

def pan_stats(g):
    f = os.path.join(W, "Group%d" % g, "panaroo_out", "gene_presence_absence.Rtab")
    df = pd.read_csv(f, sep="\t", index_col=0)
    M = (df.T.values > 0)
    M = M[:, M.sum(0) > 0]
    n = M.shape[0]; Ns = np.arange(2, n + 1)
    pan, core = curves(M)
    new = np.diff(pan, axis=1)                      # genes added at N = 2..n
    a_mean, k_mean = heaps(new.mean(0), Ns)
    al = np.array([heaps(new[p], Ns)[0] for p in range(NPERM)]); al = al[~np.isnan(al)]
    lo, hi = (np.percentile(al, [2.5, 97.5]) if len(al) else (np.nan, np.nan))
    cnt = M.sum(0)
    st = dict(group=g, genomes=n, total_genes=M.shape[1], core_genes=int((cnt == n).sum()),
              accessory_genes=int((cnt < n).sum()), unique_genes=int((cnt == 1).sum()),
              alpha=a_mean, alpha_lo=lo, alpha_hi=hi, gamma=1 - a_mean,
              pangenome="open" if a_mean < 1 else "closed")
    return st, pan, core, list(df.columns)

def parse_em(path):
    t = open(path).read()
    def get(name):
        m = re.search(r"^" + re.escape(name) + r"\s+([-\d.eE+]+)", t, re.M)
        return float(m.group(1)) if m else np.nan
    return dict(R_theta=get("R/theta"), inv_delta=get("1/delta"), nu=get("nu"))

def cfml(g):
    d = os.path.join(W, "Group%d" % g, "clonalframe")
    em = parse_em(os.path.join(d, "cfml.em.txt"))
    imp = pd.read_csv(os.path.join(d, "cfml.importation_status.txt"), sep="\t")
    imp.columns = [c.strip() for c in imp.columns]
    aln_len = sum(len(l.strip()) for l in open(os.path.join(d, "core.ACGTN.aln")).read().split(">")[1].split("\n")[1:])
    tl = open(os.path.join(d, "cfml.labelled_tree.newick")).read()
    order = re.findall(r"([A-Za-z0-9_.\-]+)(?=:)", tl)
    delta = 1.0 / em["inv_delta"] if em["inv_delta"] else np.nan
    em["delta_bp_formula"] = delta
    em["r_m"] = em["R_theta"] * delta * em["nu"]
    em["events"] = len(imp)
    em["mean_tract_bp"] = float((imp["End"] - imp["Beg"] + 1).mean()) if len(imp) else np.nan
    em["aln_len"] = aln_len
    return em, imp, order

def panel_pan(ax, g, pan, core, st):
    n = pan.shape[1]; x = np.arange(1, n + 1); c = COL[g]
    for arr, ls, lab in ((pan, "-", "Pan-genome"), (core, "--", "Core genome")):
        mu = arr.mean(0); ax.fill_between(x, np.percentile(arr, 2.5, axis=0), np.percentile(arr, 97.5, axis=0), color=c, alpha=0.15, lw=0)
        ax.plot(x, mu, ls, color=c, lw=1.8, marker="o", ms=4, label=lab)
    ax.set_xlabel("Number of genomes"); ax.set_ylabel("Number of genes"); ax.set_xticks(x)
    ax.legend(frameon=False, fontsize=8, loc="center right")
    ax.text(0.03, 0.97, u"Heaps' law α = %.2f\n%s pan-genome (%d genes, core %d)" %
            (st["alpha"], st["pangenome"], st["total_genes"], st["core_genes"]),
            transform=ax.transAxes, va="top", fontsize=8)
    for s in ("top", "right"): ax.spines[s].set_visible(False)

def panel_rec(ax, g, imp, order, em):
    tips = [o for o in order if not o.startswith("NODE")]
    nodes = [o for o in order if o in set(imp["Node"]) or o in tips]
    y = {nm: i for i, nm in enumerate(nodes)}
    for nm, i in y.items(): ax.axhline(i, color="#DDDDDD", lw=0.4, zorder=0)
    s = imp[imp["Node"].isin(y)]
    ax.barh([y[v] for v in s["Node"]], s["End"] - s["Beg"] + 1, left=s["Beg"], height=0.8, color=COL[g], lw=0)
    ax.set_yticks([y[t] for t in tips if t in y]); ax.set_yticklabels([t for t in tips if t in y], fontsize=7)
    ax.set_xlim(0, em["aln_len"]); ax.set_ylim(-0.7, len(nodes) - 0.3)
    ax.set_xlabel("Position in core alignment (bp)"); ax.set_ylabel("Branches (tips labelled)")
    ax.text(0.99, 1.01, u"r/m = %.3f   R/θ = %.3f   %d events" % (em["r_m"], em["R_theta"], em["events"]),
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)

def save(fig, name):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(os.path.join(W, "figures", "%s.%s" % (name, ext)), dpi=600, bbox_inches="tight")
    plt.close(fig)

def main():
    os.makedirs(os.path.join(W, "figures"), exist_ok=True)
    rows, data = [], {}
    for g in (1, 2, 3):
        st, pan, core, gen = pan_stats(g); em, imp, order = cfml(g)
        data[g] = (st, pan, core, em, imp, order, gen)
        rows.append({**st, **{k: v for k, v in em.items()}, "genome_list": ", ".join(gen)})
        chk = abs(em["mean_tract_bp"] - em["delta_bp_formula"]) / em["delta_bp_formula"] if em["delta_bp_formula"] else np.nan
        print("Group %d: r/m=%.4f  R/theta=%.4f  1/delta(printed)=%g  nu=%.4f  mean tract from events=%.0f bp  1/(1/delta)=%.0f bp  (relative diff %.2f)" %
              (g, em["r_m"], em["R_theta"], em["inv_delta"], em["nu"], em["mean_tract_bp"], em["delta_bp_formula"], chk))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(W, "group_pangenome_recombination_summary.tsv"), sep="\t", index=False)
    df.to_excel(os.path.join(W, "group_pangenome_recombination_summary.xlsx"), index=False)
    for g in (1, 2, 3):                                  # one figure per group
        st, pan, core, em, imp, order, gen = data[g]
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5.5), gridspec_kw={"width_ratios": [1, 1.4]})
        panel_pan(a1, g, pan, core, st); panel_rec(a2, g, imp, order, em)
        fig.suptitle("Group %d (n = %d genomes)" % (g, st["genomes"]), fontweight="bold", y=1.0)
        fig.tight_layout(); save(fig, "Fig_Group%d_pangenome_recombination" % g)
    fig, axs = plt.subplots(2, 3, figsize=(18, 10), gridspec_kw={"height_ratios": [1, 1.5]})
    for j, g in enumerate((1, 2, 3)):
        st, pan, core, em, imp, order, gen = data[g]
        panel_pan(axs[0, j], g, pan, core, st); panel_rec(axs[1, j], g, imp, order, em)
        axs[0, j].set_title("Group %d (n = %d)" % (g, st["genomes"]), fontweight="bold")
    fig.tight_layout(); save(fig, "Fig_Groups1-3_pangenome_recombination")
    fig, axs = plt.subplots(1, 3, figsize=(10, 3.6))
    for ax, key, lab in zip(axs, ("r_m", "R_theta", "mean_tract_bp"), ("r/m", u"R/θ", "Mean import length (bp)")):
        v = [data[g][3][key] for g in (1, 2, 3)]
        ax.bar(["Group 1", "Group 2", "Group 3"], v, color=[COL[g] for g in (1, 2, 3)], width=0.6)
        for i, x in enumerate(v): ax.text(i, x, ("%.0f" % x if x >= 10 else "%.3g" % x), ha="center", va="bottom", fontsize=8)
        ax.set_ylabel(lab)
        for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    fig.tight_layout(); save(fig, "Fig_Groups1-3_recombination_rates")
    print(df[["group", "genomes", "total_genes", "core_genes", "alpha", "alpha_lo", "alpha_hi", "pangenome", "R_theta", "inv_delta", "nu", "r_m", "events"]].to_string(index=False))

if __name__ == "__main__": main()
