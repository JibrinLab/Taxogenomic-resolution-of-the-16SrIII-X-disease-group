# ANALYSIS: Average nucleotide identity (FastANI, 83 curated and 299 genomes; linkage sensitivity)
# MANUSCRIPT SECTION: Average nucleotide identity
# USED FOR: ANI levels n=83 and n=299; linkage sensitivity
# FOLDER: 02_ANI (see README.md in the parent folder)
import os, re, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
ROOT = "/projects/dsn011/ppss_publication/2026-10-20_new_species_taxonomy"
XL = ROOT + "/PPS_taxonomy_package_for_Dr_Jibrin_2026-10-22/01_ANI_tables_and_figures/pruni_ANI_tables.xlsx"
O1 = ROOT + "/20_ani_linkage_sensitivity"
O2 = ROOT + "/21_ani_distribution"
norm = lambda s: re.sub(r"[^A-Za-z0-9]", "", str(s)).upper()
strip = lambda s: re.sub(r"^GC[AF]_\d+\.\d+_", "", str(s))
G = {
 "PPS group": ["PPS1A","PPS2B","PPS6A","PPS6B","PPS9B","PPS10A","CanS4","MW1","VAC","MA1","MA"],
 "CX group": ["CX","6A1","ZT3-1","PR2021","JR1","PDA15"],
 "CicWB-2022 group": ["CicWB-2022","Vc33","ChTDIII","OVSP_L01","OVSP_S01"],
 "LWB": ["LWB"],
}
GRP = {norm(n): g for g, ns in G.items() for n in ns}
def group_of(label):
    k = norm(label)
    if k not in GRP:
        print("UNMAPPED GENOME LABEL:", label); sys.exit(1)
    return GRP[k]
def partition(labels, ids):
    d = {}
    for l, i in zip(labels, ids):
        d.setdefault(i, set()).add(l)
    return set(frozenset(v) for v in d.values())
res_rows, mem = [], {}
for sheet, name in [("ANIm_matrix", "ANIm"), ("ANIb_matrix", "ANIb")]:
    M = pd.read_excel(XL, sheet_name=sheet, index_col=0)
    M.index = [strip(x) for x in M.index]; M.columns = [strip(x) for x in M.columns]
    M = M.loc[M.index, M.index]
    A = M.values.astype(float)
    if np.nanmax(A) <= 1.0001: A = A * 100.0
    A = (A + A.T) / 2.0
    labels = list(M.index)
    four = partition(labels, [group_of(l) for l in labels])
    D = 100.0 - A; np.fill_diagonal(D, 0.0)
    cond = squareform(D, checks=False)
    for method in ["average", "single", "complete"]:
        Z = linkage(cond, method=method)
        for thr in [94, 95, 96]:
            ids = fcluster(Z, t=100 - thr, criterion="distance")
            part = partition(labels, ids)
            ma = [l for l in labels if norm(l) in ("MA1", "MA")]
            ma_cl = sorted(next((c for c in part if ma and ma[0] in c), set()))
            comp = " | ".join(", ".join(sorted(c)) for c in sorted(part, key=lambda c: sorted(c)[0]))
            res_rows.append(dict(matrix=name, linkage=method, threshold=thr, n_clusters=len(part),
                                 equals_four_groups=(part == four), MA1_cluster=", ".join(ma_cl), clusters=comp))
            mem[(name, method, thr)] = dict(zip(labels, ids))
summ = pd.DataFrame(res_rows)
summ.to_csv(O1 + "/linkage_sensitivity_summary.tsv", sep="\t", index=False)
pd.DataFrame(mem).to_csv(O1 + "/linkage_sensitivity_membership.tsv", sep="\t")
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 120)
print(summ.drop(columns=["clusters"]).to_string(index=False))
P = pd.read_excel(XL, sheet_name="pairs")
P["gA"] = P["A"].map(group_of); P["gB"] = P["B"].map(group_of)
P["relation"] = np.where(P["gA"] == P["gB"], "within group", "between groups")
P.to_csv(O2 + "/pairs_with_group_labels.tsv", sep="\t", index=False)
v = np.sort(P["ANIm"].values.astype(float)); v = v[(v >= 92) & (v <= 98)]
gaps = pd.DataFrame({"from": v[:-1], "to": v[1:], "gap": np.diff(v)}).sort_values("gap", ascending=False).head(6)
print("\nLargest gaps in pair ANIm (92 to 98):"); print(gaps.to_string(index=False))
gaps.to_csv(O2 + "/largest_ANIm_gaps_92_98.tsv", sep="\t", index=False)
cols = ["A","B","gA","gB","ANIm","ANIb","min_ANIm_aligned_%"]
btw = P[(P.relation == "between groups") & (P.ANIm >= 94.5)].sort_values("ANIm", ascending=False)
wit = P[(P.relation == "within group") & (P.ANIm < 96)].sort_values("ANIm")
print("\nBetween-group pairs with ANIm >= 94.5:"); print(btw[cols].to_string(index=False))
print("\nWithin-group pairs with ANIm < 96:"); print(wit[cols].to_string(index=False))
btw[cols].to_csv(O2 + "/between_group_pairs_ge_94.5.tsv", sep="\t", index=False)
wit[cols].to_csv(O2 + "/within_group_pairs_lt_96.tsv", sep="\t", index=False)
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
bins = np.arange(88, 100.26, 0.25)
w = P[P.relation == "within group"]["ANIm"]; b = P[P.relation == "between groups"]["ANIm"]
ax[0].hist([w, b], bins=bins, stacked=True, color=["#2F6DB5", "#C2410C"], label=["within group", "between groups"])
for x in (95, 96): ax[0].axvline(x, color="#444444", ls="--", lw=1)
ax[0].set_xlabel("ANIm (%)"); ax[0].set_ylabel("Number of genome pairs"); ax[0].legend(frameon=False)
ax[0].set_title("All 231 pairs", fontsize=10)
for rel, c in (("within group", "#2F6DB5"), ("between groups", "#C2410C")):
    s = P[P.relation == rel]
    ax[1].scatter(s["ANIm"], s["min_ANIm_aligned_%"], s=22, color=c, edgecolor="white", linewidth=0.5, label=rel)
for x in (95, 96): ax[1].axvline(x, color="#444444", ls="--", lw=1)
ax[1].set_xlabel("ANIm (%)"); ax[1].set_ylabel("Minimum aligned fraction (%)"); ax[1].legend(frameon=False)
ax[1].set_title("ANIm against aligned fraction", fontsize=10)
for a in ax:
    a.spines["top"].set_visible(False); a.spines["right"].set_visible(False)
plt.tight_layout()
for ext in ("png", "pdf", "svg"): plt.savefig(O2 + "/Fig_ANIm_distribution." + ext, dpi=200)
open(O1 + "/README.txt", "w").write(
"20_ani_linkage_sensitivity\n\nQUESTION: do the four ANI groups survive other clustering rules and thresholds?\n"
"INPUT: ANIm and ANIb matrices in 01_ANI_tables_and_figures/pruni_ANI_tables.xlsx (22 genomes), symmetrised.\n"
"METHOD: hierarchical clustering (scipy) on 100 minus ANI; average, single and complete linkage, cut at 94, 95 and 96% ANI.\n"
"OUTPUT: linkage_sensitivity_summary.tsv (clusters per setting; equals_four_groups is True when the clusters are exactly the four groups), linkage_sensitivity_membership.tsv.\n"
"Script: run_ani_tests.py\n")
open(O2 + "/README.txt", "w").write(
"21_ani_distribution\n\nQUESTION: is there a gap in pairwise ANI near 95%, and do alignments support the values?\n"
"INPUT: pairs sheet of pruni_ANI_tables.xlsx (231 pairs; ANIm, ANIb, minimum aligned fraction).\n"
"OUTPUT: Fig_ANIm_distribution (histogram within vs between groups; ANIm against aligned fraction), pairs_with_group_labels.tsv, largest_ANIm_gaps_92_98.tsv, between_group_pairs_ge_94.5.tsv, within_group_pairs_lt_96.tsv.\n"
"Script: ../20_ani_linkage_sensitivity/run_ani_tests.py\n")
print("\nDone.")
