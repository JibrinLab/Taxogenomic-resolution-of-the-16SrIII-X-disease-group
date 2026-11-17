# ANALYSIS: GTDB-Tk classification of 53 genomes
# MANUSCRIPT SECTION: ANI, dDDH and GTDB clustering
# USED FOR: Tables S5.1 to S5.3
# FOLDER: 06_GTDB-Tk (see README.md in the parent folder)
import pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
t = pd.read_csv("figures/GTDBTk_53genomes_table.tsv", sep="\t")
r = pd.read_csv("ANIm_cluster95_reps.tsv", sep="\t"); r = r[r["in22"].astype(str) == "-"]
four = {"Group 1": 1, "Group 2": 2, "Group 3": 3, "Group 4": 4}
rank = {c: i + 5 for i, c in enumerate(sorted(r["cluster"]))}
t["oc"] = t["genome"].map(dict(zip(r["rep"], r["cluster"])))
t["cl"] = [four.get(g, rank.get(c)) for g, c in zip(t["ANI_group"], t["oc"])]
print("rows", len(t), "unassigned", t["cl"].isna().sum(), t.loc[t["cl"].isna(), "genome"].tolist())
t = t.dropna(subset=["cl"]); t["cl"] = t["cl"].astype(int)
t["lab"] = [g.replace("_", "-") if k <= 4 else g for g, k in zip(t["genome"], t["cl"])]
t = t.sort_values(["cl", "ANI_to_reference_pct"], ascending=[True, False])
cm = plt.get_cmap("tab20").colors + plt.get_cmap("tab20b").colors; colr = {k: cm[i] for i, k in enumerate(sorted(t["cl"].unique()))}
lo = t["ANI_to_reference_pct"].min(); xmin = int(lo) - 1 if lo == lo else 90
fig, ax = plt.subplots(figsize=(10, 14)); y = 0.0; ticks = []
for k, s in t.groupby("cl"):
    ax.text(xmin, y, (s["species"].dropna().iloc[0] if s["species"].notna().any() else "No GTDB species name") + (f" ({len(s)} genomes)" if len(s) > 1 else ""), fontsize=9, fontweight="bold", color="black", va="center", bbox=dict(facecolor="white", edgecolor="none", pad=1.5), zorder=5); y += 1.3
    for _, q in s.iterrows():
        fc = colr.get(k, "#BDBDBD"); a = q["ANI_to_reference_pct"]
        if pd.notna(a): ax.hlines(y, xmin, a, color=fc, lw=1.4, zorder=2); ax.scatter(a, y, s=70, color=fc, edgecolor="black", zorder=3); txt = f"{a:.2f}"
        else: ax.scatter(xmin, y, s=70, color=fc, edgecolor="black", marker="s", zorder=3); txt = "none"
        sp = q["species"] if isinstance(q["species"], str) else "no GTDB species name"
        ax.text(1.01, y, txt, transform=ax.get_yaxis_transform(), fontsize=8, va="center", color="black", clip_on=False)
        ticks.append((y, q["lab"])); y += 1.0
    y += 0.8
ax.axvline(95, color="black", ls="--", lw=0.8)
ax.set_yticks([a for a, _ in ticks]); ax.set_yticklabels([b for _, b in ticks], fontsize=8, color="black")
ax.set_ylim(y, -1.5); ax.set_xlim(xmin, 100.5); ax.set_xlabel("ANI to the GTDB species reference genome (%)", color="black", fontsize=9)
ax.tick_params(colors="black"); [ax.spines[s].set_visible(False) for s in ("top", "right")]
fig.subplots_adjust(left=0.22, right=0.86)
for e in ("png", "pdf"): fig.savefig(f"figures/gtdbtk_53genomes_clusters.{e}", dpi=600, bbox_inches="tight", facecolor="white")
o = t[["lab", "genome", "cl", "oc", "species", "GTDB_reference_genome", "ANI_to_reference_pct"]]
o.columns = ["Genome", "Genome ID", "ANI cluster", "ANIm 95% cluster number (83 curated genomes)", "GTDB species", "GTDB reference genome", "ANI to reference (%)"]
o.to_csv("figures/GTDBTk_53genomes_clusters_table.tsv", sep="\t", index=False); o.to_excel("figures/GTDBTk_53genomes_clusters_table.xlsx", index=False)
