# ANALYSIS: GTDB-Tk classification of 53 genomes
# MANUSCRIPT SECTION: ANI, dDDH and GTDB clustering
# USED FOR: Tables S5.1 to S5.3
# FOLDER: 06_GTDB-Tk (see README.md in the parent folder)
import pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
t = pd.read_csv("10_gtdbtk/figures/GTDBTk_22genomes_table.tsv", sep="\t")
t["genome"] = t["genome"].str.replace("_", "-")
col = {"Group 1": "#4C78A8", "Group 2": "#F58518", "Group 3": "#54A24B", "Group 4": "#B279A2"}
order = {"Phytoplasma sp000309485": 0, "Phytoplasma pruni": 1, "Phytoplasma sp001623385": 2}
t["o"] = t["species"].map(order).fillna(3)
t = t.sort_values(["o", "ANI_to_reference_pct"], ascending=[True, False])
fig, ax = plt.subplots(figsize=(8.5, 11))
y, ticks = 0.0, []
for o, s in t.groupby("o"):
    ref = s["GTDB_reference_genome"].iloc[0] if "GTDB_reference_genome" in s and isinstance(s["GTDB_reference_genome"].iloc[0], str) else ""
    head = (s.species.iloc[0] + (f" (reference {ref})" if ref else "")) if o < 3 else "No GTDB species name"
    ax.text(93.9, y, head, fontsize=10, fontweight="bold", color="black", va="center", ha="left", bbox=dict(facecolor="white", edgecolor="none", pad=2), zorder=5); y += 1.4
    for _, r in s.iterrows():
        if pd.notna(r.ANI_to_reference_pct):
            ax.hlines(y, 94.0, r.ANI_to_reference_pct, color="black", lw=0.6)
            ax.scatter(r.ANI_to_reference_pct, y, s=80, color=col[r.ANI_group], edgecolor="black", zorder=3)
            ax.text(101.7, y, f"{r.ANI_to_reference_pct:.2f}", fontsize=9, color="black", va="center", ha="right")
        else:
            ax.scatter(94.0, y, s=80, color=col[r.ANI_group], edgecolor="black", zorder=3, marker="s")
            ax.text(101.7, y, "none*", fontsize=9, color="black", va="center", ha="right")
        ticks.append((y, r.genome)); y += 1.3
    y += 1.2
ax.text(101.7, -1.6, "ANI (%)", fontsize=9, fontweight="bold", color="black", ha="right", va="center")
ax.axvline(95, color="black", ls="--", lw=0.8)
ax.set_yticks([a for a, _ in ticks]); ax.set_yticklabels([b for _, b in ticks], fontsize=9, color="black")
ax.set_ylim(y, -2.4); ax.set_xlim(93.9, 101.8); ax.set_xticks([94, 95, 96, 97, 98, 99, 100])
ax.set_xlabel("ANI to the GTDB species reference genome (%)", color="black")
h = [Line2D([], [], marker="o", ls="", markerfacecolor=c, markeredgecolor="black") for k, c in col.items() if k != "LWB"]
h += [Line2D([], [], marker="s", ls="", markerfacecolor=col["Group 4"], markeredgecolor="black"), Line2D([], [], color="black", ls="--")]
ax.legend(h, ["Group 1", "Group 2", "Group 3", "Group 4 (no ANI value)", "95% ANI"], loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3, frameon=False, fontsize=9)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
fig.text(0.01, 0.005, "*LWB was placed in the GTDB tree and has no ANI value to a GTDB reference genome.", fontsize=8, color="black")
plt.tight_layout(rect=(0, 0.02, 1, 0.97))
for e in ("png", "pdf"): fig.savefig(f"10_gtdbtk/figures/gtdbtk_22genomes_ANI_to_reference.{e}", dpi=600, facecolor="white")
