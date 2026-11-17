# ANALYSIS: GTDB-Tk classification of 53 genomes
# MANUSCRIPT SECTION: ANI, dDDH and GTDB clustering
# USED FOR: Tables S5.1 to S5.3
# FOLDER: 06_GTDB-Tk (see README.md in the parent folder)
import pandas as pd
d = pd.read_csv("10_gtdbtk/gtdbtk_out/gtdbtk.bac120.summary.tsv", sep="\t")
d["species"] = d["classification"].str.split(";s__").str[-1].replace("", "no species name")
ref = "closest_genome_reference" if "closest_genome_reference" in d.columns else None
cols = ["user_genome", "species"] + ([ref] if ref else []) + ["closest_genome_ani", "closest_genome_af", "classification_method"]
if "red_value" in d.columns: cols.append("red_value")
grp = {"Group 1": ["PPS1A","PPS2B","PPS6A","PPS6B","PPS9B","PPS10A","CanS4","MW1","VAC","MA1"],
       "Group 2": ["CX","6A1","ZT3-1","PR2021","JR1","PDA15"],
       "Group 3": ["CicWB-2022","ChTDIII","OVSP_L01","OVSP_S01","Vc33"], "Group 4": ["LWB"]}
g = {k: n for n, v in grp.items() for k in v}
d["ANI_group"] = d["user_genome"].map(g).fillna("other (31 added genomes)")
t = d[cols + ["ANI_group"]].rename(columns={"user_genome": "genome", "closest_genome_reference": "GTDB_reference_genome",
    "closest_genome_ani": "ANI_to_reference_pct", "closest_genome_af": "alignment_fraction"})
t.to_csv("10_gtdbtk/figures/GTDBTk_53genomes_table.tsv", sep="\t", index=False)
t[t.ANI_group != "other (31 added genomes)"].to_csv("10_gtdbtk/figures/GTDBTk_22genomes_table.tsv", sep="\t", index=False)
with pd.ExcelWriter("10_gtdbtk/figures/GTDBTk_tables.xlsx") as w:
    t[t.ANI_group != "other (31 added genomes)"].to_excel(w, "22_genomes", index=False)
    t.to_excel(w, "53_genomes", index=False)
print(t[t.ANI_group != "other (31 added genomes)"].sort_values(["species","ANI_to_reference_pct"], ascending=[True,False]).to_string(index=False))
