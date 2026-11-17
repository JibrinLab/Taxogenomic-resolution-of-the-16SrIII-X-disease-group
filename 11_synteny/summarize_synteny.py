# ANALYSIS: Whole-genome synteny of the PPS strains against four reference genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Figs. S1.4 to S1.7
# FOLDER: 11_synteny (see README.md in the parent folder)
import re, glob, os, pandas as pd
from openpyxl.styles import Font
rows = []
num = lambda s: float(re.match(r"[\d.]+", s).group())
pc = lambda s: float(re.search(r"\(([\d.]+)%", s).group(1))
for f in sorted(glob.glob("dnadiff/*.report")):
    ref, g = os.path.basename(f)[:-7].split("_vs_"); d = {"Reference": ref, "Genome": g}; sec = sub = None
    for l in open(f):
        t = l.split()
        if not t: continue
        if t[0].startswith("["): sec = t[0]
        if t[0] in ("1-to-1", "M-to-M"): sub = t[0].replace("-to-", "to")
        if sec == "[Sequences]" and t[0] == "TotalSeqs": d["Ref_contigs"], d["Query_contigs"] = int(t[1]), int(t[2])
        elif sec == "[Bases]" and t[0] == "TotalBases": d["Ref_bp"], d["Query_bp"] = int(t[1]), int(t[2])
        elif sec == "[Bases]" and t[0] == "AlignedBases": d["Ref_aligned_bp"], d["Query_aligned_bp"], d["Pct_ref_aligned"], d["Pct_query_aligned"] = int(num(t[1])), int(num(t[2])), pc(t[1]), pc(t[2])
        elif sec == "[Alignments]" and sub and t[0] == "AvgIdentity": d["AvgIdentity_" + sub] = float(t[1])
        elif sec == "[Alignments]" and sub and t[0] == "TotalLength": d["Query_alignedlength_" + sub] = int(t[2])
        elif sec == "[Feature" and t[0] in ("Breakpoints", "Relocations", "Translocations", "Inversions"): d[t[0] + "_query"] = int(t[2])
        elif sec == "[SNPs]" and t[0] == "TotalSNPs": d["SNPs"] = int(t[1])
    rows.append(d)
A = pd.DataFrame(rows); A.to_csv("synteny_all_pairs.tsv", sep="\t", index=False)
EO = pd.read_csv("../../06_secretome_effector_22genomes/results/effector_matrix.tsv", sep="\t").Strain.tolist(); REP = ["CanS4", "CX", "CicWB-2022", "LWB"]
pv = lambda c: A.pivot(index="Genome", columns="Reference", values=c).reindex(EO)[REP]
idm, alm, pct = pv("AvgIdentity_1to1"), pv("Query_aligned_bp"), pv("Pct_query_aligned")
cl = pd.DataFrame({"Closest_by_identity": idm.idxmax(axis=1), "Identity_pct": idm.max(axis=1), "Closest_by_aligned_bp": alm.idxmax(axis=1), "Aligned_bp": alm.max(axis=1)})
mn = lambda S: pd.DataFrame({"Mean_identity_pct": idm.loc[S].mean().round(2), "Mean_query_aligned_bp": alm.loc[S].mean().round(0), "Mean_pct_query_aligned": pct.loc[S].mean().round(1)}).reset_index()
leg = pd.DataFrame([("Method", "nucmer and dnadiff from MUMmer 4.0.0rc1. The representative is the reference and each genome is the query. Identity is dnadiff AvgIdentity of the 1-to-1 alignments; aligned bp and percent are the dnadiff AlignedBases of the query.")], columns=["Item", "Definition"])
sh = {"All_pairs": A, "Identity_1to1": idm.reset_index(), "Query_aligned_bp": alm.reset_index(), "Pct_query_aligned": pct.reset_index(), "Closest_representative": cl.reset_index(), "Mean_six_PPS": mn(EO[:6]), "Mean_four_PPS": mn(EO[:4]), "Legend": leg}
for n, d in sh.items():
    if n != "Legend": d.to_csv("Table_synteny_" + n + ".tsv", sep="\t", index=False)
with pd.ExcelWriter("Table_synteny_4reps.xlsx", engine="openpyxl") as w:
    for n, d in sh.items():
        d.to_excel(w, sheet_name=n, index=False); ws = w.sheets[n]
        for row in ws.iter_rows():
            for c in row: c.font = Font(name="Arial", color="000000", bold=(c.row == 1))
        for i in range(1, ws.max_column + 1): ws.column_dimensions[ws.cell(1, i).column_letter].width = 80 if n == "Legend" else 18
print(cl.round(2).to_string()); print(mn(EO[:6]).to_string(index=False)); print(mn(EO[:4]).to_string(index=False))
