# ANALYSIS: Secretome (SignalP, DIAMOND) and effector search in the 22 genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Fig. S6.2, Tables S8.x, S9.x, S10.x
# FOLDER: 10_secretome_effectors (see README.md in the parent folder)
import glob, os, pandas as pd
S, P = "signalp_out/", "proteomes22/"
rows, SQ = [], {}
for f in sorted(glob.glob(S + "*_signalp_summary.signalp5")):
    g = os.path.basename(f).replace("_signalp_summary.signalp5", "")
    fa = P + g + ".faa"; assert os.path.exists(fa), fa
    ln, pr, cur = {}, {}, None
    for l in open(fa):
        l = l.rstrip()
        if l.startswith(">"):
            cur = l[1:].split()[0]; pr[cur] = l.split(" ", 1)[1] if " " in l else ""; ln[cur] = 0; SQ[(g, cur)] = ""
        else: ln[cur] += len(l); SQ[(g, cur)] += l
    d = pd.read_csv(f, sep="\t", comment="#", header=None, names=["ID","Prediction","SP_prob","TAT_prob","LIPO_prob","OTHER_prob","CS_position"], dtype=str).fillna("")
    d["Genome"] = g; d["Length_aa"] = d.ID.map(ln); d["Product"] = d.ID.map(pr); rows.append(d)
A = pd.concat(rows, ignore_index=True)
from openpyxl.styles import Font
EO = pd.read_csv("results/effector_matrix.tsv", sep="\t").Strain.tolist(); PPS4 = EO[:4]; assert PPS4 == ["PPS1A", "PPS2B", "PPS6A", "PPS6B"]
cnt = A.groupby(["Genome", "Prediction"]).size().unstack(fill_value=0)
cnt.insert(0, "Proteins", A.groupby("Genome").size())
cnt["Secreted_SP_Sec_SPI"] = cnt.get("SP(Sec/SPI)", 0)
cnt["Secreted_pct_of_proteins"] = (100 * cnt.Secreted_SP_Sec_SPI / cnt.Proteins).round(1)
cnt = cnt.reindex(EO); cnt.index.name = "Genome"; cnt = cnt.reset_index()
sec = A[A.Prediction != "OTHER"].copy(); sec["o"] = sec.Genome.map({g: i for i, g in enumerate(EO)})
sec = sec.sort_values(["o", "ID"]).drop(columns="o")[["Genome","ID","Prediction","SP_prob","TAT_prob","LIPO_prob","CS_position","Length_aa","Product"]]
with open("diamond/candidates_SPSec.faa", "w") as o:
    for r in sec[sec.Prediction == "SP(Sec/SPI)"].itertuples(): o.write(f">{r.Genome}__{r.ID}\n{SQ[(r.Genome, r.ID)]}\n")
leg = pd.DataFrame([("Method", "SignalP 5.0b, organism type gram-positive (-org gram+), on the Prokka proteome of each genome. Secreted candidates are proteins with predicted class SP(Sec/SPI)."), ("Counts", "Number of proteins per SignalP class; Secreted_pct_of_proteins = SP(Sec/SPI) / all proteins."), ("Four_PPS sheets", "PPS1A, PPS2B, PPS6A, PPS6B.")], columns=["Item", "Definition"])
S4 = lambda d, c: d[d[c].isin(PPS4)]
sh = {"Counts": cnt, "Four_PPS_counts": S4(cnt, "Genome"), "Secreted_proteins": sec, "Four_PPS_secreted": S4(sec, "Genome"), "Legend": leg}
for n, d in sh.items():
    if n != "Legend": d.to_csv("results/Table_secretome_" + n + ".tsv", sep="\t", index=False)
with pd.ExcelWriter("results/Table_secretome_22genomes.xlsx", engine="openpyxl") as w:
    for n, d in sh.items():
        d.to_excel(w, sheet_name=n, index=False); ws = w.sheets[n]
        for row in ws.iter_rows():
            for c in row: c.font = Font(name="Arial", color="000000", bold=(c.row == 1))
        for i in range(1, ws.max_column + 1): ws.column_dimensions[ws.cell(1, i).column_letter].width = 70 if n == "Legend" else 16
        ws.freeze_panes = "B2"
print(cnt.to_string(index=False)); print("candidates:", int((sec.Prediction == "SP(Sec/SPI)").sum()))
