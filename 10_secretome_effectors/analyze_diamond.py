# ANALYSIS: Secretome (SignalP, DIAMOND) and effector search in the 22 genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Fig. S6.2, Tables S8.x, S9.x, S10.x
# FOLDER: 10_secretome_effectors (see README.md in the parent folder)
import pandas as pd
from openpyxl.styles import Font
H = pd.read_csv("diamond/hits.tsv", sep="\t", names=["q","s","pid","len","qlen","slen","qcov","e","bits"])
H["qg"] = H.q.str.split("__").str[0]; H["sg"] = H.s.str.split("__").str[0]; H = H[H.qg != H.sg]
EO = pd.read_csv("results/effector_matrix.tsv", sep="\t").Strain.tolist(); PPS4 = EO[:4]; PPS6 = EO[:6]; REFS = EO[6:]
assert len(REFS) == 16
sec = pd.read_csv("results/Table_secretome_Secreted_proteins.tsv", sep="\t", dtype=str)
C = sec[sec.Prediction == "SP(Sec/SPI)"][["Genome","ID","Length_aa","Product"]].copy(); C["q"] = C.Genome + "__" + C.ID
g = H.groupby("q")
C["N_other_genomes_with_hit_of_21"] = C.q.map(g.sg.nunique()).fillna(0).astype(int)
C["N_references_with_hit_of_16_excluding_self"] = C.q.map(H[H.sg.isin(REFS)].groupby("q").sg.nunique()).fillna(0).astype(int)
C["Best_identity_pct"] = C.q.map(g.pid.max()).round(1)
C["Genomes_with_hit"] = C.q.map(g.sg.agg(lambda x: ", ".join(sorted(set(x), key=EO.index))).to_dict()).fillna("")
C["Genome_order"] = C.Genome.map({x: i for i, x in enumerate(EO)}); C = C.sort_values(["Genome_order", "ID"]).drop(columns=["q", "Genome_order"])
def summ(d):
    o = d.groupby("Genome").agg(Candidates=("ID", "size"), With_hit_in_other_genome=("N_other_genomes_with_hit_of_21", lambda x: int((x > 0).sum())), With_hit_in_reference=("N_references_with_hit_of_16_excluding_self", lambda x: int((x > 0).sum()))).reindex([x for x in EO if x in set(d.Genome)])
    o["No_hit_in_any_other_genome"] = o.Candidates - o.With_hit_in_other_genome; o["No_hit_in_any_reference"] = o.Candidates - o.With_hit_in_reference
    o["Pct_with_reference_hit"] = (100 * o.With_hit_in_reference / o.Candidates).round(1); o.index.name = "Genome"; return o.reset_index()
P4 = C[C.Genome.isin(PPS4)]; P6 = C[C.Genome.isin(PPS6)]
leg = pd.DataFrame([("Method", "DIAMOND (module diamond/2.1.13) BLASTp of all SP(Sec/SPI) candidates against the 22 proteomes, e-value 1e-5, minimum identity 30%, all hits kept. Hits to the query's own genome were removed."),
 ("Reference hit", "Hit in at least one of the 16 reference genomes other than the query genome (the 16 genomes after PPS1A, PPS2B, PPS6A, PPS6B, PPS9B, PPS10A)."),
 ("Four_PPS sheets", "Candidates of PPS1A, PPS2B, PPS6A, PPS6B only; hits are still searched in all other genomes.")], columns=["Item", "Definition"])
sh = {"Per_genome": summ(C), "Four_PPS_per_genome": summ(P4), "Candidates": C, "Four_PPS_candidates": P4, "PPS6_no_reference_hit": P6[P6.N_references_with_hit_of_16_excluding_self == 0], "Four_PPS_no_reference_hit": P4[P4.N_references_with_hit_of_16_excluding_self == 0], "Legend": leg}
for n, d in sh.items():
    if n != "Legend": d.to_csv("results/Table_secretome_conservation_" + n + ".tsv", sep="\t", index=False)
with pd.ExcelWriter("results/Table_secretome_conservation.xlsx", engine="openpyxl") as w:
    for n, d in sh.items():
        d.to_excel(w, sheet_name=n, index=False); ws = w.sheets[n]
        for row in ws.iter_rows():
            for c in row: c.font = Font(name="Arial", color="000000", bold=(c.row == 1))
        for i in range(1, ws.max_column + 1): ws.column_dimensions[ws.cell(1, i).column_letter].width = 70 if n == "Legend" else 18
        ws.freeze_panes = "B2"
print(summ(C).to_string(index=False)); print(sh["PPS6_no_reference_hit"][["Genome","ID","Length_aa","N_other_genomes_with_hit_of_21","Genomes_with_hit","Product"]].to_string(index=False))
