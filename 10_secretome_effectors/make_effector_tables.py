# ANALYSIS: Secretome (SignalP, DIAMOND) and effector search in the 22 genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Fig. S6.2, Tables S8.x, S9.x, S10.x
# FOLDER: 10_secretome_effectors (see README.md in the parent folder)
import pandas as pd
from openpyxl.styles import PatternFill, Font, Alignment
R = "results/"
M = pd.read_csv(R+"effector_matrix.tsv", sep="\t", dtype=str).fillna("")
D = pd.read_csv(R+"effector_matrix_detail.tsv", sep="\t", dtype=str).fillna("")
L = pd.read_csv(R+"effector_calls_long.tsv", sep="\t", dtype=str).fillna("")
PPS4 = ["PPS1A", "PPS2B", "PPS6A", "PPS6B"]
leg = pd.DataFrame([
("Strain", "Genome name. Main sheets: the 22 pruni and pruni-related genomes. Four_PPS sheets: PPS1A, PPS2B, PPS6A, PPS6B."),
("Search", "Query proteins were searched by BLASTp (BLAST+ 2.16.0) against Prokka proteomes and by TBLASTN (genetic code 4, TGA not counted as a stop) against the assemblies. Only hits with e-value <= 1e-3 were used."),
("present", "Identity >= 60% over >= 80% of the query, one reading frame, no internal stop."),
("frameshift", "As present, but the alignment spans more than one reading frame."),
("stop codon", "Identity >= 60% over >= 50% of the query with an internal TAA or TAG, or adjacent same-frame segments separated by <= 6 nt."),
("contig-break", "Identity >= 60% over < 80% of the query, and the missing part of the query lies beyond the contig end (distance to the contig end < 3 x missing aa + 30 nt)."),
("partial", "Identity >= 60% over < 80% of the query, not explained by a contig end."),
("paralog", "Identity 35 to < 60% over >= 70% of the query."),
("absent", "No hit meeting the criteria above."),
("Detail sheets", "Call (identity of the best locus; +N call = additional loci with that call). Identity is from TBLASTN."),
("Thresholds", "All identity, coverage and distance thresholds were set by the analyst and are not taken from a publication."),
], columns=["Item", "Definition"])
C = {"present":"C6E0B4","frameshift":"FFE699","stop codon":"F4B084","contig-break":"BDD7EE","partial":"FFF2CC","paralog":"E4DFEC","absent":"D9D9D9"}
sheets = {"Matrix":M, "Matrix_detail":D, "Four_PPS_matrix":M[M.Strain.isin(PPS4)], "Four_PPS_detail":D[D.Strain.isin(PPS4)], "All_calls_long":L, "Four_PPS_calls_long":L[L.Genome.isin(PPS4)], "Legend":leg}
for n, d in sheets.items():
    if n != "Legend": d.to_csv(R+"Table_effectors_"+n+".tsv", sep="\t", index=False)
with pd.ExcelWriter(R+"Table_effectors_22genomes.xlsx", engine="openpyxl") as w:
    for n, d in sheets.items():
        d.to_excel(w, sheet_name=n, index=False); ws = w.sheets[n]
        for c in ws[1]: c.font = Font(name="Arial", bold=True, color="000000"); c.alignment = Alignment(wrap_text=True, vertical="top")
        for row in ws.iter_rows(min_row=2):
            for c in row:
                c.font = Font(name="Arial", color="000000"); v = str(c.value or "")
                if n in ("Matrix", "Matrix_detail", "Four_PPS_matrix", "Four_PPS_detail") and c.column > 1:
                    k = next((x for x in C if v.startswith(x)), None)
                    if k: c.fill = PatternFill("solid", start_color=C[k])
        for i in range(1, ws.max_column + 1): ws.column_dimensions[ws.cell(1, i).column_letter].width = 14 if "matrix" in n else (70 if n == "Legend" else 24)
        ws.freeze_panes = "B2"
print("written", R)
