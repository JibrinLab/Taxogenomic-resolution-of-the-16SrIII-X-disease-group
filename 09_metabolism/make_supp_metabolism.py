# ANALYSIS: Metabolic gene panel (66 genes, 7 pathways) in the 22 genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Fig. 8, Tables S7.1 to S7.10
# FOLDER: 09_metabolism (see README.md in the parent folder)
import csv
from openpyxl import load_workbook
from openpyxl.styles import Font, Border, Side
R = "results/"
rd = lambda p: list(csv.reader(open(p), delimiter="\t"))
GROUPS = {"Group 1": "PPS1A PPS2B PPS6A PPS6B PPS9B PPS10A CanS4 MW1 VAC MA1", "Group 2": "CX 6A1 ZT3-1 PR2021 JR1 PDA15", "Group 3": "CicWB-2022 ChTDIII OVSP-L01 OVSP-S01 Vc33", "Group 4": "LWB"}
allr = rd(R + "blastp_rescue_all.tsv"); hdr, body = allr[0], allr[1:]
mat = rd(R + "metabolic_final22_matrix.tsv"); labs = mat[0][2:]; genes = [(r[0], r[1]) for r in mat[1:]]
code = {}
for r in body:
    b = r[14].strip(); code[(r[2], r[0])] = "P" if b == "Prokka" else "B" if b == "BLASTp only" else "x" if b.startswith("BLASTp hit") else "-"
def tab(ls):
    T = [["Pathway", "Gene"] + ls + ["n Prokka (P)", "n BLASTp rescue (B)", "n P+B"]]
    for pw, g in genes:
        c = [code[(g, l)] for l in ls]; T.append([pw, g] + c + [c.count("P"), c.count("B"), c.count("P") + c.count("B")])
    T.append(["Genes detected, Prokka only (P)", ""] + [sum(1 for _, g in genes if code[(g, l)] == "P") for l in ls])
    T.append(["Genes detected, Prokka plus BLASTp (P+B)", ""] + [sum(1 for _, g in genes if code[(g, l)] in ("P", "B")) for l in ls])
    return T
GP = [["Genome", "Group", "Genes detected, Prokka (P)", "Genes detected, P+B", "Genes in panel"]]
for gr, v in GROUPS.items():
    for x in v.split(): GP.append([x, gr, sum(1 for _, g in genes if code[(g, x)] == "P"), sum(1 for _, g in genes if code[(g, x)] in ("P", "B")), len(genes)])
GL = [["Gene", "Pathway", "EC numbers", "Gene-name matches", "Product keywords"]] + [l.rstrip("\n").split("|") for l in open("scripts/metab_genes.txt") if l.strip()]
LG = [["Supplementary Table: pathway gene content of 22 Ca. Phytoplasma pruni and related genomes (66 genes, 7 pathways)"],
 ["Calls: P = detected in the Prokka v1.13 GFF by EC number, gene name or product keyword. B = not detected by Prokka; BLASTp hit (e-value <= 1e-3) to a protein annotated as a hypothetical protein. x = BLASTp hit to a protein with a different annotated function, or already assigned to another gene (not counted). - = no detection."],
 ["BLASTp: v2.16.0, -seg no -comp_based_stats 0, e-value 1e-3. Queries: the proteins of the same gene in seven annotated reference phytoplasma genomes (AT, AYWB, FD, M8 [directory OYM_Ca.P.asteris_16SrI], NZSb11, TB2022, GOE) plus UniProt SecE (Q6YQW9, Q2NJ09). Targets: the Prokka proteomes of the 22 genomes."],
 ["Quality: strong = identity >= 30% and query coverage >= 50%; weak = otherwise (criterion of this analysis). Counting rule for B: hit protein annotated as hypothetical (criterion of this analysis)."],
 ["Absence means no detection by these methods and does not prove the gene is missing. Nine genes (folB, secG, lepB, troC, troD, oppC, lysX, efrA, efrB) had no phytoplasma query protein, so only Prokka detection applies to them."]]
wb = load_workbook(R + "metabolism_22genomes_presence_absence.xlsx"); del wb["BLASTp_search_hits"]
new = [("Legend", LG), ("Gene_list", GL), ("Calls_matrix", tab(labs)), ("PPS4", tab(["PPS1A", "PPS2B", "PPS6A", "PPS6B"])), ("Genes_per_genome", GP), ("Prokka_evidence", rd(R + "metabolic_all22_evidence.tsv")), ("BLASTp_all_calls", allr), ("BLASTp_candidates", [hdr] + [r for r in body if r[3] == "0" and r[4] == "1"])]
for name, T in new:
    ws = wb.create_sheet(name)
    for r in T: ws.append(r)
    if name != "Legend": csv.writer(open(R + "sheet_" + name + ".tsv", "w"), delimiter="\t", lineterminator="\n").writerows(T)
for ws in wb.worksheets:
    for row in ws.iter_rows():
        for c in row: c.font = Font(name="Arial", size=10, bold=(c.row == 1))
    for c in ws[1]: c.border = Border(bottom=Side(style="thin"))
    if ws.title != "Legend":
        ws.freeze_panes = "C2"
        for col in ws.columns: ws.column_dimensions[col[0].column_letter].width = min(45, max(len(str(c.value)) for c in col) + 2)
wb._sheets = [wb[n] for n in ["Legend", "Gene_list", "Calls_matrix", "PPS4", "Genes_per_genome", "Presence_absence", "Set_specific_or_absent", "Summary_counts", "Prokka_evidence", "BLASTp_all_calls", "BLASTp_candidates"]]
wb.save(R + "Supplementary_Table_metabolism_22genomes.xlsx")
print([ws.title for ws in wb.worksheets])
