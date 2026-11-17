# ANALYSIS: Metabolic gene panel (66 genes, 7 pathways) in the 22 genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Fig. 8, Tables S7.1 to S7.10
# FOLDER: 09_metabolism (see README.md in the parent folder)
import csv
from openpyxl import Workbook
from openpyxl.styles import Font, Border, Side
GROUPS = {"Group 1": "PPS1A PPS2B PPS6A PPS6B PPS9B PPS10A CanS4 MW1 VAC MA1",
          "Group 2": "CX 6A1 ZT3-1 PR2021 JR1 PDA15",
          "Group 3": "CicWB-2022 ChTDIII OVSP-L01 OVSP-S01 Vc33", "Group 4": "LWB"}
SETS = dict(GROUPS); SETS["PPS4 (PPS1A PPS2B PPS6A PPS6B)"] = "PPS1A PPS2B PPS6A PPS6B"
CL = ["Specific to set (absent in all other genomes)", "Absent from set (present in other genomes)"]
rd = csv.reader(open("results/metabolic_final22_matrix.tsv"), delimiter="\t"); h = next(rd)
labs = [x.replace("_Ca.P.pruni_16SrIII", "").replace("_", "-") for x in h[2:]]
rows = [(r[0], r[1], {l: float(x) > 0 for l, x in zip(labs, r[2:])}) for r in rd]
N = len(labs)
T1 = [["Pathway", "Gene", "Class (all 22)"] + labs + ["n present"] + ["%s (n=%d)" % (k, len(v.split())) for k, v in GROUPS.items()]]
for pw, g, p in rows:
    t = sum(p.values())
    c = "Present in all 22" if t == N else "Absent in all 22" if t == 0 else "Variable"
    T1.append([pw, g, c] + [int(p[l]) for l in labs] + [t] + [sum(p[x] for x in v.split()) for v in GROUPS.values()])
T2 = [["Set", "Class", "Pathway", "Gene", "Carrier genomes", "n in set", "n outside set"]]
for name, mem in SETS.items():
    m = mem.split()
    for pw, g, p in rows:
        i = sum(p[x] for x in m); o = sum(p.values()) - i
        if i > 0 and o == 0: T2.append([name, CL[0], pw, g, ",".join(x for x in m if p[x]), i, o])
        elif i == 0 and o > 0: T2.append([name, CL[1], pw, g, ",".join(x for x in labs if p[x]), i, o])
cnt = {}
for r in T2[1:]: cnt[(r[0], r[1])] = cnt.get((r[0], r[1]), 0) + 1
T3 = [["Set", "Class", "Number of genes"]]
for name in SETS:
    for c in CL: T3.append([name, c, cnt.get((name, c), 0)])
for c in ["Present in all 22", "Absent in all 22", "Variable"]: T3.append(["All 22 genomes", c, sum(1 for r in T1[1:] if r[2] == c)])
wb = Workbook(); wb.remove(wb.active)
for title, T, fn in [("Presence_absence", T1, "presence_absence_matrix"), ("Set_specific_or_absent", T2, "set_specific_or_absent_genes"), ("Summary_counts", T3, "summary_counts")]:
    csv.writer(open("results/%s.tsv" % fn, "w"), delimiter="\t").writerows(T)
    ws = wb.create_sheet(title)
    for r in T: ws.append(r)
    for row in ws.iter_rows():
        for c in row: c.font = Font(name="Arial", size=10, bold=(c.row == 1))
    for c in ws[1]: c.border = Border(bottom=Side(style="thin"))
    ws.freeze_panes = "C2"
    for col in ws.columns: ws.column_dimensions[col[0].column_letter].width = min(45, max(len(str(c.value)) for c in col) + 2)
rr = list(csv.reader(open("results/blastp_rescue_all.tsv"), delimiter="\t"))
ws = wb.create_sheet("BLASTp_search_hits")
for r in [rr[0]] + [x for x in rr[1:] if x[3] == "0" and x[4] == "1"]: ws.append(r)
for row in ws.iter_rows():
    for c in row: c.font = Font(name="Arial", size=10, bold=(c.row == 1))
wb.save("results/metabolism_22genomes_presence_absence.xlsx")
for r in T3[1:]: print(r)
for name in SETS:
    for c in CL: print(name, "|", c, "|", [r[3] for r in T2[1:] if r[0] == name and r[1] == c])
