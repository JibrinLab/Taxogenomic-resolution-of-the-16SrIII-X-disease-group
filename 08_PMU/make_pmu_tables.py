# ANALYSIS: Potential mobile unit (PMU) genes in the 22 genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Fig. 7, Tables S6.1 to S6.5
# FOLDER: 08_PMU (see README.md in the parent folder)
import pandas as pd
R = "results/"
grp = {**{k: "Group 1" for k in "PPS1A PPS2B PPS6A PPS6B PPS9B PPS10A CanS4 MW1 VAC MA1".split()},
       **{k: "Group 2" for k in "CX 6A1 ZT3-1 PR2021 JR1 PDA15".split()},
       **{k: "Group 3" for k in "CicWB-2022 ChTDIII OVSP-L01 OVSP-S01 Vc33".split()}, "Group 4": "LWB"}
fx = lambda s: str(s).replace("OVSP_", "OVSP-")
h = pd.read_csv(R + "pmu_hybrid_cluster.tsv", sep="\t").fillna("")
h["Genome"] = h["Genome"].map(fx); h.insert(1, "Group", h["Genome"].map(grp))
h["GenesFoundOf5"] = (h[[c for c in h.columns if c.endswith("_loci")]] > 0).sum(axis=1)
w = pd.read_csv(R + "pmu_window_genes.tsv", sep="\t"); w["Genome"] = w["Genome"].map(fx)
hl = pd.read_csv(R + "pmu_hits_long.tsv", sep="\t").drop(columns=["Group"], errors="ignore"); hl["Genome"] = hl["Genome"].map(fx)
cn = pd.read_csv(R + "pmu_copy_numbers.tsv", sep="\t")
if "Genome" in cn.columns: cn["Genome"] = cn["Genome"].map(fx)
leg = pd.DataFrame({"Sheet": ["Window_summary", "Window_genes", "Prokka_name_hits", "Copy_numbers"], "Content": [
 "One row per genome. Window = smallest region on one contig holding at least one copy of sigF, ssb, hupB, ftsH and smc. Span_bp = end of last gene minus start of first. *_loci = candidate loci per gene. GenesFoundOf5 = genes with at least one locus.",
 "Contig, start, end and strand of the five genes inside each window (contig coordinates).",
 "Hits by Prokka gene name only (ftsH limited to ftsH2); use Window_summary for calls.",
 "Copy numbers per genome from the earlier name-based detection."]})
sheets = {"Legend": leg, "Window_summary": h, "Window_genes": w, "Prokka_name_hits": hl, "Copy_numbers": cn}
sheets["Window_gene_support"] = pd.read_csv(R + "pmu_window_gene_support.tsv", sep="\t")
leg.loc[len(leg)] = ["Window_gene_support", "BLASTp identity of each window protein to the same-named PPS1A window protein (BLAST+ 2.16.0, e-value 1e-3). Blank identity with a note = no annotated CDS at that position (TBLASTN hit only)."]
with pd.ExcelWriter(R + "Supplementary_Table_PMU_22genomes.xlsx", engine="openpyxl") as x:
    for n, d in sheets.items():
        d.to_excel(x, sheet_name=n, index=False)
        ws = x.sheets[n]
        for col in ws.columns:
            ws.column_dimensions[col[0].column_letter].width = min(60, max(10, max(len(str(c.value or "")) for c in col[:30]) + 2))
        for c in ws[1]: c.font = c.font.copy(bold=True)
for n, d in sheets.items():
    if n != "Legend": d.to_csv(R + "sheet_" + n + ".tsv", sep="\t", index=False)
h.to_csv(R + "pmu_window_summary_22genomes.tsv", sep="\t", index=False)
print(h[["Genome", "Group", "AllFiveOnOneContig", "MinSpan_bp", "Order", "GenesFoundOf5"]].to_string(index=False))
