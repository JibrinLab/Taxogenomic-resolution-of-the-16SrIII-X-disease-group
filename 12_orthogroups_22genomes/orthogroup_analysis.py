# ANALYSIS: OrthoFinder orthogroups of the 22 genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Figs. S6.3 to S6.5
# FOLDER: 12_orthogroups_22genomes (see README.md in the parent folder)
import glob, pandas as pd
from openpyxl.styles import Font
B = "/projects/dsn011/ppss_publication/2026-10-27_22genome_reanalysis/"
EO = pd.read_csv(B + "06_secretome_effector_22genomes/results/effector_matrix.tsv", sep="\t").Strain.tolist(); PPS6, REF = EO[:6], EO[6:]
prod = {}
for f in glob.glob("proteins22/*.faa"):
    for l in open(f):
        if l.startswith(">"): k, _, v = l[1:].rstrip().partition(" "); prod[k] = v
def analyse(run, PPS, REF):
    d = pd.read_csv(f"{run}/Orthogroups/Orthogroups.tsv", sep="\t", dtype=str).fillna("").set_index("Orthogroup"); G = PPS + REF; P = d[G] != ""
    nP, nR = P[PPS].sum(axis=1), P[REF].sum(axis=1); cat = pd.Series("other", index=d.index)
    cat[(nP > 0) & (nR > 0)] = "shared (PPS and reference)"; cat[(nP == len(PPS)) & (nR == len(REF))] = "present in all genomes"
    cat[(nP > 0) & (nR == 0)] = "PPS only"; cat[(nP == 0) & (nR > 0)] = "reference only"
    return d, P, nP, nR, cat
def tables(run, PPS, REF):
    d, P, nP, nR, cat = analyse(run, PPS, REF); G = PPS + REF; s = cat.value_counts()
    summ = pd.DataFrame({"Category": s.index, "Orthogroups": s.values})
    extra = [("PPS only, present in all PPS genomes", int(((cat == "PPS only") & (nP == len(PPS))).sum())), ("reference only, present in all references", int(((cat == "reference only") & (nR == len(REF))).sum()))]
    summ = pd.concat([summ, pd.DataFrame(extra, columns=["Category", "Orthogroups"])], ignore_index=True)
    pg = pd.DataFrame([(g, c, int(((cat == c) & P[g]).sum())) for g in G for c in ["present in all genomes", "shared (PPS and reference)", "PPS only", "reference only"]], columns=["Genome", "Category", "Orthogroups"])
    def ex(c, S):
        r = []
        for og in cat.index[cat == c]:
            have = [g for g in S if d.at[og, g]]; first = d.at[og, have[0]].split(", ")[0]
            r.append((og, len(have), "; ".join(f"{g} {len(d.at[og, g].split(', '))}" for g in have), prod.get(first, "")))
        return pd.DataFrame(r, columns=["Orthogroup", "Genomes_present", "Genes_per_genome", "Product_of_first_gene"]).sort_values("Genomes_present", ascending=False)
    return summ, pg, ex("PPS only", PPS), ex("reference only", REF)
R22 = tables("proteins22/OrthoFinder/Results_r22", PPS6, REF); R4 = tables("proteins20/OrthoFinder/Results_r4pps", PPS6[:4], REF)
leg = pd.DataFrame([("Run", "OrthoFinder 3.1.5 (-S diamond -M msa) on the 22 proteomes (main sheets) and on the 20 proteomes of the four complete PPS genomes plus the 16 references (Four_PPS sheets). Orthogroup IDs are not comparable between the two runs."), ("Categories", "present in all genomes; shared (at least one PPS and one reference genome); PPS only (at least one PPS genome, no reference); reference only (at least one reference, no PPS). Present means at least one gene in the orthogroup."), ("Products", "Prokka product of the first gene of the orthogroup in the first genome that has it.")], columns=["Item", "Definition"])
sh = {"Summary": R22[0], "Per_genome": R22[1], "PPS_only": R22[2], "Reference_only": R22[3], "Four_PPS_summary": R4[0], "Four_PPS_per_genome": R4[1], "Four_PPS_only": R4[2], "Four_PPS_reference_only": R4[3], "Legend": leg}
REFA = [r.replace("-", "_") for r in REF if r != "Vc33"]
RA = tables("/projects/dsn011/ppss_publication/2026-10-20_new_species_taxonomy/06_core_genome_tree/proteins/OrthoFinder/Results_Oct20", PPS6, REFA)
RA = list(RA); RA[0] = RA[0][RA[0].iloc[:, 0] != "other"].reset_index(drop=True)
leg.loc[len(leg)] = ["Curated_84 sheets", "Earlier OrthoFinder 3.1.5 run on 84 curated proteomes (-S diamond -M msa), summarised for the six PPS genomes against the 15 reference genomes it contains. Products are available only for genes whose IDs match the 22-genome proteomes."]
sh.update({"Curated84_summary": RA[0], "Curated84_per_genome": RA[1], "Curated84_PPS_only": RA[2], "Curated84_reference_only": RA[3]})
print(RA[0].to_string(index=False))
with pd.ExcelWriter("Table_orthogroups.xlsx", engine="openpyxl") as w:
    for n, d in sh.items():
        d.to_excel(w, sheet_name=n, index=False); ws = w.sheets[n]
        for row in ws.iter_rows():
            for c in row: c.font = Font(name="Arial", color="000000", bold=(c.row == 1))
        for i in range(1, ws.max_column + 1): ws.column_dimensions[ws.cell(1, i).column_letter].width = 80 if n == "Legend" else 26
        if n != "Legend": d.to_csv(f"Table_orthogroups_{n}.tsv", sep="\t", index=False)
print(R22[0].to_string(index=False)); print(R4[0].to_string(index=False))
