# ANALYSIS: Secretome (SignalP, DIAMOND) and effector search in the 22 genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Fig. S6.2, Tables S8.x, S9.x, S10.x
# FOLDER: 10_secretome_effectors (see README.md in the parent folder)
import pandas as pd
GR = {"G1": ["PPS1A","PPS2B","PPS6A","PPS6B","PPS9B","PPS10A","CanS4","MW1","VAC","MA1"], "G2": ["CX","6A1","ZT3-1","PR2021","JR1","PDA15"], "G3": ["CicWB-2022","ChTDIII","OVSP-L01","OVSP-S01","Vc33"], "G4": ["LWB"]}
gof = {g: k for k, v in GR.items() for g in v}
C = pd.read_csv("results/Table_secretome_conservation_Candidates.tsv", sep="\t", dtype=str).fillna("")
hs = C.Genomes_with_hit.apply(lambda s: [x.strip() for x in s.split(",") if x.strip()])
for k in GR: C["Hit_in_" + k] = hs.apply(lambda l: any(gof[x] == k for x in l))
C["QueryGroup"] = C.Genome.map(gof)
t = C.groupby("QueryGroup")[["Hit_in_" + k for k in GR]].agg(["sum"]).droplevel(1, axis=1); t.insert(0, "Candidates", C.groupby("QueryGroup").size())
for k in GR: t["Pct_hit_in_" + k] = (100 * t["Hit_in_" + k] / t.Candidates).round(1)
t.to_csv("results/Table_secretome_conservation_by_group.tsv", sep="\t"); print(t.to_string())
