# ANALYSIS: Marker genes, orthology and phylogenetic trees
# MANUSCRIPT SECTION: Multilocus sequence comparison; Phylogenetic reconstructions
# USED FOR: Fig. S2.8, Figs. S4.1 to S4.6
# FOLDER: 04_markers_orthology_trees (see README.md in the parent folder)
import csv,collections
R=["PPS1A","OVSP_L01","CicWB_2022","LWB","CX"]
D=collections.defaultdict(dict); O=collections.defaultdict(lambda:collections.defaultdict(float)); sites={}
for r in csv.DictReader(open("marker_identity_long.tsv"),delimiter="\t"):
    if r["identity"]=="": continue
    m,a,g,gr,i=r["marker"],r["reference"],r["genome"],r["group"],float(r["identity"])
    if g in R: D[m][(a,g)]=i; sites[m]=int(r["aligned_sites"])
    if gr=="other": O[m][a]=max(O[m][a],i)
w=csv.writer(open("marker_identity_summary.tsv","w"),delimiter="\t")
w.writerow(["marker","sites_PPS1A_vs_OVSPL01","lowest_cross_reference_identity","highest_identity_to_non_study_species","gap","reference_with_highest_non_study_identity"])
for m in sorted(D):
    pr=[v for (a,g),v in D[m].items()]
    lo=min(pr); hi=max(O[m].values()) if O[m] else None; ref=max(O[m],key=O[m].get) if O[m] else ""
    w.writerow([m,D[m].get(("PPS1A","OVSP_L01"),""),lo,hi,round(lo-hi,2) if hi else "",ref])
for l in open("marker_identity_summary.tsv"): print(l.rstrip())
