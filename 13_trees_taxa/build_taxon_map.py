#!/usr/bin/env python3
# ANALYSIS: Trees with the four species marked
# MANUSCRIPT SECTION: Figures
# USED FOR: Tree figures
# FOLDER: 13_trees_taxa (see README.md in the parent folder)
"""Species (taxon) of every tip, from the Curated83 sheet.
usage: build_taxon_map.py [--pick PICK.tsv] [--groups GROUPS.tsv] CURATED_XLSX GENOME_LABELS_TSV OUT_MAP.tsv
One name per cluster. When NCBI gives several species names in one cluster, the name carried by most genomes of the
cluster is used, unless PICK.tsv (cluster<TAB>name, no header) says otherwise.
Clusters with no species name are called 'Phytoplasma sp. (16SrX)' from GROUPS.tsv (genome label<TAB>group, no header)
or, if the sheet has a column with the 16Sr group, from that column; clusters without a group are listed at the end.
Each genome belongs to an ANIm 95% cluster (column ANIm_cluster95). A cluster is named
  - by the proposed species when most of its members are in our four groups,
  - else by the NCBI organism name when the cluster members carry a 'Candidatus Phytoplasma <epithet>' name,
  - else 'Phytoplasma sp. (ANI cluster N)'.
Writes OUT_MAP.tsv (tip, taxon) and OUT_MAP.review.tsv (every cluster, its members and NCBI names) for checking."""
import sys, re
import pandas as pd
A = sys.argv[1:]; pick = {}; grp = {}
while A and A[0] in ("--pick", "--groups"):
    d = pick if A[0] == "--pick" else grp
    for l in open(A[1]):
        p = l.rstrip("\n").split("\t")
        if len(p) >= 2 and p[0] and not p[0].lower().startswith(("cluster", "genome", "tip")): d[p[0].strip()] = p[1].strip()
    A = A[2:]
sys.argv = [sys.argv[0]] + A
GROUPS = {"Ca. P. caryae": ["PPS1A","PPS2B","PPS6A","PPS6B","PPS9B","PPS10A","CanS4","MW1","VAC","MA1"],
          "Ca. P. pruni": ["CX","6A1","ZT3-1","PR2021","JR1","PDA15"],
          "Ca. P. manihotii": ["CicWB-2022","Vc33","ChTDIII","OVSP-L01","OVSP-S01"],
          "Ca. P. lactucae": ["LWB"]}
key = lambda s: re.sub(r"[^A-Za-z0-9]", "", str(s)).upper()
OURS = {key(n): t for t, m in GROUPS.items() for n in m}
def short(org):
    m = re.search(r"Candidatus Phytoplasma ([A-Za-z\-]+)", str(org))
    if not m or m.group(1).rstrip("'").lower() in ("sp", "spp"): return None
    return "Ca. P. " + m.group(1).rstrip("'")
x = pd.read_excel(sys.argv[1], sheet_name="Curated83")
lab = dict(l.rstrip("\n").split("\t")[:2] for i, l in enumerate(open(sys.argv[2])) if i)
x["label"] = x["Genome"].map(lambda g: lab.get(g, re.sub(r"^GC[AF]_\d+\.\d_", "", g)))
x["ours"] = x["label"].map(lambda l: OURS.get(key(l)))
x["cand"] = x["Organism name"].map(short)
gcol = [c for c in x.columns if re.search(r"16\s*sr|group", str(c), re.I) and "cluster" not in str(c).lower()]
print("columns in sheet:", list(x.columns)); print("possible 16Sr group columns:", gcol)
need = []; ties = []
names = {}; rev = []
for c, d in x.groupby("ANIm_cluster95"):
    o = d["ours"].dropna()
    if len(o) and len(o) >= len(d) / 2: nm = o.value_counts().idxmax()
    else:
        vc = d["cand"].dropna().value_counts()
        if str(c) in pick: nm = pick[str(c)]
        elif len(vc):
            tied = sorted(vc[vc == vc.max()].index); nm = tied[0]
            if len(tied) > 1: ties.append((c, nm, tied))
        else:
            g = [grp[l] for l in d["label"] if l in grp]
            for col in gcol:
                g += [str(v) for v in d[col].dropna() if str(v).strip() and str(v).lower() != "nan"]
            if g: nm = "Phytoplasma sp. (%s)" % pd.Series(g).value_counts().idxmax()
            else: nm = "Phytoplasma sp. (strain %s)" % ", ".join(d["label"]); need.append((c, ", ".join(d["label"])))
    names[c] = nm
    rev.append((c, nm, len(d), ", ".join(d["label"]), "; ".join(sorted(set(d["Organism name"].astype(str))))))
with open(sys.argv[3], "w") as f:
    f.write("tip\ttaxon\n")
    for _, r in x.iterrows():
        f.write("%s\t%s\n" % (r["label"], names[r["ANIm_cluster95"]]))
        f.write("%s\t%s\n" % (r["Genome"], names[r["ANIm_cluster95"]]))
    f.write("A_laidlawii\tAcholeplasma laidlawii\nPG-8A\tAcholeplasma laidlawii\n")
    for t, m in GROUPS.items():
        for n in m: f.write("%s\t%s\n" % (n, t))
pd.DataFrame(rev, columns=["ANIm_cluster95", "taxon_assigned", "n_genomes", "genomes", "NCBI_organism_names"]).to_csv(sys.argv[3].replace(".tsv", ".review.tsv"), sep="\t", index=False)
multi = [(r[0], r[1]) for r in rev if len(set(re.findall(r"Candidatus Phytoplasma ([A-Za-z\-]+)", r[4]))) > 1]
print("\nclusters with several NCBI names (name used shown):"); [print(" ", m) for m in multi]
print("\nTIES (names carried by the same number of genomes; first alphabetically was used, needs a decision):"); [print(" ", t) for t in ties]
print("\nclusters still needing a 16Sr group (add to GROUPS.tsv):"); [print(" ", n) for n in need]
print(pd.DataFrame(rev, columns=["cluster", "taxon", "n", "genomes", "ncbi"]).drop(columns=["ncbi"]).to_string(max_colwidth=70))
