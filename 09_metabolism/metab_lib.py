# ANALYSIS: Metabolic gene panel (66 genes, 7 pathways) in the 22 genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Fig. 8, Tables S7.1 to S7.10
# FOLDER: 09_metabolism (see README.md in the parent folder)
import os, re, csv
P = "/projects/dsn011/ppss_publication/02_annotation/prokka"
R = "/projects/dsn011/ppss_publication/02_annotation/prokka_refs"
OUT = "/projects/dsn011/ppss_publication/2026-10-27_22genome_reanalysis/05_metabolism_22genomes/results"
GENEFILE = "/projects/dsn011/ppss_publication/2026-10-27_22genome_reanalysis/05_metabolism_22genomes/scripts/metab_genes.txt"

PPS = ["PPS1A","PPS2B","PPS6A","PPS6B","PPS9B","PPS10A"]
REFS = ["CanS4_Ca.P.pruni_16SrIII","CX_Ca.P.pruni_16SrIII","6A1_Ca.P.pruni_16SrIII",
        "PR2021_Ca.P.pruni_16SrIII","ZT3-1_Ca.P.pruni_16SrIII","ChTDIII_Ca.P.pruni_16SrIII"]

G = {}
order = []
with open(GENEFILE) as fh:
    for line in fh:
        line = line.rstrip("\n")
        if not line:
            continue
        parts = line.split("|")
        gene, pathway, ecs, names, prods = parts
        ecs = [e for e in ecs.split(",") if e]
        names = [n for n in names.split(",") if n]
        prods = [p for p in prods.split(",") if p]
        G[gene] = (pathway, ecs, names, prods)
        order.append(gene)

AR = re.compile(r'(\w+)=([^;]*)')

def parse(path):
    feats = []
    if not os.path.exists(path):
        return feats
    for line in open(path):
        if line.startswith("#"):
            if line.startswith("##FASTA"):
                break
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 9 or c[2] != "CDS":
            continue
        a = dict(AR.findall(c[8]))
        feats.append({
            "id": a.get("ID",""),
            "gene": (a.get("gene") or a.get("Name") or "").strip(),
            "ec": (a.get("eC_number") or "").strip(),
            "product": (a.get("product") or "").strip(),
        })
    return feats

def base(g):
    return re.sub(r'_\d+$', '', g)

def detect(feats, ecs, names, prods):
    hits = []
    ln = [n.lower() for n in names]
    for f in feats:
        how = None
        if ecs and f["ec"]:
            for e in ecs:
                if f["ec"] == e:
                    how = "EC:" + f["ec"]
                    break
        if not how and f["gene"] and base(f["gene"]).lower() in ln:
            how = "name:" + f["gene"]
        if not how and f["product"]:
            pl = f["product"].lower()
            for kw in prods:
                if kw and kw in pl:
                    how = "product:" + f["product"][:45]
                    break
        if how:
            hits.append((f["id"], f["gene"], how))
    return hits

def run(labs, base_dir, tag):
    os.makedirs(OUT, exist_ok=True)
    long_rows, evid = [], []
    matrix = {g: {} for g in order}
    for L in labs:
        feats = parse(os.path.join(base_dir, L, L + ".gff"))
        if not feats:
            print("  WARNING no CDS parsed: " + L)
        for g in order:
            pw, ecs, names, prods = G[g]
            hits = detect(feats, ecs, names, prods)
            matrix[g][L] = len(hits)
            long_rows.append([L, pw, g, 1 if hits else 0, len(hits)])
            for hid, hg, how in hits:
                evid.append([L, pw, g, hid, hg, how])
        print("  %s: %d CDS, %d/%d genes detected" % (
            L, len(feats), sum(1 for g in order if matrix[g][L]), len(order)))

    w = csv.writer(open("%s/metabolic_%s_long.tsv" % (OUT, tag), "w"), delimiter="\t")
    w.writerow(["Sample","Pathway","Gene","Present","CopyNumber"])
    w.writerows(long_rows)

    w = csv.writer(open("%s/metabolic_%s_matrix.tsv" % (OUT, tag), "w"), delimiter="\t")
    w.writerow(["Pathway","Gene"] + labs)
    for g in order:
        w.writerow([G[g][0], g] + [matrix[g][l] for l in labs])

    w = csv.writer(open("%s/metabolic_%s_evidence.tsv" % (OUT, tag), "w"), delimiter="\t")
    w.writerow(["Sample","Pathway","Gene","LocusTag","ProkkaName","MatchedBy"])
    w.writerows(evid)
    return matrix

