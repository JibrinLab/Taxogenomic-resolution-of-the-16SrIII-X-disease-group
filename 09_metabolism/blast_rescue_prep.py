# ANALYSIS: Metabolic gene panel (66 genes, 7 pathways) in the 22 genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Fig. 8, Tables S7.1 to S7.10
# FOLDER: 09_metabolism (see README.md in the parent folder)
import os, sys, glob
sys.path.insert(0, "scripts"); import metab_lib as M
W = os.getcwd(); GT = W + "/gff_tree"; OUTD = W + "/blast_rescue"; os.makedirs(OUTD, exist_ok=True)
RP = "/projects/dsn011/ppss_publication/02_annotation/prokka_refs"
REFQ = ["AT_Ca.P.mali_16SrX", "AYWB_Ca.P.asteris_16SrI", "FD_Ca.P.ulmi_16SrV", "OYM_Ca.P.asteris_16SrI", "AWB_Ca.P.australiense_16SrXII", "PWB_Ca.P.aurantifolia_16SrII", "STOL_Ca.P.solani_16SrXII"]
def readfa(p):
    d, k = {}, None
    for l in open(p):
        l = l.rstrip("\n")
        if l.startswith(">"): k = l[1:].split()[0]; d[k] = [l, []]
        elif k: d[k][1].append(l)
    return d
n = {}
with open(OUTD + "/proteomes22.faa", "w") as out:
    for L in sorted(os.listdir(GT)):
        d = os.path.dirname(os.path.realpath("%s/%s/%s.gff" % (GT, L, L)))
        f = glob.glob(d + "/*.faa")
        if len(f) != 1: print("WARNING", L, d, f); continue
        s = readfa(f[0]); n[L] = len(s)
        for k, (h, q) in s.items(): out.write(">%s|%s\n%s\n" % (L, h[1:], "".join(q)))
print("proteomes:", len(n), n)
nq = {}
with open(OUTD + "/queries.faa", "w") as out:
    for rd in REFQ:
        gff = glob.glob("%s/%s/*.gff" % (RP, rd)); faa = glob.glob("%s/%s/*.faa" % (RP, rd))
        if len(gff) != 1 or len(faa) != 1: print("WARNING ref", rd, gff, faa); continue
        feats = M.parse(gff[0]); s = readfa(faa[0])
        for g in M.order:
            pw, ecs, names, prods = M.G[g]
            for hid, hg, how in M.detect(feats, ecs, names, prods):
                if hid in s: out.write(">%s|%s|%s\n%s\n" % (g, rd, hid, "".join(s[hid][1]))); nq[g] = nq.get(g, 0) + 1
print("genes with no reference query:", [g for g in M.order if g not in nq])
print("queries per gene:", nq)
