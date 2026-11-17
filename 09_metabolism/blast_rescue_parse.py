# ANALYSIS: Metabolic gene panel (66 genes, 7 pathways) in the 22 genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Fig. 8, Tables S7.1 to S7.10
# FOLDER: 09_metabolism (see README.md in the parent folder)
import csv, collections
rd = csv.reader(open("results/metabolic_all22_matrix.tsv"), delimiter="\t"); h = next(rd); labs = h[2:]
pw, gffp, order = {}, {}, []
for r in rd:
    pw[r[1]] = r[0]; order.append(r[1])
    for L, x in zip(labs, r[2:]): gffp[(r[1], L)] = float(x) > 0
ev = collections.defaultdict(set)
for r in csv.DictReader(open("results/metabolic_all22_evidence.tsv"), delimiter="\t"): ev[(r["Sample"], r["LocusTag"])].add(r["Gene"])
best = {}
for l in open("blast_rescue/blastp_all.tsv"):
    f = l.rstrip("\n").split("\t"); g, ref, _ = f[0].split("|", 2); L, loc = f[1].split("|", 1)
    pid, ql, qs, qe, e, bs = float(f[2]), int(f[4]), int(f[5]), int(f[6]), float(f[7]), float(f[8])
    prod = f[9].split(" ", 1)[1] if " " in f[9] else ""
    k = (g, L)
    if k not in best or (e, -bs) < (best[k][0], -best[k][1]): best[k] = (e, bs, ref, pid, (qe - qs + 1) / ql, loc, prod)
rows, fin = [], {}
for g in order:
    for L in labs:
        p = gffp[(g, L)]; b = best.get((g, L))
        oth = sorted(ev.get((L, b[5]), set()) - {g}) if b else []
        q = "" if not b else ("strong" if b[3] >= 30 and b[4] >= 0.5 else "weak")
        hyp = bool(b) and "hypothetical" in b[6].lower() and not oth
        ok = bool(p or hyp); fin[(g, L)] = int(ok)
        basis = "Prokka" if p else ("BLASTp only" if hyp else ("BLASTp hit is a protein assigned to another gene" if oth else ("BLASTp hit to a protein with a different annotated function (not counted)" if b else "absent")))
        rows.append([L, pw[g], g, int(p), int(bool(b)), b[2] if b else "", "%.1f" % b[3] if b else "", "%.2f" % b[4] if b else "", "%.2g" % b[0] if b else "", b[5] if b else "", b[6] if b else "", q, ",".join(oth), int(ok), basis])
w = csv.writer(open("results/blastp_rescue_all.tsv", "w"), delimiter="\t")
w.writerow(["Genome", "Pathway", "Gene", "Prokka", "BLASTp", "Query", "Pident", "QCov", "Evalue", "SubjectLocus", "SubjectProduct", "Quality", "AssignedToOtherGene", "Final", "Basis"]); w.writerows(rows)
w = csv.writer(open("results/metabolic_final22_matrix.tsv", "w"), delimiter="\t")
w.writerow(["Pathway", "Gene"] + labs)
for g in order: w.writerow([pw[g], g] + [fin[(g, L)] for L in labs])
print("gene | Prokka n | rescued strong | rescued weak | hit belongs to other gene | final n")
for g in order:
    rr = [x for x in rows if x[2] == g]; s = sum(x[3] for x in rr)
    st = sum(1 for x in rr if x[14] == "BLASTp only" and x[11] == "strong"); wk = sum(1 for x in rr if x[14] == "BLASTp only" and x[11] == "weak")
    ex = sum(1 for x in rr if x[14].startswith("BLASTp hit "))
    if st or wk or ex: print(g, s, st, wk, ex, s + st + wk)
for x in rows:
    if x[14].startswith("BLASTp hit is") and x[0] == "PPS1A": print("OTHER", x[2], "->", x[12], "|", x[10])
