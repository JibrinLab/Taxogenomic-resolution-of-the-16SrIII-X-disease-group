# ANALYSIS: Secretome (SignalP, DIAMOND) and effector search in the 22 genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Fig. S6.2, Tables S8.x, S9.x, S10.x
# FOLDER: 10_secretome_effectors (see README.md in the parent folder)
import csv, collections
W = "results/"
Q = {"SAP05":"SAP05_AYWB","SAP06/48-like":"SAP06like_Kai573","SAP11":"SAP11_AYWB","SAP54":"SAP54_AYWB","SAP54-like":"SAP54like_Kai704","TENGU":"TENGU_OYM","SVM_705":"SVM_Kai705","SVM_1019":"SVM_Kai1019","SVM_1100":"SVM_Kai1100","SVM_1102":"SVM_Kai1102","SVM_1103":"SVM_Kai1103"}
ORDER = "PPS1A PPS2B PPS6A PPS6B PPS9B PPS10A CanS4 MW1 VAC MA1 CX 6A1 ZT3-1 PR2021 JR1 PDA15 CicWB-2022 Vc33 ChTDIII OVSP-L01 OVSP-S01 LWB".split()
def sp(s): return s.replace("lcl|", "").split("__", 1)
prot = {}
for r in csv.reader(open("blast/blastp.tsv"), delimiter="\t"):
    g, p = sp(r[1]); pid = float(r[2]); cov = 100*(int(r[7])-int(r[6])+1)/int(r[4])
    if (r[0], g) not in prot or cov*pid > prot[(r[0], g)][3]: prot[(r[0], g)] = (p, pid, cov, cov*pid)
H = collections.defaultdict(list)
for fn in ("blast/tblastn.tsv", "blast/tblastn_tengu.tsv"):
    for r in csv.reader(open(fn), delimiter="\t"):
        if fn.endswith("tengu.tsv") != r[0].startswith("TENGU"): continue
        g, ct = sp(r[1]); ss, se = int(r[8]), int(r[9])
        H[(r[0], g, ct, "+" if ss < se else "-")].append(dict(pid=float(r[2]), n=int(r[3]), qlen=int(r[4]), slen=int(r[5]), qs=int(r[6]), qe=int(r[7]), lo=min(ss, se), hi=max(ss, se), fr=r[10], star=r[13].count("*")))
def joins(a, h, st):
    dq = (h["qs"]-a["qe"]) if st == "+" else (a["qs"]-h["qe"])
    return -15 <= h["lo"]-a["hi"] <= 150 and -10 <= dq <= 15
loci = collections.defaultdict(list)
for (q, g, ct, st), L in H.items():
    L.sort(key=lambda x: (x["lo"], x["hi"])); groups = []
    for h in L:
        for grp in groups:
            if joins(grp[-1], h, st): grp.append(h); break
        else: groups.append([h])
    for c in groups:
        cov = set()
        for h in c: cov.update(range(h["qs"], h["qe"]+1))
        tot = sum(h["n"] for h in c); sg = sum(1 for a, b in zip(c, c[1:]) if 0 < b["lo"]-a["hi"] <= 6 and a["fr"] == b["fr"])
        loci[(q, g)].append(dict(ct=ct, st=st, lo=min(h["lo"] for h in c), hi=max(h["hi"] for h in c), qlen=c[0]["qlen"], slen=c[0]["slen"],
            pid=sum(h["pid"]*h["n"] for h in c)/tot, qcov=100*len(cov)/c[0]["qlen"], frames=len({h["fr"] for h in c}),
            stars=sum(h["star"] for h in c)+sg, qs=min(h["qs"] for h in c), qe=max(h["qe"] for h in c)))
def broken(o):
    mN, mC = o["qs"]-1, o["qlen"]-o["qe"]; L, R = o["lo"]-1, o["slen"]-o["hi"]
    dN, dC = (L, R) if o["st"] == "+" else (R, L)
    return (mN >= 3 and dN < 3*mN+30) or (mC >= 3 and dC < 3*mC+30)
def lcall(o):
    if o["pid"] >= 60:
        if o["stars"] > 0 and o["qcov"] >= 50: return "stop codon"
        if o["qcov"] >= 80: return "frameshift" if o["frames"] > 1 else "present"
        return "contig-break" if broken(o) else "partial"
    return "paralog" if o["pid"] >= 35 and o["qcov"] >= 70 else None
RK = {"present":0, "frameshift":1, "stop codon":2, "contig-break":3, "partial":4, "paralog":5}
def rank(c, o): return 6 if c == "partial" and o["qcov"] < 50 else RK[c]
def call(loc):
    cl = [(lcall(o), o) for o in loc]; cl = [(c, o) for c, o in cl if c]
    if not cl: return "absent", None, ""
    cl.sort(key=lambda x: (rank(*x), -x[1]["pid"], -x[1]["qcov"])); c, b = cl[0]
    return c, b, ", ".join(f"{n} {k}" for k, n in collections.Counter(x[0] for x in cl[1:]).items())
code = {"present":"P", "absent":"-", "paralog":"para", "contig-break":"CB", "stop codon":"STOP", "frameshift":"FS", "partial":"part"}
lg = open(W+"effector_calls_long.tsv", "w")
lg.write("Genome\tEffector\tQuery\tCall\tIdentity\tQueryCoverage\tContig\tStart\tEnd\tStrand\tAdditionalLoci\tProteinID\tProteinIdentity\tProteinCoverage\n")
M, D = {}, {}
for g in ORDER:
    for e, q in Q.items():
        c, b, ex = call(loci.get((q, g), [])); p = prot.get((q, g)); M[(g, e)] = c
        D[(g, e)] = c + (f" ({b['pid']:.0f}%" + (f"; +{ex}" if ex else "") + ")" if b else "")
        lg.write("\t".join(str(x) for x in (g, e, q, c, f"{b['pid']:.1f}" if b else "", f"{b['qcov']:.0f}" if b else "", b["ct"] if b else "", b["lo"] if b else "", b["hi"] if b else "", b["st"] if b else "", ex, p[0] if p else "", f"{p[1]:.1f}" if p else "", f"{p[2]:.0f}" if p else "")) + "\n")
lg.close()
for fn, X in (("effector_matrix.tsv", M), ("effector_matrix_detail.tsv", D)):
    with open(W+fn, "w") as m:
        m.write("Strain\t" + "\t".join(Q) + "\n")
        for g in ORDER: m.write(g + "\t" + "\t".join(X[(g, e)] for e in Q) + "\n")
print("Strain".ljust(11) + "".join(e[:8].ljust(9) for e in Q))
for g in ORDER: print(g.ljust(11) + "".join((code[M[(g, e)]] + ("+" if "+" in D[(g, e)] else "")).ljust(9) for e in Q))
