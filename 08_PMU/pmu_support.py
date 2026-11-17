# ANALYSIS: Potential mobile unit (PMU) genes in the 22 genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Fig. 7, Tables S6.1 to S6.5
# FOLDER: 08_PMU (see README.md in the parent folder)
import csv, re, subprocess
W = "results/pmu_window_genes.tsv"; GFF = "../05_metabolism_22genomes/gff_tree/"; FA = "../06_secretome_effector_22genomes/blast/all_proteomes.faa"
prot = {}; k = None
for l in open(FA):
    if l[0] == ">": k = l[1:].split()[0]; prot[k] = []
    else: prot[k].append(l.strip())
prot = {a: "".join(b) for a, b in prot.items()}
def cds_id(lab, s, e, st):
    best = (0, None)
    for l in open(f"{GFF}{lab}/{lab}.gff"):
        if l.startswith("##FASTA"): break
        c = l.rstrip("\n").split("\t")
        if len(c) > 8 and c[2] == "CDS" and c[6] == st:
            a, b = int(c[3]), int(c[4]); ov = min(b, e) - max(a, s) + 1; j = ov / (max(b, e) - min(a, s) + 1)
            if j > best[0] and j >= 0.7:
                m = re.search(r"ID=([^;]+)", c[8]); best = (j, m.group(1) if m else None)
    return best[1]
rows = list(csv.DictReader(open(W), delimiter="\t")); seqs = {}
for r in rows:
    lab = r["Genome"].replace("OVSP_", "OVSP-"); i = cds_id(lab, int(r["Start"]), int(r["End"]), r["Strand"])
    if i and f"{lab}__{i}" in prot: seqs[f"{lab}|{r['Gene']}"] = prot[f"{lab}__{i}"]
open("support_q.faa", "w").write("".join(f">{a}\n{b}\n" for a, b in seqs.items()))
open("support_ref.faa", "w").write("".join(f">{a}\n{b}\n" for a, b in seqs.items() if a.startswith("PPS1A|")))
subprocess.run("makeblastdb -in support_ref.faa -dbtype prot -out support_ref > /dev/null", shell=True, check=True)
o = subprocess.run('blastp -query support_q.faa -db support_ref -evalue 1e-3 -seg no -comp_based_stats 0 -outfmt "6 qseqid sseqid pident length qlen slen bitscore"', shell=True, capture_output=True, text=True).stdout
best = {}
for l in o.splitlines():
    q, s, pid, n, ql, sl, b = l.split("\t"); g = q.split("|")[1]
    if s == f"PPS1A|{g}" and (q not in best or float(b) > best[q][5]): best[q] = (float(pid), int(n), int(ql), int(sl), g, float(b))
out = open("results/pmu_window_gene_support.tsv", "w"); out.write("Genome\tGene\tProtein_aa\tPPS1A_protein_aa\tIdentity_to_PPS1A_pct\tAlignment_aa\n")
for r in rows:
    lab = r["Genome"].replace("OVSP_", "OVSP-"); q = f"{lab}|{r['Gene']}"
    if q not in seqs: out.write(f"{lab}\t{r['Gene']}\t\t\tno annotated CDS (TBLASTN hit)\t\n"); continue
    h = best.get(q); ref = len(seqs.get(f"PPS1A|{r['Gene']}", ""))
    out.write(f"{lab}\t{r['Gene']}\t{len(seqs[q])}\t{ref}\t{h[0] if h else 'no hit'}\t{h[1] if h else ''}\n")
out.close(); print(open("results/pmu_window_gene_support.tsv").read())
