# ANALYSIS: 16S rRNA gene analysis
# MANUSCRIPT SECTION: Multilocus sequence comparison and 16S classification
# USED FOR: Figs. S2.1 to S2.7, Tables S2.1 to S2.6
# FOLDER: 03_16S (see README.md in the parent folder)
import re,time,json,urllib.request,urllib.parse
from Bio.Align import PairwiseAligner
D="/projects/dsn011/ppss_publication/2026-10-20_new_species_taxonomy"
def readfa(p):
    s={};h=None
    for l in open(p):
        l=l.strip()
        if l.startswith(">"): h=l[1:]; s[h]=[]
        elif h: s[h].append(l.upper())
    return {k:"".join(v) for k,v in s.items()}
aln=readfa(D+"/04_16S/16S_allcopies.aln.fasta")
raw={k:v.replace("-","") for k,v in aln.items()}
lab=lambda h:h.rsplit("_c",1)[0]
def colid(a,b):
    n=m=0
    for x,y in zip(a,b):
        if x in "ACGT" and y in "ACGT": n+=1; m+=x!=y
    return (100*(1-m/n),n) if n else (0.0,0)
al=PairwiseAligner(); al.mode="global"; al.match_score=2; al.mismatch_score=-3
al.open_gap_score=-6; al.extend_gap_score=-1; al.target_end_gap_score=0; al.query_end_gap_score=0
rc=lambda s:s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
def pid(a,b):
    best=None
    for strand,bb in (("+",b),("-",rc(b))):
        x=al.align(a,bb)[0]
        if best is None or x.score>best[0]: best=(x.score,x,bb,strand)
    _,x,bb,strand=best
    n=m=0
    for (ts,te),(qs,qe) in zip(*x.aligned):
        for p,q in zip(a[ts:te],bb[qs:qe]):
            if p in "ACGT" and q in "ACGT": n+=1; m+=p!=q
    return (100*(1-m/n) if n else 0.0),n,strand
print("== validation: PPS1A_c1 vs PR2021_c2 (alignment gave 99.08% over 1516 sites)")
i,n,s=pid(raw["PPS1A_c1"],raw["PR2021_c2"]); print("   pairwise: %.2f%% over %d sites, strand %s"%(i,n,s))
print("\n== group 3 copies vs all copies of the six pruni genomes (existing alignment)")
pr=[h for h in aln if lab(h) in ("CX","6A1","ZT3_1","PR2021","JR1","PDA15")]
print("   pruni copies:",len(pr))
for g in ("ChTDIII","CicWB_2022","Vc33","OVSP_L01","OVSP_S01"):
    for h in sorted(k for k in aln if lab(k)==g):
        v=[colid(aln[h],aln[p])[0] for p in pr]
        i1,n1=colid(aln["PPS1A_c1"],aln[h])
        print("   %-16s min %.2f  max %.2f vs pruni | vs PPS1A_c1 %.2f (%d sites)%s"%(h,min(v),max(v),i1,n1,"   <98.65 vs pruni" if min(v)<98.65 else ""))
print("\n== five queries vs the fetched described-species sequences (pairwise)")
new=readfa(D+"/04_16S/described_missing_16S.fasta")
out=open(D+"/04_16S/identity_vs_described_species_extra.tsv","w"); out.write("query\tspecies\taccession\tidentity_pct\tsites\tstrand\tlength\n")
for q in ("PPS1A_c1","CX_c1","OVSP_L01_c1","CicWB_2022_c1","LWB_c1"):
    res=[]
    for h,sq in new.items():
        i,n,s=pid(raw[q],sq); sp,ac=h.split("|"); res.append((i,n,s,sp,ac,len(sq)))
        out.write("%s\t%s\t%s\t%.2f\t%d\t%s\t%d\n"%(lab(q),sp,ac,i,n,s,len(sq)))
    res.sort(reverse=True)
    print("%s: >=98.65%%: %d of %d"%(q,sum(r[0]>=98.65 for r in res),len(res)))
    for i,n,s,sp,ac,L in res[:5]: print("   %-16s %-11s %.2f%% (%d sites, strand %s, len %d)%s"%(sp,ac,i,n,s,L,"  CHECK" if n<1000 or s=="-" else ""))
out.close()
print("\n== candidates for the 10 species with no sequence")
B="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
def get(u):
    for t in range(6):
        try:
            time.sleep(0.5); return json.loads(urllib.request.urlopen(u,timeout=90).read().decode())
        except Exception as e: print("retry",t,e,flush=True); time.sleep(3*(t+1))
    raise SystemExit("failed "+u)
tax={l.split("\t")[1].replace("Candidatus Phytoplasma ","").strip("'\" "):l.split("\t")[0] for l in open(D+"/04_16S/ncbi_species_list.tsv").read().splitlines()[1:]}
for sp in "allocasuarinae australamericanum brasiliense cirsii cocostanzaniae dypsidis gossypii hispanola japonicum periplocae".split():
    ids=get(B+"esearch.fcgi?db=nuccore&retmode=json&retmax=8&term="+urllib.parse.quote("txid%s[Organism:exp] AND (16S OR ribosomal)"%tax[sp]))["esearchresult"]
    print("\n%s: %s hits"%(sp,ids["count"]))
    if ids["idlist"]:
        r=get(B+"esummary.fcgi?db=nuccore&retmode=json&id="+",".join(ids["idlist"]))["result"]
        for k in r["uids"]: print("   %-12s %6s bp  %s"%(r[k]["accessionversion"],r[k]["slen"],r[k]["title"][:95]))
