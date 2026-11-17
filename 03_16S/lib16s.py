# ANALYSIS: 16S rRNA gene analysis
# MANUSCRIPT SECTION: Multilocus sequence comparison and 16S classification
# USED FOR: Figs. S2.1 to S2.7, Tables S2.1 to S2.6
# FOLDER: 03_16S (see README.md in the parent folder)
import time,json,urllib.request
from Bio.Align import PairwiseAligner
D="/projects/dsn011/ppss_publication/2026-10-20_new_species_taxonomy"
B="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
def readfa(p):
    s={};h=None
    for l in open(p):
        l=l.strip()
        if l.startswith(">"): h=l[1:]; s[h]=[]
        elif h: s[h].append(l.upper())
    return {k:"".join(v) for k,v in s.items()}
aln=readfa(D+"/04_16S/16S_allcopies.aln.fasta"); raw={k:v.replace("-","") for k,v in aln.items()}
lab=lambda h:h.rsplit("_c",1)[0]
rc=lambda s:s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
al=PairwiseAligner(); al.mode="global"; al.match_score=2; al.mismatch_score=-3
al.open_gap_score=-6; al.extend_gap_score=-1; al.target_end_gap_score=0; al.query_end_gap_score=0
def pid(a,b):
    best=None
    for st,bb in (("+",b),("-",rc(b))):
        x=al.align(a,bb)[0]
        if best is None or x.score>best[0]: best=(x.score,x,bb,st)
    _,x,bb,st=best; n=m=0
    for (ts,te),(qs,qe) in zip(*x.aligned):
        for p,q in zip(a[ts:te],bb[qs:qe]):
            if p in "ACGT" and q in "ACGT": n+=1; m+=p!=q
    return (100*(1-m/n) if n else 0.0),n,st
def get(u,raw_=False):
    for t in range(6):
        try:
            time.sleep(0.5); d=urllib.request.urlopen(u,timeout=120).read().decode()
            return d if raw_ else json.loads(d)
        except Exception as e: print("retry",t,e,flush=True); time.sleep(3*(t+1))
    raise SystemExit("failed "+u)
def fasta(acc):
    t=get(B+"efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id="+",".join(acc),raw_=True)
    s={};h=None
    for l in t.splitlines():
        if l.startswith(">"): h=l[1:].split()[0]; s[h]=[]
        elif h: s[h].append(l.strip().upper())
    return {k:"".join(v) for k,v in s.items()}
