# ANALYSIS: 16S rRNA gene analysis
# MANUSCRIPT SECTION: Multilocus sequence comparison and 16S classification
# USED FOR: Figs. S2.1 to S2.7, Tables S2.1 to S2.6
# FOLDER: 03_16S (see README.md in the parent folder)
from lib16s import *
Q=("PPS1A_c1","CX_c1","OVSP_L01_c1","CicWB_2022_c1","LWB_c1")
pick={"allocasuarinae":("AY725236.1",""),"brasiliense":("PP581850.1","clone chosen arbitrarily"),"cirsii":("KR869147.1","16S+tRNA-Ile"),
 "dypsidis":("MT536196.1","16S+ISR"),"japonicum":("AB010425.2",""),"periplocae":("GQ249159.1","title names host not species"),
 "australamericanum":("GU292081.1","record titled sudamericanum")}
seq=fasta([v[0] for v in pick.values()])
out=open(D+"/04_16S/identity_vs_described_species_extra2.tsv","w")
out.write("query\tspecies\taccession\tidentity_pct\tsites\tstrand\tlength\tnote\n")
for sp,(ac,note) in pick.items():
    k=[x for x in seq if x.startswith(ac.split(".")[0])]
    if not k: print("missing",sp); continue
    r=[]
    for q in Q:
        i,n,s=pid(raw[q],seq[k[0]]); out.write("%s\t%s\t%s\t%.2f\t%d\t%s\t%d\t%s\n"%(lab(q),sp,ac,i,n,s,len(seq[k[0]]),note)); r.append("%s %.2f%%(%d)"%(lab(q),i,n))
    print("%-18s %-11s %s"%(sp,ac," | ".join(r)))
out.close()
rows=[l.split("\t")[:5]+["genome set"] for l in open(D+"/04_16S/identity_vs_described_species.tsv").read().splitlines()[1:]]
for fn in ("extra","extra2"):
    rows+=[l.split("\t")[:5]+[fn] for l in open(D+"/04_16S/identity_vs_described_species_%s.tsv"%fn).read().splitlines()[1:]]
open(D+"/04_16S/all_described_species_identity_combined.tsv","w").write("query\tspecies\tsubject\tidentity_pct\tsites\tsource\n"+"\n".join("\t".join(r) for r in rows)+"\n")
for q in ("PPS1A","CX","OVSP_L01","CicWB_2022","LWB"):
    sp={r[1] for r in rows if r[0]==q}; hi=sorted({r[1] for r in rows if r[0]==q and float(r[3])>=98.65})
    print("%-11s species: %d | at or above 98.65%%: %s"%(q,len(sp),hi))
