# ANALYSIS: 16S rRNA gene analysis
# MANUSCRIPT SECTION: Multilocus sequence comparison and 16S classification
# USED FOR: Figs. S2.1 to S2.7, Tables S2.1 to S2.6
# FOLDER: 03_16S (see README.md in the parent folder)
from lib16s import *
import shutil,subprocess,urllib.parse
tid=[l.split("\t")[0] for l in open(D+"/04_16S/ncbi_species_list.tsv") if "gossypii" in l][0]
ids=get(B+"esearch.fcgi?db=nuccore&retmode=json&retmax=200&term="+urllib.parse.quote("txid%s[Organism:exp] AND JAOSIT01[All Fields]"%tid))["esearchresult"]["idlist"]
g={}
for i in range(0,len(ids),40): g.update(fasta(ids[i:i+40]))
g={k:v for k,v in g.items() if k.startswith("NZ_JAOSIT01")}
print("contigs:",len(g),"bp:",sum(map(len,g.values())))
open("gossypii.fa","w").write("".join(">%s\n%s\n"%kv for kv in g.items()))
open("q16S.fa","w").write(">q\n"+raw["PPS1A_c1"]+"\n")
print("blastn:",shutil.which("blastn"))
r=subprocess.run(["blastn","-query","q16S.fa","-subject","gossypii.fa","-outfmt","6 sseqid sstart send length pident"],capture_output=True,text=True)
hits=sorted([l.split("\t") for l in r.stdout.strip().splitlines()],key=lambda x:-int(x[3]))[:4]
for h in hits: print("hit",h)
c,s,e=hits[0][0],int(hits[0][1]),int(hits[0][2]); sq=g[c][min(s,e)-1:max(s,e)]; sq=rc(sq) if s>e else sq
out=open(D+"/04_16S/all_described_species_identity_combined.tsv","a")
for q in ("PPS1A_c1","CX_c1","OVSP_L01_c1","CicWB_2022_c1","LWB_c1"):
    i,n,st=pid(raw[q],sq); out.write("%s\tgossypii\t%s\t%.2f\t%d\textra3\n"%(lab(q),c,i,n)); print(lab(q),"%.2f%% (%d sites) extracted length %d"%(i,n,len(sq)))
out.close()
