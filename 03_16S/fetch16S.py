# ANALYSIS: 16S rRNA gene analysis
# MANUSCRIPT SECTION: Multilocus sequence comparison and 16S classification
# USED FOR: Figs. S2.1 to S2.7, Tables S2.1 to S2.6
# FOLDER: 03_16S (see README.md in the parent folder)
import urllib.request,urllib.parse,json,time,sys
D="/projects/dsn011/ppss_publication/2026-10-20_new_species_taxonomy"
B="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
def get(u,raw=False):
    for t in range(6):
        try:
            time.sleep(0.5)
            d=urllib.request.urlopen(u,timeout=90).read().decode()
            return d if raw else json.loads(d)
        except Exception as e:
            print("retry",t,e,flush=True); time.sleep(3*(t+1))
    raise SystemExit("failed: "+u)
miss=[l.rstrip("\n").split("\t") for l in open(D+"/04_16S/ncbi_species_list.tsv").read().splitlines()[1:]]
miss=[m for m in miss if m[2]=="MISSING"]
print("species to query:",len(miss),flush=True)
cand=open(D+"/04_16S/described_missing_candidates.tsv","w")
cand.write("species\ttaxid\taccession\tlength\ttitle\texact_species_title\tchosen\n")
chosen={}
for tid,name,_ in miss:
    sp=name.replace("Candidatus Phytoplasma ","").strip("'\" ")
    print("esearch",sp,flush=True)
    term='txid%s[Organism:exp] AND 16S AND 1200:1800[SLEN]'%tid
    ids=get(B+"esearch.fcgi?db=nuccore&retmode=json&retmax=100&term="+urllib.parse.quote(term))["esearchresult"]["idlist"]
    if not ids:
        print("  no hits",flush=True); cand.write("%s\t%s\t-\t-\tNO HITS\t\t\n"%(sp,tid)); continue
    print("  esummary",len(ids),"ids",flush=True)
    r=get(B+"esummary.fcgi?db=nuccore&retmode=json&id="+",".join(ids))["result"]
    c=[]
    for k in r["uids"]:
        t=r[k]["title"]; tl=t.lower()
        if "16s ribosomal rna" not in tl or "uncultured" in tl or "23s" in tl or "spacer" in tl: continue
        c.append((r[k]["accessionversion"],int(r[k]["slen"]),t,("phytoplasma "+sp) in tl))
    if not c:
        print("  none passed title filter",flush=True); cand.write("%s\t%s\t-\t-\tNONE PASSED TITLE FILTER\t\t\n"%(sp,tid)); continue
    best=max(c,key=lambda x:(x[3],x[1])); chosen[sp]=best
    for a,l,t,e in c: cand.write("%s\t%s\t%s\t%d\t%s\t%s\t%s\n"%(sp,tid,a,l,t,"Y" if e else "N","YES" if a==best[0] else ""))
    print("  chosen %-12s %4d bp exact-name:%s (%d candidates)"%(best[0],best[1],"Y" if best[3] else "N",len(c)),flush=True)
cand.close()
print("efetch",len(chosen),"sequences",flush=True)
fa=get(B+"efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id="+",".join(v[0] for v in chosen.values()),raw=True)
out=[]
for l in fa.splitlines():
    if l.startswith(">"):
        a=l[1:].split()[0]; sp=[k for k,v in chosen.items() if v[0]==a][0]; out.append(">%s|%s"%(sp,a))
    else: out.append(l)
open(D+"/04_16S/described_missing_16S.fasta","w").write("\n".join(out)+"\n")
print("fasta records:",sum(1 for l in out if l.startswith(">")),"of",len(miss),flush=True)
