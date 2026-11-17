# ANALYSIS: 16S rRNA gene analysis
# MANUSCRIPT SECTION: Multilocus sequence comparison and 16S classification
# USED FOR: Figs. S2.1 to S2.7, Tables S2.1 to S2.6
# FOLDER: 03_16S (see README.md in the parent folder)
import csv,openpyxl
from openpyxl.styles import Font,PatternFill,Alignment
D="/projects/dsn011/ppss_publication/2026-10-20_new_species_taxonomy"
Q=[("PPS1A","caryae (PPS1A)"),("OVSP_L01","manihotii (OVSP-L01)"),("CicWB_2022","manihotii (CicWB-2022)"),("LWB","lactucae (LWB)"),("CX","pruni (CX)")]
best={}
for r in csv.DictReader(open(D+"/04_16S/all_described_species_identity_combined.tsv"),delimiter="\t"):
    if r["species"]=="pruni" and r["subject"]!="JQ044393.1": continue
    k=(r["species"],r["query"]); v=(float(r["identity_pct"]),int(r["sites"]),r["subject"])
    if k not in best or v[0]>best[k][0]: best[k]=v
rows=[]
for l in open(D+"/04_16S/ncbi_species_list.tsv").read().splitlines()[1:]:
    tid,name,_=l.split("\t"); sp=name.replace("Candidatus Phytoplasma ","").strip("'\" ")
    vals=[best.get((sp,q)) for q,_ in Q]
    rows.append((sp,vals))
rows.sort(key=lambda r:-max([v[0] for v in r[1][:4] if v] or [0]))
wb=openpyxl.Workbook(); ws=wb.active; ws.title="16S identity"
ws.append(["Described species (Ca. Phytoplasma)","16S sequence used"]+[n+" identity (%)" for _,n in Q]+["Aligned sites (PPS1A)"])
yel=PatternFill("solid",fgColor="FFFF00")
for sp,vals in rows:
    first=next((v for v in vals if v),None)
    ws.append([sp,first[2] if first else "no sequence found"]+[(round(v[0],2) if v else ("self" if (q=="CX" and sp=="pruni") else "")) for (q,_),v in zip(Q,vals)]+[vals[0][1] if vals[0] else ""])
for row in ws.iter_rows(min_row=2):
    for c in row[2:7]:
        if isinstance(c.value,float) and c.value>=98.65: c.fill=yel
for c in ws[1]: c.font=Font(bold=True); c.alignment=Alignment(wrap_text=True)
for i,w in enumerate([22,24,16,18,18,16,14,14],1): ws.column_dimensions[openpyxl.utils.get_column_letter(i)].width=w
ws2=wb.create_sheet("Notes")
for t in ["Identity = matching bases / aligned positions with A, C, G or T in both sequences. Query = first 16S copy of one genome per taxon.","Highlighted cells are 98.65% or higher.","Species are the 62 species-rank names in NCBI Taxonomy (12 November 2026). The GenBank sequence used for each species is not confirmed to be its type or reference strain.","Short sequences (about 1,100 to 1,200 bp) cover less of the gene; Ca. P. gossypii is 1,023 bp from a draft genome."]: ws2.append([t])
ws2.column_dimensions["A"].width=140
wb.save("/projects/dsn011/ppss_publication/2026-10-20_new_species_taxonomy/26_jibrin_update_2026-11-12/results/01_16S_identity_table.xlsx")
n=sum(1 for _,v in rows if any(v)); print("species with a value:",n,"of",len(rows))
for i,(q,nm) in enumerate(Q):
    hi=sorted([(v[i][0],sp) for sp,v in rows if v[i] and sp!="pruni"],reverse=True)
    print("%-24s highest other than pruni: %s | species >=98.65: %s"%(nm,"%s %.2f"%(hi[0][1],hi[0][0]) if hi else "-",[sp for sp,v in rows if v[i] and v[i][0]>=98.65]))
