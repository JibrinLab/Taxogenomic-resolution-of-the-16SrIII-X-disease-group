# ANALYSIS: GTDB-Tk classification of 53 genomes
# MANUSCRIPT SECTION: ANI, dDDH and GTDB clustering
# USED FOR: Tables S5.1 to S5.3
# FOLDER: 06_GTDB-Tk (see README.md in the parent folder)
import openpyxl, sys, re
X="/projects/dsn011/ppss_publication/2026-10-20_new_species_taxonomy/03_ani_curated/Curated83_Supplementary.xlsx"
Q=set("PPS1A PPS2B PPS6A PPS6B PPS9B PPS10A CanS4 MW1 VAC MA1 CX 6A1 ZT3-1 PR2021 JR1 PDA15 CicWB-2022 Vc33 ChTDIII OVSP_L01 OVSP_S01 LWB".split())
ws=openpyxl.load_workbook(X,read_only=True)["Curated83"]
rows=list(ws.iter_rows(values_only=True)); h=rows[0]; ix={k:i for i,k in enumerate(h)}
for col in ("ANIm_cluster95","ANIm_cluster96"):
    cl={}
    for r in rows[1:]:
        g=r[0]; lab=g.split("_",2)[-1] if g.startswith("GC") else g
        cl.setdefault(r[ix[col]],[]).append((g,lab,r[1],r[2]))
    out=open(col+"_reps.tsv","w"); out.write("cluster\tn\tin22\trep\tcompleteness\tcontamination\n")
    for c,m in sorted(cl.items(),key=lambda x:str(x[0])):
        q=[x[1] for x in m if x[1] in Q]
        b=sorted(m,key=lambda x:(-x[2],x[3]))[0]
        out.write(f"{c}\t{len(m)}\t{','.join(q) or '-'}\t{b[0]}\t{b[2]}\t{b[3]}\n")
    out.close(); print(col,len(cl),"clusters")
