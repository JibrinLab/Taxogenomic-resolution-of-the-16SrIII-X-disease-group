# ANALYSIS: Genome tables
# MANUSCRIPT SECTION: Genome sequencing; Additional genomes used
# USED FOR: Table 1 and the master table of the 22 genomes
# FOLDER: 00_genome_tables (see README.md in the parent folder)
import pandas as pd
R='/projects/dsn011/ppss_publication'; A=R+'/2026-10-27_22genome_reanalysis'
S=['PPS1A','PPS2B','PPS6A','PPS6B','PPS9B','PPS10A']
cm=pd.read_csv(R+'/2026-10-20_new_species_taxonomy/01_genome_quality/checkm_lineage_wf_summary.tsv',sep='\t').set_index('Bin Id')
rows=[]
for s in S:
    L=[]; gc=0
    for l in open(R+'/2026-10-20_new_species_taxonomy/16_ddh_upload_PPS1A/'+('PPS1A_query' if s=='PPS1A' else s)+'.fna'):
        if l[0]=='>': L.append(0)
        else:
            q=l.strip().upper(); L[-1]+=len(q); gc+=q.count('G')+q.count('C')
    tot=sum(L); c=0; n50=0
    for x in sorted(L,reverse=True):
        c+=x
        if c>=tot/2: n50=x; break
    cds=trna=0
    for l in open(f'{A}/09_core_tree/pps_pruni22/inputs/{s}.gff'):
        if l[0]=='>': break
        f=l.split('\t')
        if len(f)>2: cds+=f[2]=='CDS'; trna+=f[2]=='tRNA'
    n=d=0
    for l in open(f'{R}/pangenome_v2/coverage/{s}_depth.txt'): n+=1; d+=int(l.split()[2])
    dep=round(d/n,1) if n==tot else None
    if dep is None: print('WARNING depth file for',s,'covers',n,'bp, assembly is',tot)
    rows.append([s,tot,len(L),n50,round(100*gc/tot,2),cds,trna,dep,cm.loc[s,'Completeness'],cm.loc[s,'Contamination']])
df=pd.DataFrame(rows,columns=['Strain','Genome size (bp)','Contigs','N50 (bp)','GC (%)','CDS','tRNA','Mean depth (x)','CheckM completeness (%)','CheckM contamination (%)'])
o=A+'/00_master_genome_table/Table_1_assembly_statistics'
df.to_csv(o+'.tsv',sep='\t',index=False); df.to_excel(o+'.xlsx',index=False)
print(df.to_string(index=False))
