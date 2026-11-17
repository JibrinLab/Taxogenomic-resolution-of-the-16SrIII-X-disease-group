# ANALYSIS: Genome tables
# MANUSCRIPT SECTION: Genome sequencing; Additional genomes used
# USED FOR: Table 1 and the master table of the 22 genomes
# FOLDER: 00_genome_tables (see README.md in the parent folder)
import os, re, glob
ROOT='/projects/dsn011/ppss_publication'
TAX=ROOT+'/2026-10-20_new_species_taxonomy'
LABELS=['PPS1A','PPS2B','PPS6A','PPS6B','PPS9B','PPS10A','CanS4','MW1','VAC','MA1','CX','6A1','ZT3-1','PR2021','JR1','PDA15','CicWB-2022','Vc33','ChTDIII','OVSP_L01','OVSP_S01','LWB','A_laidlawii']
ALIAS={'A_laidlawii':['PG8A','PG-8A']}
SKIP={'ani_all299','ani_23_with_PG8A','2026-10-27_22genome_reanalysis','__pycache__','check_tar_urbello','removed_from_package_1025','PPS_taxonomy_package_for_Dr_Jibrin_2026-10-22'}
EXT=('.faa','.gff','.gbk','.gbff','.fna','.fasta','.fa','.ffn')
files=[]
for dp,dn,fn in os.walk(ROOT):
    dn[:]=[d for d in dn if d not in SKIP]
    for f in fn:
        if f.endswith(EXT):
            files.append(os.path.join(dp,f))
print('files scanned:',len(files))
out=open('all_matches.tsv','w')
out.write('label\text\tpath\n')
cnt={}
for lab in LABELS:
    names=[lab]+ALIAS.get(lab,[])
    pat=re.compile('|'.join('(?<![A-Za-z0-9])'+re.escape(x)+'(?![A-Za-z0-9])' for x in names))
    for p in files:
        key=os.path.basename(os.path.dirname(p))+'/'+os.path.basename(p)
        if pat.search(key):
            e=os.path.splitext(p)[1]
            out.write(lab+'\t'+e+'\t'+p+'\n')
            cnt[(lab,e)]=cnt.get((lab,e),0)+1
out.close()
exts=['.faa','.gff','.gbk','.gbff','.fna','.fasta','.fa','.ffn']
print('label\t'+'\t'.join(exts))
for lab in LABELS:
    print(lab+'\t'+'\t'.join(str(cnt.get((lab,e),0)) for e in exts))
def acc(lab):
    if lab=='A_laidlawii':
        return 'NC_010163.1'
    g=glob.glob(ROOT+'/ani_all299/genomes/*_'+lab+'.fna')
    return os.path.basename(g[0])[:-len('_'+lab+'.fna')] if g else 'NA'
def fa(lab):
    if lab=='A_laidlawii':
        return ROOT+'/00_genomes/gtdbtk_batches/all_genomes_fixed/PG8A_Acholeplasma.laidlawii.fasta'
    n='PPS1A_query' if lab=='PPS1A' else lab
    return TAX+'/16_ddh_upload_PPS1A/'+n+'.fna'
def stats(p):
    L=[]; cur=0; gc=0
    for l in open(p):
        if l.startswith('>'):
            if cur: L.append(cur)
            cur=0
        else:
            s=l.strip().upper(); cur+=len(s); gc+=s.count('G')+s.count('C')
    if cur: L.append(cur)
    L.sort(reverse=True); tot=sum(L); run=0; n50=0
    for x in L:
        run+=x
        if run>=tot/2:
            n50=x; break
    return tot,len(L),100.0*gc/tot,n50
t=open('master_genome_table.tsv','w')
t.write('label\taccession\tfasta\tsize_bp\tcontigs\tGC_percent\tN50\n')
for lab in LABELS:
    p=fa(lab)
    if os.path.exists(p):
        a,b,c,d=stats(p)
        t.write('%s\t%s\t%s\t%d\t%d\t%.2f\t%d\n'%(lab,acc(lab),p,a,b,c,d))
    else:
        t.write('%s\t%s\tMISSING\t\t\t\t\n'%(lab,acc(lab)))
t.close()
print(open('master_genome_table.tsv').read())
