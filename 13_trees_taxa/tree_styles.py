#!/usr/bin/env python3
# ANALYSIS: Trees with the four species marked
# MANUSCRIPT SECTION: Figures
# USED FOR: Tree figures
# FOLDER: 13_trees_taxa (see README.md in the parent folder)
"""Draw each Newick tree twice (horizontal and circular), one tree per file.
usage: tree_styles.py OUTDIR TREE1 [TREE2 ...]   (names for output = file stem)"""
import sys,os,re,math,matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
matplotlib.rcParams.update({"font.family":"DejaVu Sans","pdf.fonttype":42,"svg.fonttype":"none"})
GROUPS=[("Group 1","#C2410C",["PPS1A","PPS2B","PPS6A","PPS6B","PPS9B","PPS10A","CanS4","MW1","VAC","MA1"]),
        ("Group 2","#2F6DB5",["CX","6A1","ZT3-1","PR2021","JR1","PDA15"]),
        ("Group 3","#2A8C6E",["CicWB-2022","Vc33","ChTDIII","OVSP-L01","OVSP-S01"]),
        ("Group 4","#7E3F98",["LWB"])]
OUTG=re.compile(r"laidlawii|PG[-_]?8A",re.I)
def key(s): return re.sub(r"[^A-Za-z0-9]","",s).upper()
LK={key(n):(g,c,n) for g,c,m in GROUPS for n in m}
def match(name):
    """exact normalised match first, then boundary match of a known name inside the tip label"""
    k=key(name)
    if k in LK: return LK[k]
    n=name.replace("-","_")
    best=None
    for kk,v in LK.items():
        nm=v[2].replace("-","_")
        if re.search(r"(?<![A-Za-z0-9])"+re.escape(nm)+r"(?![A-Za-z0-9])",n):
            if best is None or len(nm)>len(best[2]): best=v
    return best
# ---------------- Newick ----------------
class N:
    __slots__=("dc","ch","name","bl","sup","d","y","ang","par","nt")
    def __init__(s): s.ch=[];s.name="";s.bl=0.0;s.sup="";s.par=None
def parse(t):
    t=t.strip(); t=re.sub(r"\[[^\]]*\]","",t).rstrip(";").strip(); i=0
    def node():
        nonlocal i
        n=N()
        if t[i]=="(":
            i+=1
            while True:
                c=node(); c.par=n; n.ch.append(c)
                if t[i]==",": i+=1; continue
                if t[i]==")": i+=1; break
                raise ValueError("bad newick at %d"%i)
        if i<len(t) and t[i]=="'":
            j=t.index("'",i+1); lab=t[i+1:j]; i=j+1
        else:
            m=re.match(r"[^:,();]*",t[i:]); lab=m.group(0); i+=len(lab)
        if n.ch: n.sup=lab.strip()
        else: n.name=lab.strip()
        if i<len(t) and t[i]==":":
            m=re.match(r":([0-9eE.+-]+)",t[i:]); n.bl=float(m.group(1)); i+=len(m.group(0))
        return n
    return node()
def leaves(n):
    if not n.ch: return [n]
    out=[]
    for c in n.ch: out+=leaves(c)
    return out
# ---------------- rooting (display only) ----------------
def to_graph(root):
    nodes=[];idx={}
    def walk(n):
        idx[id(n)]=len(nodes); nodes.append(n)
        for c in n.ch: walk(c)
    walk(root)
    adj={i:[] for i in range(len(nodes))}
    for n in nodes:
        if n.par is not None:
            a=idx[id(n)];b=idx[id(n.par)]
            adj[a].append((b,n.bl,n.sup)); adj[b].append((a,n.bl,n.sup))
    return nodes,adj
def reroot(root):
    """root on outgroup (PG-8A / laidlawii) if present, else midpoint. returns new root"""
    nodes,adj=to_graph(root)
    tips=[i for i,n in enumerate(nodes) if not n.ch]
    if len(tips)<3: return root
    def dist_from(s):
        d={s:0.0};prev={s:None};st=[s]
        while st:
            u=st.pop()
            for v,l,_ in adj[u]:
                if v not in d: d[v]=d[u]+l;prev[v]=u;st.append(v)
        return d,prev
    og=[i for i in tips if OUTG.search(nodes[i].name)]
    if og:
        o=og[0]; (v,l,lab),=adj[o][:1]; a,b,L,lab=o,v,l,lab; pos=L/2.0 if L>0 else 0.0
        mode="outgroup"
    else:
        d0,_=dist_from(tips[0]); a0=max(tips,key=lambda t:d0[t])
        d1,pv=dist_from(a0); b0=max(tips,key=lambda t:d1[t]); half=d1[b0]/2.0
        path=[b0]
        while path[-1]!=a0: path.append(pv[path[-1]])
        path=path[::-1]  # a0 ... b0
        acc=0.0; a=b=None
        for u,v in zip(path[:-1],path[1:]):
            l=[x for x in adj[u] if x[0]==v][0]
            if acc+l[1]>=half: a,b,L,lab=u,v,l[1],l[2]; pos=half-acc; break
            acc+=l[1]
        if a is None: return root
        mode="midpoint"
    # build new rooted tree with root on edge (a,b), a-side length pos
    R=N()
    def sub(u,excl):
        n=N()
        kids=[x for x in adj[u] if x[0]!=excl]
        if not kids: n.name=nodes[u].name; return n
        for v,l,lb in kids:
            c=sub(v,u); c.bl=l; c.sup=lb if c.ch else ""; c.par=n; n.ch.append(c)
        return n
    ca=sub(a,b); ca.bl=pos; ca.sup=lab if ca.ch else ""; ca.par=R
    cb=sub(b,a); cb.bl=L-pos; cb.sup=lab if cb.ch else ""; cb.par=R
    R.ch=[ca,cb]
    return R
# ---------------- layout ----------------
def layout(root):
    lv=leaves(root)
    def dist(n,d):
        n.d=d+n.bl
        for c in n.ch: dist(c,n.d)
    root.bl=0.0; dist(root,0.0)
    for k,l in enumerate(lv): l.y=k
    def ypos(n):
        if n.ch:
            for c in n.ch: ypos(c)
            n.y=(n.ch[0].y+n.ch[-1].y)/2
    ypos(root)
    return lv
def clean_sup(s):
    s=s.strip()
    if not s: return ""
    try:
        parts=s.split("/")
        return "/".join(("%g"%float(p)) for p in parts)
    except Exception: return s
def sup_first(s):
    try: return float(s.split("/")[-1] if "/" in s else s)
    except Exception: return None
def internal(root):
    out=[];st=[root]
    while st:
        n=st.pop()
        if n.ch: out.append(n); st+=n.ch
    return out
def legend_handles(lv):
    gs=[];seen=set(); other=False
    for l in lv:
        m=match(l.name)
        if m: seen.add(m[0])
        else: other=True
    h=[Patch(fc=c,label=g) for g,c,_ in GROUPS if g in seen]
    if other: h.append(Patch(fc="#777777",label="Other genomes"))
    return h
def dot(ax,x,y,sf,ms):
    if sf is None: ax.plot([x],[y],marker="o",ms=ms,color="black",zorder=4)
    elif sf>=95: ax.plot([x],[y],marker="o",ms=ms,mfc="black",mec="black",zorder=4)
    elif sf>=70: ax.plot([x],[y],marker="o",ms=ms,mfc="#999999",mec="black",mew=0.6,zorder=4)
    else: ax.plot([x],[y],marker="o",ms=ms,mfc="white",mec="black",mew=0.8,zorder=4)
def sup_legend(fig,loc,bbox):
    from matplotlib.lines import Line2D
    h=[Line2D([],[],marker="o",ls="",mfc="black",mec="black",ms=5,label="support \u226595"),
       Line2D([],[],marker="o",ls="",mfc="#999999",mec="black",ms=5,label="70 to 94"),
       Line2D([],[],marker="o",ls="",mfc="white",mec="black",ms=5,label="below 70")]
    fig.legend(handles=h,loc=loc,bbox_to_anchor=bbox,ncol=3,frameon=False,fontsize=6.5)
def bl_txt(x):
    return ("%.4f"%x) if x>=0.001 else ("%.1e"%x)
# ---------------- horizontal ----------------
def draw_h(root,lv,title,base,thr,mode):
    n=len(lv); H=max(4.5,0.30*n+2.2); W=9.5 if n<=40 else 11
    fig=plt.figure(figsize=(W,H)); ax=fig.add_axes([0.04,0.08,0.80,0.84])
    mx=max(l.d for l in lv); lab_x=mx*1.04; show_bl=n<=30
    fs=7 if n<=40 else 5.5
    def draw(nd):
        for c in nd.ch:
            ax.plot([nd.d,nd.d],[nd.y,c.y],color="black",lw=0.9,solid_capstyle="butt")
            ax.plot([nd.d,c.d],[c.y,c.y],color="black",lw=0.9,solid_capstyle="butt")
            draw(c)
    draw(root)
    for nd in internal(root):
        v=clean_sup(nd.sup)
        if nd is root or not v: continue
        sf=sup_first(nd.sup)
        if thr and sf is not None and sf<thr: continue
        dot(ax,nd.d,nd.y,sf,3.6 if n<=40 else 3.0)
        ax.text(nd.d-mx*0.008,nd.y-0.18,v,ha="right",va="bottom",fontsize=6.2 if n<=40 else 4.8,color="black",zorder=5,bbox=dict(fc="white",ec="none",pad=0.5))
    for l in lv:
        m=match(l.name); c=m[1] if m else "#777777"; disp=m[2] if m else l.name
        if lab_x-l.d>mx*0.01: ax.plot([l.d,lab_x],[l.y,l.y],color="#BDBDBD",lw=0.5,zorder=0)
        ax.scatter([lab_x],[l.y],s=14,color=c,zorder=3)
        ax.text(lab_x+mx*0.012,l.y,disp,ha="left",va="center",fontsize=fs,color="black")
    ax.set_ylim(n+0.9,-0.8); ax.set_xlim(-mx*0.02,mx*1.35); ax.axis("off")
    # scale bar
    sb=10**math.floor(math.log10(mx/4)); 
    ax.plot([0,sb],[n-0.2,n-0.2],color="black",lw=1.2); ax.text(sb/2,n+0.1,"%g substitutions per site"%sb,ha="center",va="top",fontsize=6.5,color="black")
    fig.text(0.04,0.975,title,fontsize=9,fontweight="bold",va="top")
    hh=legend_handles(lv)
    if hh: fig.legend(handles=hh,loc="upper right",bbox_to_anchor=(0.98,0.985),ncol=len(hh),frameon=False,fontsize=7)
    note="Branch lengths are drawn to scale (substitutions per site); %s. Node values are branch support%s."%(mode,", values below %g not shown"%thr if thr else "")
    fig.text(0.04,0.012,note,fontsize=5.8,color="black")
    sup_legend(fig,"lower right",(0.98,0.03)); save(fig,base+"_horizontal")
# ---------------- circular ----------------
def draw_c(root,lv,title,base,thr,mode):
    n=len(lv); size=8.5 if n<=40 else 11
    fig=plt.figure(figsize=(size,size)); ax=fig.add_axes([0.04,0.06,0.92,0.88],projection="polar"); ax.axis("off")
    for k,l in enumerate(lv): l.ang=2*math.pi*k/n
    def setang(nd):
        if nd.ch:
            for c in nd.ch: setang(c)
            nd.ang=(nd.ch[0].ang+nd.ch[-1].ang)/2
    setang(root)
    ds=sorted(l.d for l in lv); cap=ds[-1]; trunc=False
    if len(ds)>2 and ds[-1]>1.4*ds[-2]: cap=ds[-2]*1.05; trunc=True
    def setdc(nd):
        nd.dc=min(nd.d,cap)
        for c in nd.ch: setdc(c)
    setdc(root); mx=cap
    ax.set_ylim(0,mx*1.65); fs=7 if n<=40 else 4.8
    def draw(nd):
        for c in nd.ch:
            ax.plot([c.ang,c.ang],[nd.dc,c.dc],color="black",lw=0.9,solid_capstyle="butt")
        if nd.ch:
            a=[c.ang for c in nd.ch]; th=np.linspace(min(a),max(a),80); ax.plot(th,[nd.dc]*80,color="black",lw=0.9)
            for c in nd.ch: draw(c)
    draw(root)
    for nd in internal(root):
        v=clean_sup(nd.sup)
        if nd is root or not v: continue
        sf=sup_first(nd.sup)
        if thr and sf is not None and sf<thr: continue
        dot(ax,nd.ang,nd.dc,sf,3.4 if n<=40 else 2.8)
        ax.text(nd.ang,nd.dc+mx*0.035,v,ha="center",va="center",fontsize=5.6 if n<=40 else 4.2,color="black",zorder=5,bbox=dict(fc="white",ec="none",pad=0.4))
    for l in lv:
        m=match(l.name); c=m[1] if m else "#777777"; disp=m[2] if m else l.name
        ax.plot([l.ang,l.ang],[l.dc,mx*1.05],color="#BDBDBD",lw=0.5,zorder=0)
        ax.scatter([l.ang],[mx*1.07],s=14,color=c,zorder=3)
        rot=math.degrees(l.ang); flip=90<rot<270
        ax.text(l.ang,mx*1.11,disp,rotation=rot+180 if flip else rot,rotation_mode="anchor",ha="right" if flip else "left",va="center",fontsize=fs,color="black")
    fig.text(0.04,0.975,title,fontsize=9,fontweight="bold",va="top")
    hh=legend_handles(lv)
    if hh: fig.legend(handles=hh,loc="lower center",bbox_to_anchor=(0.5,0.015),ncol=len(hh),frameon=False,fontsize=7.5)
    sb=10**math.floor(math.log10(mx/4))
    from matplotlib.lines import Line2D as _L
    ln=sb/(mx*1.65)*0.46
    fig.add_artist(_L([0.05,0.05+ln],[0.085,0.085],transform=fig.transFigure,color="black",lw=1.6))
    fig.text(0.05,0.077,"%g substitutions per site"%sb,ha="left",va="top",fontsize=6.5)
    note="Branch lengths are drawn to scale (substitutions per site); %s. Node values are branch support%s%s."%(mode,", values below %g not shown"%thr if thr else "","; outgroup branch truncated" if trunc else "")
    fig.text(0.04,0.003,note,fontsize=5.8,color="black")
    sup_legend(fig,"upper right",(0.98,0.985)); save(fig,base+"_circular")
def save(fig,base):
    for e in ("pdf","svg"): fig.savefig(base+"."+e)
    fig.savefig(base+".png",dpi=300); plt.close(fig)
def main():
    out=sys.argv[1]; os.makedirs(out,exist_ok=True); ok=0
    for p in sys.argv[2:]:
        try:
            stem=os.path.splitext(os.path.basename(p))[0]
            if len(sys.argv)>3 or True:
                parent=os.path.basename(os.path.dirname(p))
                stem=(parent+"__"+stem) if parent and parent not in ("",".") else stem
            root=parse(open(p).read().split(";")[0]+";")
            lv0=leaves(root)
            if len(lv0)<3: print("SKIP (<3 tips):",p); continue
            if len(lv0)>400: print("SKIP (>400 tips):",p); continue
            og=[l for l in lv0 if OUTG.search(l.name)]
            root=reroot(root); mode="rooted on the outgroup" if og else "midpoint-rooted for display"
            lv=layout(root)
            thr=70 if len(lv)>60 else 0
            title=stem.replace("__"," / ").replace("_"," ")
            base=os.path.join(out,stem)
            draw_h(root,lv,title,base,thr,mode); draw_c(root,lv,title,base,thr,mode)
            print("ok  %-60s tips=%d  matched22=%d"%(stem,len(lv),sum(1 for l in lv if match(l.name)))); ok+=1
        except Exception as e:
            print("FAILED",p,repr(e))
    print("done:",ok,"trees ->",out)
if __name__=="__main__": main()
