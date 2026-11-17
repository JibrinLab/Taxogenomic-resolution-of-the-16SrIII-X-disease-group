#!/usr/bin/env python3
# ANALYSIS: Trees with the four species marked
# MANUSCRIPT SECTION: Figures
# USED FOR: Tree figures
# FOLDER: 13_trees_taxa (see README.md in the parent folder)
"""Rectangular tree with the four proposed taxa assigned (no title, no legend, no branch lengths).
Tip circles are coloured by taxon; when all tips of a taxon form one clade a bracket with the species name
is drawn on the right. Other taxa are grey. Support values only. Writes PNG (600 dpi), PDF, SVG and
taxa_assignment_summary.tsv (is each taxon one clade in this tree, and its support).
Same folder as tree_styles.py.
usage: tree_taxa.py [--map MAP.tsv] [--groups-only] OUTDIR TREE [TREE ...]
--groups-only: keep only the 22 genomes of the four groups (plus the outgroup) and label the four taxa only.
Child branches are flipped (topology unchanged) so that tips of one taxon sit together; a taxon that is not one
clade gets a bracket over each block of neighbouring tips, and taxa_assignment_summary.tsv says one_clade=no.
MAP.tsv (optional, tab-separated, header tip<TAB>taxon): species name of every other tip, so that the
referenced strains also get their taxon bracket. Tips not in the map stay unlabelled."""
import sys, os, math, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tree_styles as T
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
TAXA = {"Group 1": ("Ca. P. caryae", "#4F7FB8"), "Group 2": ("Ca. P. pruni", "#F28E2B"),
        "Group 3": ("Ca. P. manihotii", "#59A14F"), "Group 4": ("Ca. P. lactucae", "#B07AA1")}
def lca(n, S):
    for c in n.ch:
        if S <= {id(x) for x in T.leaves(c)}: return lca(c, S)
    return n
def load_map(p):
    m = {}
    for i, l in enumerate(open(p)):
        a = l.rstrip("\n").split("\t")
        if i == 0 and a[0].lower() == "tip": continue
        if len(a) >= 2 and a[1].strip(): m[T.key(a[0])] = re.sub(r"\bconii\b", "manihotii", a[1].strip())
    return m

def prune(root, keep):
    """drop leaves for which keep(name) is False; collapse single-child nodes (branch lengths added)"""
    def rec(n):
        if not n.ch: return n if keep(n.name) else None
        kids = [k for k in (rec(c) for c in n.ch) if k is not None]
        if not kids: return None
        if len(kids) == 1:
            k = kids[0]; k.bl += n.bl; return k
        n.ch = kids
        for k in kids: k.par = n
        return n
    r = rec(root); r.par = None; return r
def arrange(root, lab):
    """flip child order (topology unchanged) so that tips of the same taxon sit next to each other:
    hill-climb on the number of neighbouring tip pairs with different labels"""
    import random
    nodes = [n for n in T.internal(root) if len(n.ch) > 1]
    def breaks():
        lv = T.leaves(root); L = [lab(x) for x in lv]
        return sum(1 for a, b in zip(L, L[1:]) if a != b or a is None)
    best = breaks(); random.seed(1)
    for _ in range(40):
        improved = False
        order = nodes[:]; random.shuffle(order)
        for nd in order:
            nd.ch.reverse(); b = breaks()
            if b < best: best = b; improved = True
            else: nd.ch.reverse()
        if not improved: break
    for nd in nodes:  # tie-break: larger clade first
        pass
    return best
def draw(path, outdir, rows, cmap, groups_only=False):
    root = T.parse(open(path).read().split(";")[0] + ";")
    if groups_only:
        root = prune(root, lambda nm: bool(T.match(nm)) or bool(T.OUTG.search(nm)))
    root = T.reroot(root)
    def lab(l):
        m = T.match(l.name)
        if m: return TAXA[m[0]][0]
        k0 = T.key(l.name); return None if groups_only else (cmap.get(k0) or cmap.get(re.sub(r"C\d+$", "", k0)))
    nbreaks = arrange(root, lab)
    lv = T.layout(root); n = len(lv)
    mx = max(l.d for l in lv) or 1.0; fs = 8 if n <= 40 else 5.5
    fig = plt.figure(figsize=(8 if n <= 40 else 9.5, max(4, 0.27 * n + 1.2)))
    ax = fig.add_axes([0.03, 0.04, 0.80, 0.92])
    def rec(nd):
        for c in nd.ch:
            ax.plot([nd.d, nd.d], [nd.y, c.y], color="black", lw=1.0, solid_capstyle="butt")
            ax.plot([nd.d, c.d], [c.y, c.y], color="black", lw=1.0, solid_capstyle="butt")
            rec(c)
    rec(root)
    for nd in T.internal(root):
        v = T.clean_sup(nd.sup)
        if nd is root or not v: continue
        ax.text(nd.d - mx * 0.008, nd.y - 0.15, v, ha="right", va="bottom", fontsize=fs - 1.5, color="black")
    gof = {}; texts = []; other = {}
    for l in lv:
        m = T.match(l.name); g = m[0] if m else None; disp = m[2] if m else l.name
        if m and l.name != m[2]: disp = l.name
        gof[id(l)] = g
        if not g:
            k0 = T.key(l.name); other[id(l)] = cmap.get(k0) or cmap.get(re.sub(r"C\d+$", "", k0))
        ax.scatter([l.d], [l.y], s=40 if n <= 40 else 14, color=TAXA[g][1] if g else "#B0B0B0", edgecolor="black", linewidth=0.6, zorder=3)
        texts.append(ax.text(l.d + mx * 0.015, l.y, disp, ha="left", va="center", fontsize=fs, color="black",
                             fontweight="bold" if disp.startswith("PPS") else "normal"))
    ax.set_ylim(n + 0.9, -0.8); ax.set_xlim(-mx * 0.12, mx * 1.28); ax.axis("off")
    fig.canvas.draw(); rend = fig.canvas.get_renderer(); inv = ax.transData.inverted()
    xr = max(inv.transform((t.get_window_extent(rend).x1, 0))[0] for t in texts) + mx * 0.03
    tn = os.path.basename(path)
    # taxon of every tip: the four proposed taxa first, then the map (other referenced strains)
    tax = {}
    for l in lv:
        g = gof[id(l)]
        tax[id(l)] = TAXA[g][0] if g else other.get(id(l))
    names = []
    for l in lv:
        t = tax[id(l)]
        if t and t not in names: names.append(t)
    own = {v[0] for v in TAXA.values()}
    for t in names:
        S = {id(l) for l in lv if tax[id(l)] == t}
        node = lca(root, S); inside = {id(x) for x in T.leaves(node)}
        mono = inside == S
        # blocks of neighbouring tips with the same taxon (after the branch flip)
        blocks = []; cur = []
        for l in lv:
            if tax[id(l)] == t: cur.append(l.y)
            elif cur: blocks.append(cur); cur = []
        if cur: blocks.append(cur)
        rows.append((tn, t, len(S), "yes" if mono else "no", T.clean_sup(node.sup) if (mono and node.ch) else "",
                     len(inside) - len(S) if not mono else 0, len(blocks)))
        for ys in blocks:
            solid = t in own
            ax.plot([xr, xr], [min(ys) - 0.35, max(ys) + 0.35], color="black" if solid else "#606060",
                    lw=1.4 if solid else 1.0, solid_capstyle="butt", clip_on=False)
            ax.text(xr + mx * 0.02, (min(ys) + max(ys)) / 2, t, ha="left", va="center",
                    fontsize=fs + 1 if solid else fs, style="italic", color="black", clip_on=False)
    ax.set_xlim(-mx * 0.12, xr + mx * 0.30)
    sb = 10 ** math.floor(math.log10(mx / 4))
    ax.plot([0, sb], [n + 0.3, n + 0.3], color="black", lw=1.2)
    ax.text(sb / 2, n + 0.55, "%g" % sb, ha="center", va="top", fontsize=fs - 1, color="black")
    base = os.path.join(outdir, os.path.splitext(os.path.basename(path))[0] + "_taxa")
    for ext, kw in (("png", {"dpi": 600}), ("pdf", {}), ("svg", {})):
        fig.savefig(base + "." + ext, bbox_inches="tight", facecolor="white", **kw)
    plt.close(fig); print("ok", base)
if __name__ == "__main__":
    a = sys.argv[1:]; cmap = {}; go = False
    while a and a[0] in ("--map", "--groups-only"):
        if a[0] == "--map": cmap = load_map(a[1]); a = a[2:]
        else: go = True; a = a[1:]
    sys.argv = [sys.argv[0]] + a
    os.makedirs(sys.argv[1], exist_ok=True); rows = []
    for p in sys.argv[2:]:
        try: draw(p, sys.argv[1], rows, cmap, go)
        except Exception as e: print("FAILED", p, repr(e))
    with open(os.path.join(sys.argv[1], "taxa_assignment_summary.tsv"), "w") as f:
        f.write("tree\ttaxon\ttips\tone_clade\tclade_support\textra_tips_inside_span\tdisplay_blocks\n")
        for r in rows: f.write("\t".join(str(x) for x in r) + "\n")
    print(open(os.path.join(sys.argv[1], "taxa_assignment_summary.tsv")).read())
