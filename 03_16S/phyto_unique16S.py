#!/usr/bin/env python3
# ANALYSIS: 16S rRNA gene analysis
# MANUSCRIPT SECTION: Multilocus sequence comparison and 16S classification
# USED FOR: Figs. S2.1 to S2.7, Tables S2.1 to S2.6
# FOLDER: 03_16S (see README.md in the parent folder)
"""Unique 16S rRNA gene regions for four 'Candidatus Phytoplasma' species.
1) downloads 16S sequences from NCBI (target species + all other Ca. Phytoplasma species)
2) per target species: finds k-mers present in >= FRAC of its sequences and absent from every sequence
   of the other species, merges them into regions, and reports the position in the reference sequence,
   the sequence, and the fewest mismatches to any other species.
usage: python3 phyto_unique16S.py OUTDIR [--skip-fetch] [--k 15] [--frac 0.9]"""
import sys, os, re, time, json, argparse, urllib.request, urllib.parse
from collections import defaultdict
import numpy as np

TARGETS = {"caryae": "Ca. Phytoplasma caryae", "pruni": "Ca. Phytoplasma pruni",
           "conii": "Ca. Phytoplasma conii", "lactucae": "Ca. Phytoplasma lactucae"}
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

def get(url, tries=4):
    for i in range(tries):
        try:
            return urllib.request.urlopen(url, timeout=120).read().decode()
        except Exception as e:
            err = e; time.sleep(2 * (i + 1))
    raise err

def esearch(term, retmax):
    u = EU + "esearch.fcgi?" + urllib.parse.urlencode({"db": "nuccore", "term": term, "retmax": retmax, "retmode": "json"})
    return json.loads(get(u))["esearchresult"]["idlist"]

def efetch(ids):
    out = []
    for i in range(0, len(ids), 100):
        u = EU + "efetch.fcgi?" + urllib.parse.urlencode({"db": "nuccore", "id": ",".join(ids[i:i + 100]), "rettype": "fasta", "retmode": "text"})
        out.append(get(u)); time.sleep(0.4)
    return "".join(out)

def read_fasta(txt):
    recs = []; h = None; s = []
    for line in txt.splitlines():
        if line.startswith(">"):
            if h: recs.append((h, "".join(s).upper().replace("U", "T")))
            h = line[1:].strip(); s = []
        elif line.strip(): s.append(line.strip())
    if h: recs.append((h, "".join(s).upper().replace("U", "T")))
    return recs

def species_of(h):
    m = re.search(r"Phytoplasma\s+([a-z][a-z\-]+)", h)
    return m.group(1) if m and m.group(1) not in ("sp", "species") else None

def fetch_all(out):
    d = os.path.join(out, "seqs"); os.makedirs(d, exist_ok=True)
    base = " AND (16S[Title] OR 16S ribosomal RNA[Title]) AND 1400:1900[SLEN]"
    for sp in TARGETS:
        ids = esearch('"Candidatus Phytoplasma %s"[Organism]%s' % (sp, base), 500)
        print(sp, len(ids), "records"); open(os.path.join(d, sp + ".fasta"), "w").write(efetch(ids))
    ids = esearch('"Candidatus Phytoplasma"[Organism]' + base, 4000)
    print("all Candidatus Phytoplasma", len(ids), "records")
    open(os.path.join(d, "_all_candidatus.fasta"), "w").write(efetch(ids))

def load(out, cap):
    d = os.path.join(out, "seqs"); S = defaultdict(list)
    for sp in TARGETS:
        f = os.path.join(d, sp + ".fasta")
        if os.path.exists(f): S[sp] += read_fasta(open(f).read())
    f = os.path.join(d, "_all_candidatus.fasta")
    if os.path.exists(f):
        for h, s in read_fasta(open(f).read()):
            sp = species_of(h)
            if sp and sp not in TARGETS and len([1 for x in S[sp]]) < cap: S[sp].append((h, s))
    for sp in list(S):                                   # drop target records not named as that species, short and ambiguous reads
        S[sp] = [(h, s) for h, s in S[sp] if len(s) >= 1400 and sum(c not in "ACGT" for c in s) <= 5]
    # target files may include records of other species or duplicates: keep only same-species names
    for sp in TARGETS:
        keep = [(h, s) for h, s in S[sp] if species_of(h) in (sp, None)]
        S[sp] = keep
    return S

def kmers(s, k): return {s[i:i + k] for i in range(len(s) - k + 1) if set(s[i:i + k]) <= set("ACGT")}

def min_mismatch(region, others):
    L = len(region); r = np.frombuffer(region.encode(), dtype=np.uint8); best = (L, "")
    for sp, seqs in others.items():
        for h, s in seqs:
            if len(s) < L: continue
            a = np.frombuffer(s.encode(), dtype=np.uint8)
            w = np.lib.stride_tricks.sliding_window_view(a, L)
            mm = (w != r).sum(1); j = int(mm.argmin())
            if mm[j] < best[0]: best = (int(mm[j]), "%s|%s" % (sp, h.split()[0]))
    return best

def regions_for(sp, S, k, frac, pool):
    seqs = S[sp]; n = len(seqs)
    if n == 0: return None, []
    ref_h, ref = max(seqs, key=lambda x: len(x[1]))
    cnt = defaultdict(int)
    for h, s in seqs:
        for km in kmers(s, k): cnt[km] += 1
    other = set()
    for o, lst in pool.items():
        for h, s in lst: other |= kmers(s, k)
    ok = {km for km, c in cnt.items() if c >= frac * n and km not in other}
    pos = [i for i in range(len(ref) - k + 1) if ref[i:i + k] in ok]
    regs = []
    for i in pos:
        if regs and i <= regs[-1][1]: regs[-1][1] = i + k
        else: regs.append([i, i + k])
    out = []
    for a, b in regs:
        r = ref[a:b]; nstr = sum(r in s for h, s in seqs)
        out.append((a + 1, b, b - a, r, nstr, n))
    return (ref_h, ref), out

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("out"); ap.add_argument("--skip-fetch", action="store_true")
    ap.add_argument("--k", type=int, default=15); ap.add_argument("--frac", type=float, default=0.9); ap.add_argument("--cap", type=int, default=30)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    if not a.skip_fetch: fetch_all(a.out)
    S = load(a.out, a.cap)
    # sequence inventory
    with open(os.path.join(a.out, "sequences_used.tsv"), "w") as f:
        f.write("species\taccession\tlength\theader\n")
        for sp in sorted(S):
            for h, s in S[sp]: f.write("%s\t%s\t%d\t%s\n" % (sp, h.split()[0], len(s), h))
    rep = []; refs = {}
    for sp, name in TARGETS.items():
        print("%s: %d sequences" % (name, len(S[sp])))
        pool_all = {o: v for o, v in S.items() if o != sp}
        pool_t = {o: S[o] for o in TARGETS if o != sp}
        for tag, pool in (("vs_other_three", pool_t), ("vs_all_other_species", pool_all)):
            ref, regs = regions_for(sp, S, a.k, a.frac, pool)
            if ref is None: print("  no sequences"); continue
            refs[sp] = ref
            for (st, en, L, seq, ns, n) in regs:
                mm, who = min_mismatch(seq, pool)
                rep.append(dict(species=name, comparison=tag, reference=ref[0].split()[0], start=st, end=en, length=L,
                                sequence=seq, strains_with_region="%d of %d" % (ns, n), min_mismatches_to_others=mm, closest_other=who))
            print("  %s: %d unique regions" % (tag, len(regs)))
    import csv
    rep.sort(key=lambda r: (r["species"], r["comparison"], -r["min_mismatches_to_others"], -r["length"]))
    cols = ["species", "comparison", "reference", "start", "end", "length", "sequence", "strains_with_region", "min_mismatches_to_others", "closest_other"]
    with open(os.path.join(a.out, "unique_16S_regions.tsv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(rep)
    with open(os.path.join(a.out, "reference_16S_per_species.fasta"), "w") as f:
        for sp, (h, s) in refs.items(): f.write(">%s\n%s\n" % (h, s))
    print("done ->", a.out)
if __name__ == "__main__": main()
