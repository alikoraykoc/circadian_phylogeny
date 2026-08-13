#!/usr/bin/env python3
"""Count the convergent substitutions that ACTUALLY occurred, by parsimony.

transition_opportunity.py estimates opportunity from branch lengths, which
assumes a uniform per-site rate. That is pessimistic: most columns in these
proteins are invariant, so averaging over them understates how much freedom the
sites that CAN change actually have. This script removes the assumption and
counts directly.

Fitch parsimony assigns an amino acid to every internal node of each gene tree,
per column. A substitution on a branch is simply parent state != child state.
That gives, with no evolutionary model and no rate assumption:

  opportunity  how many sites changed on the transition branch of at least two
               independent diurnal events. Convergence is undefined below two,
               so this is the honest denominator for the whole study.
  convergence  of those, how many changed to the SAME residue in two or more
               independent events. This is the textbook definition of convergent
               substitution, and it is what PCOC's OneChange component keys on.

Parallel changes to the same residue are also counted against a null: the same
tally computed after re-placing the diurnal clades at random on the same tree,
which preserves phylogenetic clumping (see model_free_convergence_screen.py for
why a plain label shuffle would be wrong).

Output: results/pcoc_sim/observed_convergent_substitutions.csv
        results/pcoc_sim/observed_convergent_sites.csv  (the actual hits)
"""
import csv
import os
import random
from collections import defaultdict

from ete3 import Tree

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TRIM = os.path.join(PROJ, "results", "trim")
TREES = os.path.join(PROJ, "results", "branchlengths", "pcoc")
SCEN = os.path.join(PROJ, "results", "pcoc", "scenarios", "ER_gain")
OUT = os.path.join(PROJ, "results", "pcoc_sim")

GENES = ("CLOCK NPAS2 ARNTL PER1 PER2 PER3 CRY1 CRY2 NR1D1 NR1D2 "
         "RORA RORB RORC CSNK1D CSNK1E FBXL3 BHLHE40 BHLHE41").split()

GAPS = set("-X?*BZJU")
N_PERM = 200
rng = random.Random(20260813)


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        elif line:
            buf.append(line)
    if name:
        seqs[name] = "".join(buf)
    return seqs


def numbered(tree_file):
    """PCOC node numbering: postorder from 0 (events_placing.py:init_tree)."""
    t = Tree(tree_file, format=1)
    by_id = {}
    for i, n in enumerate(t.traverse("postorder")):
        n.add_features(ND=i)
        by_id[i] = n
    return t, by_id


def fitch(tree, col_state):
    """Fitch parsimony. Returns {node: assigned residue or None}.

    Gaps and ambiguous characters are treated as missing: they contribute
    nothing to the state sets rather than counting as a 21st residue, so a
    gapped clade cannot manufacture a spurious substitution.
    """
    up = {}
    for n in tree.traverse("postorder"):
        if n.is_leaf():
            ch = col_state.get(n.name)
            up[n] = {ch} if ch and ch not in GAPS else set()
        else:
            sets = [up[c] for c in n.children if up[c]]
            if not sets:
                up[n] = set()
            else:
                inter = set.intersection(*sets)
                up[n] = inter if inter else set().union(*sets)
    out = {}
    for n in tree.traverse("preorder"):
        s = up[n]
        if not s:
            out[n] = None
            continue
        parent = n.up
        if parent is not None and out.get(parent) in s:
            out[n] = out[parent]
        else:
            out[n] = sorted(s)[0]
        # keep assignment deterministic for reproducibility
    return out


def tally(tree, by_id, events, seqs, ncol, order):
    """Per site: which events changed on their transition branch, and to what."""
    opp = 0
    conv_sites = []
    for c in range(ncol):
        col = {nm: seqs[nm][c] for nm in order}
        assign = fitch(tree, col)
        derived = {}
        for ei, ids in enumerate(events):
            node = by_id[ids[0]]
            parent = node.up
            if parent is None:
                continue
            a, b = assign.get(parent), assign.get(node)
            if a and b and a != b:
                derived[ei] = b
        if len(derived) >= 2:
            opp += 1
            counts = defaultdict(list)
            for ei, aa in derived.items():
                counts[aa].append(ei)
            for aa, evs in counts.items():
                if len(evs) >= 2:
                    conv_sites.append((c + 1, aa, len(evs), sorted(evs)))
                    break
    return opp, conv_sites


def main():
    rows, hits = [], []
    for g in GENES:
        aln = os.path.join(TRIM, f"{g}.trim.fa")
        tf = os.path.join(TREES, f"{g}.treefile")
        sf = os.path.join(SCEN, f"{g}.scenario")
        if not all(os.path.exists(p) for p in (aln, tf, sf)):
            continue
        seqs = read_fasta(aln)
        tree, by_id = numbered(tf)
        order = [n for n in tree.get_leaf_names() if n in seqs]
        ncol = len(seqs[order[0]])
        events = [[int(x) for x in e.split(",")] for e in open(sf).read().strip().split("/")]

        opp, conv = tally(tree, by_id, events, seqs, ncol, order)
        for site, aa, nev, evs in conv:
            hits.append(dict(gene=g, site=site, residue=aa, n_events=nev,
                             events=";".join(map(str, evs))))

        # Null: re-place the same number of events at random internal nodes of
        # matched clade size, preserving phylogenetic clumping.
        sizes = [len(by_id[e[0]]) for e in events]
        pool = defaultdict(list)
        for n in tree.traverse():
            if n.up is not None:
                pool[len(n)].append(n)
        null_conv, null_opp = [], []
        for _ in range(N_PERM):
            picked, used, ok = [], set(), True
            for s in sorted(sizes, reverse=True):
                cands = [n for n in pool.get(s, ()) if not (set(n.get_leaf_names()) & used)]
                if not cands:
                    ok = False
                    break
                p = rng.choice(cands)
                picked.append([p.ND])
                used |= set(p.get_leaf_names())
            if not ok or len(picked) < 2:
                continue
            o, cv = tally(tree, by_id, picked, seqs, ncol, order)
            null_opp.append(o)
            null_conv.append(len(cv))

        nc = sum(null_conv) / len(null_conv) if null_conv else float("nan")
        no = sum(null_opp) / len(null_opp) if null_opp else float("nan")
        p_count = ((sum(1 for v in null_conv if v >= len(conv)) + 1) / (len(null_conv) + 1)
                   if null_conv else float("nan"))

        # The raw count is confounded. Real transition branches carry far more
        # substitutions than random branches of the same clade size (NR1D1: 58
        # opportunity sites against a null of 9.6), because the ASR places
        # transitions on real, often long branches while the null draws branches
        # matched only on the number of descendants. More substitutions of any
        # kind mechanically produce more same-residue coincidences, so a raw
        # excess says nothing about convergence.
        #
        # Conditioning on opportunity removes it: of the sites that DID change
        # in two or more independent lineages, what fraction landed on the same
        # residue? That ratio is what convergence would inflate.
        obs_rate = len(conv) / opp if opp else float("nan")
        null_rates = [c / o for c, o in zip(null_conv, null_opp) if o > 0]
        nr = sum(null_rates) / len(null_rates) if null_rates else float("nan")
        p_rate = ((sum(1 for v in null_rates if v >= obs_rate) + 1) / (len(null_rates) + 1)
                  if null_rates and opp else float("nan"))

        rows.append(dict(gene=g, n_sites=ncol, n_events=len(events),
                         sites_2plus_events_changed=opp,
                         null_mean_opportunity=round(no, 1),
                         convergent_sites=len(conv),
                         null_mean_convergent=round(nc, 2),
                         p_count_confounded=round(p_count, 4),
                         obs_same_residue_rate=(round(obs_rate, 4) if opp else None),
                         null_same_residue_rate=round(nr, 4),
                         p_rate=(round(p_rate, 4) if opp else None)))
        rate_s = f"{obs_rate:.3f}" if opp else "  n/a"
        p_s = f"{p_rate:.3f}" if opp else "  n/a"
        print(f"{g:10s} opp={opp:4d} (null {no:6.1f})  conv={len(conv):3d} (null {nc:5.2f})  "
              f"rate={rate_s} (null {nr:.3f}, p={p_s})", flush=True)

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "observed_convergent_substitutions.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    if hits:
        with open(os.path.join(OUT, "observed_convergent_sites.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(hits[0].keys()))
            w.writeheader()
            w.writerows(hits)
    # Pooled across genes. Per-gene counts are small, so the aggregate rate is
    # the better-powered comparison, and it is the number that answers the
    # question the whole study asks.
    tot_opp = sum(r["sites_2plus_events_changed"] for r in rows)
    tot_conv = sum(r["convergent_sites"] for r in rows)
    null_opp_tot = sum(r["null_mean_opportunity"] for r in rows)
    null_conv_tot = sum(r["null_mean_convergent"] for r in rows)
    obs_rate = tot_conv / tot_opp if tot_opp else float("nan")
    null_rate = null_conv_tot / null_opp_tot if null_opp_tot else float("nan")
    print(f"\nPOOLED across {len(rows)} genes")
    print(f"  opportunity (sites changing in >=2 independent diurnal lineages): "
          f"{tot_opp}  (null {null_opp_tot:.1f})")
    print(f"  of those, same residue in >=2 lineages: {tot_conv}  (null {null_conv_tot:.1f})")
    print(f"  same-residue RATE: {obs_rate:.4f} observed vs {null_rate:.4f} null "
          f"({'excess' if obs_rate > null_rate else 'no excess'})")


if __name__ == "__main__":
    main()
