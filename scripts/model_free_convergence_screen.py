#!/usr/bin/env python3
"""Model-free screen for diurnal-specific residues, as a cross-check on PCOC.

The ER_gain PCOC sweep returned zero convergent sites. That could mean there is
no convergence, or that PCOC's particular model missed it. This screen asks the
same question without any evolutionary model: is there any alignment column
where the diurnal species share an amino acid that the nocturnal species lack?

For each column the statistic is

    score = max over residues a of  ( freq(a | diurnal) - freq(a | nocturnal) )

which is 1.0 for a perfectly phenotype-diagnostic column and near 0 for a column
whose residues ignore the phenotype.

The null is NOT a label shuffle. Diurnality arose in 10 clustered events, so
shuffling tip labels destroys the phylogenetic clumping and would make almost
any clade-restricted residue look significant. Instead the null re-places the
observed diurnal clades at random points on the same gene tree, drawing clades
of the same sizes, so the null carries the same autocorrelation structure and
only the PHENOTYPE ASSIGNMENT is randomized.

Inputs : results/trim/<gene>.trim.fa, results/branchlengths/pcoc/<gene>.treefile,
         data/diel_activity.csv
Outputs: results/pcoc_sim/model_free_screen.csv  (per-site scores above cutoff)
         results/pcoc_sim/model_free_summary.csv (per-gene observed vs null)
"""
import csv
import os
import random

import numpy as np
from ete3 import Tree

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TRIM = os.environ.get("TRIM_DIR", os.path.join(PROJ, "results", "trim"))
TREES = os.path.join(PROJ, "results", "branchlengths", "pcoc")
OUT = os.environ.get("OUT_DIR", os.path.join(PROJ, "results", "pcoc_sim"))

GENES = ("CLOCK NPAS2 ARNTL PER1 PER2 PER3 CRY1 CRY2 NR1D1 NR1D2 "
         "RORA RORB RORC CSNK1D CSNK1E FBXL3 BHLHE40 BHLHE41").split()

AAS = "ACDEFGHIKLMNPQRSTVWY"
N_PERM = 1000
MIN_PER_GROUP = 5     # columns with fewer ungapped residues in a group are skipped
# Report any column at least this diagnostic. 0.6 means the top residue is at
# least 60 percentage points commoner in diurnal than nocturnal species, e.g.
# present in 20 of 30 diurnal and 2 of 30 nocturnal.
REPORT_CUTOFF = 0.6
rng = random.Random(20260813)


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
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


def encode(seqs, order):
    """One-hot encode the alignment as (nseq, ncol, 20) uint8 over the 20 AAs.

    Anything not in AAS (gaps, X, ambiguity) encodes as all-zero, so it drops out
    of the per-column counts instead of being treated as a residue.
    """
    ncol = len(seqs[order[0]])
    arr = np.zeros((len(order), ncol, len(AAS)), dtype=np.uint8)
    idx = {a: i for i, a in enumerate(AAS)}
    for r, name in enumerate(order):
        for c, ch in enumerate(seqs[name]):
            j = idx.get(ch)
            if j is not None:
                arr[r, c, j] = 1
    return arr


def scores(onehot, total, mask):
    """Vectorised per-column max frequency gap between mask and its complement.

    Returns (score, residue_index) arrays of length ncol. Only counts_a is
    computed; the complement falls out of the precomputed column totals.
    """
    ca = onehot[mask].sum(axis=0).astype(np.float64)      # (ncol, 20)
    cb = total - ca
    na = ca.sum(axis=1)
    nb = cb.sum(axis=1)
    ok = (na >= MIN_PER_GROUP) & (nb >= MIN_PER_GROUP)
    with np.errstate(invalid="ignore", divide="ignore"):
        diff = ca / na[:, None] - cb / nb[:, None]
    diff[~ok] = 0.0
    diff = np.nan_to_num(diff, nan=0.0)
    best = diff.argmax(axis=1)
    return diff[np.arange(diff.shape[0]), best], best, ca, cb


def maximal_diurnal_clades(tree, lab):
    """Sizes of the maximal all-diurnal clades, and the tip sets behind them."""
    sizes = []
    covered = set()
    for n in tree.traverse("levelorder"):
        leaves = set(n.get_leaf_names()) & set(lab)
        if not leaves or leaves & covered:
            continue
        if all(lab[t] == "diurnal" for t in leaves):
            sizes.append(len(leaves))
            covered |= leaves
    return sizes


def random_masks(tree, sizes, order, n_perm):
    """Draw phenotype masks by re-placing clades of the observed sizes at random.

    Node leaf sets are precomputed once; ete3 traversal inside the permutation
    loop is what made the first version of this script unusable.
    """
    pos = {name: i for i, name in enumerate(order)}
    by_size = {}
    for n in tree.traverse():
        if n.is_root():
            continue
        leaves = [pos[x] for x in n.get_leaf_names() if x in pos]
        if leaves:
            by_size.setdefault(len(leaves), []).append(frozenset(leaves))
    out = []
    for _ in range(n_perm):
        for _attempt in range(50):
            chosen, used, ok = set(), set(), True
            for s in sorted(sizes, reverse=True):
                cands = [c for c in by_size.get(s, ()) if not (c & used)]
                if not cands:
                    ok = False
                    break
                pick = rng.choice(cands)
                chosen |= pick
                used |= pick
            if ok and chosen:
                m = np.zeros(len(order), dtype=bool)
                m[list(chosen)] = True
                out.append(m)
                break
    return out


def main():
    labels = {}
    with open(os.path.join(PROJ, "data", "diel_activity.csv")) as fh:
        for r in csv.DictReader(fh):
            labels[r["species"]] = r["activity"]

    hits, summary = [], []
    for g in GENES:
        aln = os.path.join(TRIM, f"{g}.trim.fa")
        twk = os.path.join(TREES, f"{g}.treefile")
        if not (os.path.exists(aln) and os.path.exists(twk)):
            print(f"skip {g}: missing input")
            continue
        seqs = read_fasta(aln)
        tree = Tree(twk, format=1)
        order = [t for t in tree.get_leaf_names() if t in seqs and t in labels]
        lab = {t: labels[t] for t in order}

        onehot = encode(seqs, order)
        total = onehot.sum(axis=0).astype(np.float64)
        obs_mask = np.array([lab[t] == "diurnal" for t in order])

        obs, best, ca, cb = scores(onehot, total, obs_mask)
        obs_max = float(obs.max()) if obs.size else 0.0
        over = np.where(obs >= REPORT_CUTOFF)[0]
        for i in over:
            hits.append(dict(gene=g, site=int(i) + 1, score=round(float(obs[i]), 3),
                             residue=AAS[best[i]],
                             n_diurnal=int(ca[i, best[i]]),
                             n_nocturnal=int(cb[i, best[i]])))

        sizes = maximal_diurnal_clades(tree, lab)
        masks = random_masks(tree, sizes, order, N_PERM)
        null_max, null_cnt = [], []
        for m in masks:
            s, _, _, _ = scores(onehot, total, m)
            null_max.append(float(s.max()) if s.size else 0.0)
            null_cnt.append(int((s >= REPORT_CUTOFF).sum()))

        if null_max:
            p_max = (sum(1 for v in null_max if v >= obs_max) + 1) / (len(null_max) + 1)
            p_cnt = (sum(1 for v in null_cnt if v >= len(over)) + 1) / (len(null_cnt) + 1)
            nm, nc = float(np.mean(null_max)), float(np.mean(null_cnt))
        else:
            p_max = p_cnt = nm = nc = float("nan")

        summary.append(dict(gene=g, n_sites=int(obs.size), n_diurnal_clades=len(sizes),
                            n_perm=len(masks), obs_max_score=round(obs_max, 3),
                            null_mean_max=round(nm, 3), p_max=round(p_max, 4),
                            obs_n_over_cutoff=int(len(over)),
                            null_mean_n_over=round(nc, 2), p_count=round(p_cnt, 4)))
        print(f"{g:10s} sites={obs.size:5d} max={obs_max:.3f} (null {nm:.3f}, p={p_max:.3f})  "
              f"n>={REPORT_CUTOFF}: {len(over)} (null {nc:.2f}, p={p_cnt:.3f})", flush=True)

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "model_free_summary.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary[0].keys()))
        w.writeheader()
        w.writerows(summary)
    if hits:
        with open(os.path.join(OUT, "model_free_screen.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(hits[0].keys()))
            w.writeheader()
            w.writerows(hits)
    print(f"\n{len(hits)} sites at or above score {REPORT_CUTOFF}; wrote summaries to {OUT}")


if __name__ == "__main__":
    main()
