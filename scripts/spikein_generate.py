#!/usr/bin/env python3
"""Plant known convergent sites in the real alignments, for an end-to-end control.

Why this exists
---------------
Nine silent defects have been found in this pipeline. Not one raised an error;
every one produced a confident, plausible, wrong number, and several sat in code
that had already been reviewed. They were all interface errors: ape numbering
against ete3 numbering, species-tree clades against gene-tree clades under a
different rooting, node labels in memory against node labels on disk, a CLI flag
order against a tree's root state, alignment columns before pruning against
after.

Fixing components one at a time does not establish that the assembled pipeline
works. `pcoc_sim` does not either: it validates PCOC against PCOC's own model,
using its own scenario handling, so it cannot catch a scenario built wrongly.

This does. It plants convergence in the REAL alignments at known positions, then
the whole pipeline runs unchanged. If the planted sites come back, the plumbing
carries signal from end to end. If they do not, the null means nothing until we
find where the signal dies.

Design
------
Convergence is planted by writing the SAME residue into every species descending
from a chosen set of diurnal transition events, at a chosen column. The residue
is picked to be absent, or nearly absent, at that column so the planted signal is
unambiguous, and columns are chosen to be variable but not saturated.

Sites are planted at three difficulty levels, because the k-curve
(PROGRESS_REPORT.md 0A.5b) showed PCOC only detects near-universal convergence:

  all      all 10 gain events converge      PCOC should find these
  most     7 of 10 converge                 borderline for PCOC
  few      3 of 10 converge                 PCOC should miss, model-free should not

That spread turns the control into a measurement rather than a pass/fail: it
reports which methods recover which difficulty level, and it measures TDG09's
false positive rate directly, since TDG09 currently calls 23 percent of variable
sites significant and nothing corroborates it.

The answer key is written to a separate file. Evaluate with
spikein_evaluate.py, which is the only thing that should read it.

Outputs: results/spikein/trim/<gene>.trim.fa   spiked alignments
         results/spikein/ANSWER_KEY.csv        planted sites (do not peek)
"""
import csv
import os
import random
from collections import Counter

from ete3 import Tree

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TRIM = os.path.join(PROJ, "results", "trim")
TREES = os.path.join(PROJ, "results", "branchlengths", "pcoc")
SCEN = os.path.join(PROJ, "results", "pcoc", "scenarios", "ER_gain")
OUT = os.path.join(PROJ, "results", "spikein")

GENES = ("CLOCK NPAS2 ARNTL PER1 PER2 PER3 CRY1 CRY2 NR1D1 NR1D2 "
         "RORA RORB RORC CSNK1D CSNK1E FBXL3 BHLHE40 BHLHE41").split()

GAPS = set("-X?*BZJU")
AAS = "ACDEFGHIKLMNPQRSTVWY"

# Sites planted per gene at each difficulty level.
LEVELS = {"all": 2, "most": 1, "few": 2}
EVENTS_PER_LEVEL = {"all": None, "most": 7, "few": 3}   # None means every event

SEED = 20260814
rng = random.Random(SEED)


def read_fasta(path):
    seqs, name, buf, order = {}, None, [], []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name = line[1:].split()[0]
            order.append(name)
            buf = []
        elif line:
            buf.append(line)
    if name:
        seqs[name] = "".join(buf)
    return seqs, order


def numbered(tree_file):
    t = Tree(tree_file, format=1)
    by_id = {}
    for i, n in enumerate(t.traverse("postorder")):
        n.add_features(ND=i)
        by_id[i] = n
    return t, by_id


def main():
    os.makedirs(os.path.join(OUT, "trim"), exist_ok=True)
    key = []

    for g in GENES:
        aln = os.path.join(TRIM, f"{g}.trim.fa")
        tf = os.path.join(TREES, f"{g}.treefile")
        sf = os.path.join(SCEN, f"{g}.scenario")
        if not all(os.path.exists(p) for p in (aln, tf, sf)):
            print(f"{g}: missing input, skipped")
            continue

        seqs, order = read_fasta(aln)
        _, by_id = numbered(tf)
        events = [[int(x) for x in e.split(",")] for e in
                  open(sf).read().strip().split("/")]

        # Species descending from each event, i.e. the lineages that "converged".
        event_tips = []
        for ev in events:
            tips = set()
            for nid in ev:
                tips |= set(by_id[nid].get_leaf_names())
            event_tips.append(tips & set(seqs))
        usable = [i for i, t in enumerate(event_tips) if t]

        ncol = len(seqs[order[0]])
        # Columns worth planting in: variable, not saturated, few gaps. Planting
        # into an invariant column would create a signal no method could miss and
        # would not test anything interesting.
        candidates = []
        for c in range(ncol):
            col = [seqs[s][c] for s in order]
            ngap = sum(1 for ch in col if ch in GAPS)
            resid = Counter(ch for ch in col if ch not in GAPS)
            if ngap > 0.1 * len(col) or len(resid) < 2 or len(resid) > 6:
                continue
            candidates.append(c)
        rng.shuffle(candidates)

        chosen = []
        seqs_out = dict(seqs)
        for level, n_sites in LEVELS.items():
            for _ in range(n_sites):
                if not candidates:
                    break
                c = candidates.pop()
                col = [seqs[s][c] for s in order]
                present = {ch for ch in col if ch not in GAPS}
                # A residue absent at this column, so the planted state cannot be
                # confused with standing variation.
                pool = [a for a in AAS if a not in present]
                if not pool:
                    continue
                aa = rng.choice(pool)

                k = EVENTS_PER_LEVEL[level]
                idx = usable if k is None else sorted(rng.sample(usable, min(k, len(usable))))
                targets = set()
                for i in idx:
                    targets |= event_tips[i]
                if len(targets) < 2:
                    continue

                for s in targets:
                    if seqs_out[s][c] in GAPS:
                        continue        # never overwrite a gap: that invents data
                    seqs_out[s] = seqs_out[s][:c] + aa + seqs_out[s][c + 1:]

                chosen.append(c)
                key.append(dict(gene=g, trimmed_col=c + 1, planted_residue=aa,
                                level=level, n_events=len(idx),
                                events=";".join(map(str, idx)),
                                n_species_changed=len(targets)))

        with open(os.path.join(OUT, "trim", f"{g}.trim.fa"), "w") as fh:
            for s in order:
                fh.write(f">{s}\n{seqs_out[s]}\n")
        print(f"{g:9s} {ncol:5d} cols, {len(candidates)} usable, planted {len(chosen)}")

    with open(os.path.join(OUT, "ANSWER_KEY.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(key[0].keys()))
        w.writeheader()
        w.writerows(key)

    lv = Counter(k["level"] for k in key)
    print(f"\nplanted {len(key)} sites across {len({k['gene'] for k in key})} genes")
    print(f"  by level: {dict(lv)}")
    print(f"spiked alignments: {OUT}/trim/")
    print(f"answer key:        {OUT}/ANSWER_KEY.csv   (do not read before evaluating)")


if __name__ == "__main__":
    main()
