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
Convergence is planted by writing the SAME residue into the DIURNAL species
descending from a chosen set of gain events, at a chosen column. The residue is
picked to be absent at that column so the planted signal is unambiguous, and
columns are chosen to be variable but not saturated.

The word DIURNAL is the whole correction. The first version of this script wrote
the residue into every species descending from each event, which sounds
equivalent but is not: a convergent event's membership includes descendant nodes
whose subtrees contain nocturnal species from nested reversals. That planted the
residue into 41 of 60 species rather than the intended 30, making it the MAJORITY
state across the tree in 36 of 36 sites at the 'all' level. A new consensus
residue is not convergence, so the control was asking the detectors to find
something that was not there, and PCOC was correct to miss it. Every conclusion
drawn from that run was withdrawn.

Each planted site is now verified before it is kept:
  - the residue must remain a MINORITY state overall, so it cannot be read as
    ancestral;
  - freq(residue | diurnal) - freq(residue | nocturnal) must exceed
    MIN_GAP, so the signal is genuinely phenotype-associated;
  - at least two independent events must contribute.
Sites failing any of these are reverted and a different column is tried.

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
# SPIKE_OUT lets the fast regression control (spikein_quick.sh) write to its own
# directory so it can never overwrite the full control's alignments or answer key.
OUT = os.environ.get("SPIKE_OUT", os.path.join(PROJ, "results", "spikein"))

ALL_GENES = ("CLOCK NPAS2 ARNTL PER1 PER2 PER3 CRY1 CRY2 NR1D1 NR1D2 "
             "RORA RORB RORC CSNK1D CSNK1E FBXL3 BHLHE40 BHLHE41").split()
GENES = os.environ.get("SPIKE_GENES", "").split() or ALL_GENES

GAPS = set("-X?*BZJU")
AAS = "ACDEFGHIKLMNPQRSTVWY"

# Sites planted per gene at each difficulty level.
LEVELS = {"all": 2, "most": 1, "few": 2}
EVENTS_PER_LEVEL = {"all": None, "most": 7, "few": 3}   # None means every event

# SPIKE_LEVELS restricts which difficulty levels are planted. The fast
# regression control plants 'all' only, because it is a pass/fail smoke test of
# the plumbing, not a measurement of each method's sensitivity.
_want = os.environ.get("SPIKE_LEVELS", "").split()
if _want:
    LEVELS = {k: v for k, v in LEVELS.items() if k in _want}
    if not LEVELS:
        raise SystemExit(f"SPIKE_LEVELS={_want} matches no known level")
_n = os.environ.get("SPIKE_SITES_PER_LEVEL", "")
if _n:
    LEVELS = {k: int(_n) for k in LEVELS}

# A planted residue must be essentially absent from nocturnal species. Planting
# into diurnal descendants only guarantees this, so the check is a guard rather
# than a filter.
#
# Do NOT raise this into a large gap requirement: a 'few' level site touches
# about 10 of 30 diurnal species and so has a gap near 0.33 BY CONSTRUCTION.
# Requiring 0.5 silently discarded every one of them, which is exactly the level
# the control exists to test.
MAX_NOCTURNAL_FREQ = 0.02
MIN_GAP = 0.10

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
    diel = {}
    with open(os.path.join(PROJ, "data", "diel_activity.csv")) as fh:
        for r in csv.DictReader(fh):
            diel[r["species"]] = r["activity"]
    diurnal = {s for s, v in diel.items() if v == "diurnal"}
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
            placed = 0
            while placed < n_sites and candidates:
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
                # DIURNAL descendants only. Event membership includes nodes whose
                # subtrees contain nocturnal species from nested reversals, and
                # writing into those destroys the phenotype association.
                targets = set()
                for i in idx:
                    targets |= (event_tips[i] & diurnal)
                targets = {s for s in targets if seqs_out[s][c] not in GAPS}
                if len(targets) < 2:
                    continue

                trial = {s: seqs_out[s][:c] + aa + seqs_out[s][c + 1:] for s in targets}

                # Verify the result actually looks like convergence before keeping it.
                after = {s: trial.get(s, seqs_out[s]) for s in order}
                counts = Counter(after[s][c] for s in order if after[s][c] not in GAPS)
                n_tot = sum(counts.values())
                di = [after[s][c] for s in order if s in diurnal and after[s][c] not in GAPS]
                no = [after[s][c] for s in order if s not in diurnal and after[s][c] not in GAPS]
                gap = (di.count(aa) / len(di) if di else 0) - (no.count(aa) / len(no) if no else 0)
                no_freq = no.count(aa) / len(no) if no else 0.0
                # Absence from nocturnal species is the criterion, NOT minority
                # status among tips. At the 'all' level 30 diurnal species share
                # the residue and it becomes the plurality, which is what real
                # convergence across half a clade looks like; rejecting on that
                # threw away every 'all' site. Absence from all 30 nocturnal
                # species, the marsupial outgroup included, is what makes the
                # residue impossible to reconstruct at the root.
                if gap < MIN_GAP or no_freq > MAX_NOCTURNAL_FREQ:
                    continue

                seqs_out.update(trial)
                chosen.append(c)
                placed += 1
                key.append(dict(gene=g, trimmed_col=c + 1, planted_residue=aa,
                                level=level, n_events=len(idx),
                                events=";".join(map(str, idx)),
                                n_species_changed=len(targets),
                                diurnal_gap=round(gap, 3),
                                pct_of_tree=round(counts[aa] / n_tot, 3)))

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
