#!/usr/bin/env python3
"""How much evolutionary opportunity was there for convergence to arise at all?

This is the necessary companion to the pcoc_sim power calibration, and it asks a
different question.

pcoc_sim simulates the convergent shift using Bio++'s OneChange model on every
transition branch (bpp_lib.py: MODEL_OC='OneChange(model=$(MODEL_C))'), which
CONDITIONS on at least one substitution occurring there. The simulated change is
therefore guaranteed to happen regardless of branch length. That is a fair test
of the detector, and pcoc_det fits the same OneChange model, so simulation and
detection are matched. But it means the reported power is CONDITIONAL: given
that a convergent substitution occurred in each lineage, PCOC finds it.

It says nothing about whether there was time for those substitutions to occur.
ARNTL scores power 1.000 while its 52 convergent branches total only 0.097
substitutions per site, so at any single site a substitution on a given
transition branch is a ~0.2 percent event. Convergence requires the same site to
change independently in several lineages, so for such a gene the null is close to
guaranteed by branch length alone, whatever the biology.

This script quantifies that per gene. For each convergent event it computes

    p_event = 1 - exp(-L)

the probability of at least one substitution at a given site somewhere in that
event, under a Poisson process with rate 1 per unit branch length. Summing over
events gives the expected number of independent lineages that could show a
change at a site. Convergence needs at least two, so a gene whose expected count
is below 2 has essentially no opportunity and its null is uninformative.

Two versions of L are reported:
  transition  only the transition branch itself, which is what PCOC's OneChange
              component keys on;
  clade       the whole convergent event (transition branch plus the descendant
              branches that stay in the derived state), the upper bound on where
              a derived profile could reveal itself.

Output: results/pcoc_sim/transition_opportunity.csv
"""
import csv
import math
import os

from ete3 import Tree

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TREES = os.path.join(PROJ, "results", "branchlengths", "pcoc")
SCEN = os.path.join(PROJ, "results", "pcoc", "scenarios", "ER_gain")
OUT = os.path.join(PROJ, "results", "pcoc_sim")

GENES = ("CLOCK NPAS2 ARNTL PER1 PER2 PER3 CRY1 CRY2 NR1D1 NR1D2 "
         "RORA RORB RORC CSNK1D CSNK1E FBXL3 BHLHE40 BHLHE41").split()

# Convergence needs the same site to change in at least this many independent
# lineages, so it is the minimum for the event to be definable at all.
MIN_LINEAGES = 2

# A gene is judged to have adequate opportunity on the expected number of SITES
# reaching MIN_LINEAGES, not on the expected lineages per site. The per-site
# expectation is below 2 in almost every gene, which sounds damning but is not:
# it is an average over hundreds of columns, and a gene with 600 sites and a
# per-site probability of 0.05 still expects 30 usable sites. Judging on the
# per-site average alone would wrongly condemn nearly the whole dataset.
#
# This whole script is the branch-length ESTIMATE of opportunity. It assumes a
# uniform per-site rate and so overestimates: observed_convergent_substitutions.py
# counts the same quantity directly by parsimony with no rate assumption, and
# gets substantially lower numbers (PER2: 196 observed against 719 predicted
# here). Prefer the empirical count; this one shows which genes are constrained
# by branch length alone.
MIN_SITES = 10


def numbered(tree_file):
    """Reproduce PCOC's node numbering: postorder traversal from 0.

    Matches events_placing.py:init_tree, so the integers in the scenario files
    index the same branches here as they do inside PCOC.
    """
    t = Tree(tree_file, format=1)
    by_id = {}
    for i, n in enumerate(t.traverse("postorder")):
        n.add_features(ND=i)
        by_id[i] = n
    return t, by_id


def prob_at_least(probs, k):
    """Poisson-binomial tail: P(at least k of these independent events occur).

    The events have different branch lengths and therefore different
    probabilities, so a plain binomial would be wrong.
    """
    dist = [1.0]
    for p in probs:
        nxt = [0.0] * (len(dist) + 1)
        for i, d in enumerate(dist):
            nxt[i] += d * (1 - p)
            nxt[i + 1] += d * p
        dist = nxt
    return sum(dist[k:])


def n_sites(gene):
    """Column count of the trimmed alignment actually fed to PCOC."""
    path = os.path.join(PROJ, "results", "trim", f"{gene}.trim.fa")
    if not os.path.exists(path):
        return 0
    n, started = 0, False
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if started:
                break
            started = True
        elif started:
            n += len(line)
    return n


def main():
    rows = []
    for g in GENES:
        tf = os.path.join(TREES, f"{g}.treefile")
        sf = os.path.join(SCEN, f"{g}.scenario")
        if not (os.path.exists(tf) and os.path.exists(sf)):
            continue
        _, by_id = numbered(tf)
        events = [e.split(",") for e in open(sf).read().strip().split("/")]

        p_trans, p_clade, l_trans_all, l_clade_all = [], [], 0.0, 0.0
        for ev in events:
            ids = [int(x) for x in ev]
            lt = by_id[ids[0]].dist
            lc = sum(by_id[i].dist for i in ids)
            l_trans_all += lt
            l_clade_all += lc
            p_trans.append(1 - math.exp(-lt))
            p_clade.append(1 - math.exp(-lc))

        exp_trans = sum(p_trans)
        exp_clade = sum(p_clade)
        # Expected number of SITES, not lineages, that could look convergent by
        # chance. The events have different branch lengths, so the count of
        # events changing at a site is Poisson-binomial rather than binomial.
        p2 = prob_at_least(p_clade, MIN_LINEAGES)
        nsites = n_sites(g)
        rows.append(dict(
            gene=g, n_events=len(events), n_sites=nsites,
            total_transition_bl=round(l_trans_all, 4),
            total_clade_bl=round(l_clade_all, 4),
            exp_lineages_changed_transition=round(exp_trans, 3),
            exp_lineages_changed_clade=round(exp_clade, 3),
            max_single_event_p=round(max(p_clade), 3),
            p_site_2plus_lineages=round(p2, 5),
            exp_sites_2plus_lineages=(round(p2 * nsites, 1) if nsites else None),
            opportunity=("adequate" if p2 * nsites >= MIN_SITES else "LOW"),
        ))

    rows.sort(key=lambda r: r["exp_lineages_changed_clade"])
    hdr = (f"{'gene':10s} {'sites':>6s} {'cladeBL':>8s} {'E[lin|clade]':>12s} "
           f"{'P(site>=2)':>11s} {'E[sites>=2]':>12s} {'verdict':>9s}")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        print(f"{r['gene']:10s} {r['n_sites']:6d} {r['total_clade_bl']:8.4f} "
              f"{r['exp_lineages_changed_clade']:12.3f} {r['p_site_2plus_lineages']:11.5f} "
              f"{r['exp_sites_2plus_lineages']:12.1f} {r['opportunity']:>9s}")

    os.makedirs(OUT, exist_ok=True)
    out = os.path.join(OUT, "transition_opportunity.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    low = [r["gene"] for r in rows if r["opportunity"] == "LOW"]
    total = sum(r["exp_sites_2plus_lineages"] or 0 for r in rows)
    print(f"\n{len(low)} of {len(rows)} genes have LOW opportunity "
          f"(fewer than {MIN_SITES} sites expected to change in >= {MIN_LINEAGES} "
          f"lineages): {', '.join(low) if low else 'none'}")
    print(f"Expected sites with opportunity, summed over genes: {total:.0f} "
          f"(branch-length estimate; see observed_convergent_substitutions.py "
          f"for the direct parsimony count)")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
