#!/usr/bin/env python3
"""Count the convergent substitutions that ACTUALLY occurred, by parsimony.

This is the study's primary evidence. PCOC only has power when nearly all 10
diurnal lineages converge together (see PROGRESS_REPORT.md 0A.5b), so a signal
confined to a few lineages would be invisible to it. This test declares no
convergent set at all: it counts substitutions on every candidate transition
branch and asks how often independent lineages landed on the SAME residue.

Fitch parsimony assigns an amino acid to every internal node per column, so a
substitution on a branch is simply parent state != child state. No evolutionary
model, no rate assumption.

Four things this version does that a first pass should not be trusted without
-----------------------------------------------------------------------------
1. **Conditions on opportunity.** Raw counts are confounded: real transition
   branches carry far more substitutions of ANY kind than random branches of
   matched clade size (NR1D1: 58 opportunity sites against a null of 9.6),
   because the ASR places transitions on real, often long branches. More
   substitutions mechanically produce more same-residue coincidences. The test
   is therefore the FRACTION of sites changing in 2+ lineages that landed on the
   same residue.

2. **Matches the null on branch length, not only clade size.** Conditioning on
   opportunity corrects the symptom; drawing length-matched null branches
   removes the cause. Both are reported, and they should agree.

3. **Propagates ancestral-state uncertainty.** Fitch leaves ties, and resolving
   them deterministically (the previous behaviour) silently commits to one
   history. Each site is resolved N_RESOLUTIONS times with random tie-breaking,
   so the reported statistic carries the spread across equally parsimonious
   reconstructions.

4. **Tests individual sites, not just the aggregate.** The pooled rate is well
   powered against a pervasive signal but blind to a handful of real sites,
   which is the more likely biology. Every site gets a permutation p-value and
   Benjamini-Hochberg q-value across all genes.

5. **Requires the convergent residue to be phenotype-specific.** Points 1 to 4
   alone are not enough, and the first run of this script proved it: 12 sites
   reached q <= 0.05, and on inspection 9 of them carried a residue that is just
   as common in nocturnal species as in diurnal ones. CLOCK site 675 was the
   clearest, residue V in 2 diurnal and 2 nocturnal species, a frequency gap of
   -0.002, two of the four carriers being nocturnal marsupials.

   The reason is that the permutation test asks only whether 2+ TRANSITION
   branches changed to the same residue more often than randomly placed branches
   would. At a tolerant site a residue arises repeatedly all over the tree, some
   of those changes land on transition branches by chance, and because the real
   transitions are a phylogenetically clustered set while the null scatters
   branches everywhere, the coincidence looks unusual. That is homoplasy
   unrelated to diel activity, not convergence.

   Every candidate site therefore also carries `diurnal_gap`, the frequency of
   the convergent residue among diurnal species minus among nocturnal species. A
   site is only reported as convergent if it is BOTH statistically unusual and
   phenotype-specific.

Speed note: the parsimony reconstruction does not depend on which branches are
called transitions, so it is computed ONCE per site and the permutations become
lookups. That is what makes 2000 permutations affordable.

Outputs: results/pcoc_sim/observed_convergent_substitutions.csv  per gene
         results/pcoc_sim/observed_convergent_sites.csv          per site, with q
"""
import csv
import os
import random
from collections import defaultdict

from ete3 import Tree

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TRIM = os.environ.get("TRIM_DIR", os.path.join(PROJ, "results", "trim"))
TREES = os.path.join(PROJ, "results", "branchlengths", "pcoc")
SCEN = os.path.join(PROJ, "results", "pcoc", "scenarios", "ER_gain")
OUT = os.environ.get("OUT_DIR", os.path.join(PROJ, "results", "pcoc_sim"))

GENES = ("CLOCK NPAS2 ARNTL PER1 PER2 PER3 CRY1 CRY2 NR1D1 NR1D2 "
         "RORA RORB RORC CSNK1D CSNK1E FBXL3 BHLHE40 BHLHE41").split()

GAPS = set("-X?*BZJU")
N_PERM = 2000
N_RESOLUTIONS = 10        # equally parsimonious reconstructions per site
BL_TOLERANCE = 2.0        # null branch must be within this factor of the real one
FDR = 0.05
# Minimum enrichment of the convergent residue in diurnal over nocturnal species.
# Phenotype-associated convergence must leave the residue commoner in the
# lineages that have the phenotype; without this, homoplasy at tolerant sites is
# indistinguishable from convergence (see point 5 in the docstring).
MIN_DIURNAL_GAP = 0.25
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


def fitch_resolutions(tree, col_state, n_res, order_nodes, parent_id):
    """N_RESOLUTIONS equally parsimonious assignments for one column.

    Returns a list of lists indexed by node id, each holding a residue or None.
    Gaps and ambiguous characters are missing data: they contribute nothing to
    the state sets rather than acting as a 21st residue, so a gapped clade
    cannot manufacture a spurious substitution.

    Ties are broken at random rather than deterministically. A single
    deterministic resolution commits to one history among several equally
    supported ones, which understates uncertainty in exactly the statistic this
    script reports.
    """
    up = {}
    for n in tree.traverse("postorder"):
        if n.is_leaf():
            ch = col_state.get(n.name)
            up[n.ND] = {ch} if ch and ch not in GAPS else set()
        else:
            sets = [up[c.ND] for c in n.children if up[c.ND]]
            if not sets:
                up[n.ND] = set()
            else:
                inter = set.intersection(*sets)
                up[n.ND] = inter if inter else set().union(*sets)

    out = []
    ambiguous = any(len(s) > 1 for s in up.values())
    reps = n_res if ambiguous else 1
    for _ in range(reps):
        assign = [None] * (max(up) + 1)
        for n in order_nodes:                    # preorder
            s = up[n]
            if not s:
                continue
            p = parent_id[n]
            if p is not None and assign[p] in s:
                assign[n] = assign[p]
            elif len(s) == 1:
                assign[n] = next(iter(s))
            else:
                assign[n] = rng.choice(sorted(s))
        out.append(assign)
    return out


def tally(resolutions, trans_ids, parent_id):
    """Per site: (had opportunity, converged, residue, events) for one event set.

    Averaged over the parsimony resolutions: a site counts as convergent if the
    majority of equally parsimonious histories make it so.
    """
    opp_votes = conv_votes = 0
    best = None
    for assign in resolutions:
        derived = {}
        for ei, nid in enumerate(trans_ids):
            p = parent_id[nid]
            if p is None:
                continue
            a, b = assign[p], assign[nid]
            if a and b and a != b:
                derived[ei] = b
        if len(derived) >= 2:
            opp_votes += 1
            counts = defaultdict(list)
            for ei, aa in derived.items():
                counts[aa].append(ei)
            hit = max(counts.items(), key=lambda kv: len(kv[1]))
            if len(hit[1]) >= 2:
                conv_votes += 1
                if best is None or len(hit[1]) > best[1]:
                    best = (hit[0], len(hit[1]), sorted(hit[1]))
    n = len(resolutions)
    return opp_votes > n / 2, conv_votes > n / 2, best


def bh(pvals):
    n = len(pvals)
    if not n:
        return []
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        k = n - rank + 1
        prev = min(prev, pvals[i] * n / k)
        q[i] = prev
    return q


def draw_null_events(pool_by_size, real, rng_):
    """Random transition set matched on BOTH clade size and branch length.

    Matching only on clade size lets the null draw short branches where the real
    transitions are long, which is precisely the confound that made the raw
    counts look significant. Tolerance widens only if nothing matches, and the
    caller is told how often that happened.
    """
    picked, used, relaxed = [], set(), 0
    # Largest clades first. Placed in event order instead, the big clades find no
    # unused node of their size left and the whole draw is discarded; PER1 and
    # PER2 lost every one of 2000 permutations that way.
    for size, bl, tips in sorted(real, key=lambda r: -r[0]):
        for factor in (BL_TOLERANCE, BL_TOLERANCE ** 2, None):
            cands = [n for n in pool_by_size.get(size, ())
                     if not (n["tips"] & used)
                     and (factor is None
                          or (bl / factor <= n["bl"] <= bl * factor)
                          or (bl == 0 and n["bl"] == 0))]
            if cands:
                if factor is not None and factor != BL_TOLERANCE:
                    relaxed += 1
                elif factor is None:
                    relaxed += 1
                c = rng_.choice(cands)
                picked.append(c["nd"])
                used |= c["tips"]
                break
        else:
            return None, relaxed
    return picked, relaxed


def diel_states():
    with open(os.path.join(PROJ, "data", "diel_activity.csv")) as fh:
        return {r["species"]: r["activity"] for r in csv.DictReader(fh)}


def residue_gap(seqs, col, residue, diel):
    """freq(residue | diurnal) - freq(residue | nocturnal), gaps excluded.

    Positive means the residue is commoner in the lineages that have the
    phenotype, which is the minimum any claim of phenotype-associated
    convergence has to satisfy.
    """
    d = n = dt = nt = 0
    for sp, s in seqs.items():
        ch = s[col]
        if ch in GAPS:
            continue
        if diel.get(sp) == "diurnal":
            dt += 1
            d += ch == residue
        else:
            nt += 1
            n += ch == residue
    return (d / dt if dt else 0.0) - (n / nt if nt else 0.0)


def main():
    diel = diel_states()
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
        events = [[int(x) for x in e.split(",")] for e in
                  open(sf).read().strip().split("/")]
        trans_ids = [e[0] for e in events]

        nnode = max(by_id) + 1
        parent_id = [None] * nnode
        for nid, node in by_id.items():
            parent_id[nid] = node.up.ND if node.up is not None else None
        preorder = [n.ND for n in tree.traverse("preorder")]

        # candidate pool for the null, with branch length and tip set
        pool_by_size = defaultdict(list)
        for nid, node in by_id.items():
            if node.up is None:
                continue
            pool_by_size[len(node)].append(
                dict(nd=nid, bl=node.dist, tips=frozenset(node.get_leaf_names())))
        real = [(len(by_id[i]), by_id[i].dist,
                 frozenset(by_id[i].get_leaf_names())) for i in trans_ids]

        # ---- reconstruct every site once, then reuse ------------------------
        per_site = []
        for c in range(ncol):
            col = {nm: seqs[nm][c] for nm in order}
            per_site.append(fitch_resolutions(tree, col, N_RESOLUTIONS,
                                              preorder, parent_id))

        obs_opp, obs_conv, site_hit = 0, 0, {}
        for c, res in enumerate(per_site):
            o, cv, best = tally(res, trans_ids, parent_id)
            if o:
                obs_opp += 1
            if cv and best:
                obs_conv += 1
                site_hit[c] = best

        # ---- permutations: reuse the reconstructions -------------------------
        null_opp, null_conv, relaxed_total = [], [], 0
        site_null = defaultdict(int)
        done = 0
        for _ in range(N_PERM):
            picked, relaxed = draw_null_events(pool_by_size, real, rng)
            if picked is None or len(picked) < 2:
                continue
            relaxed_total += relaxed
            done += 1
            o_n = c_n = 0
            for c, res in enumerate(per_site):
                o, cv, best = tally(res, picked, parent_id)
                if o:
                    o_n += 1
                if cv and best:
                    c_n += 1
                    # per-site null: did THIS site converge as strongly?
                    if c in site_hit and best[1] >= site_hit[c][1]:
                        site_null[c] += 1
            null_opp.append(o_n)
            null_conv.append(c_n)

        nc = sum(null_conv) / len(null_conv) if null_conv else float("nan")
        no = sum(null_opp) / len(null_opp) if null_opp else float("nan")
        obs_rate = obs_conv / obs_opp if obs_opp else float("nan")
        null_rates = [c / o for c, o in zip(null_conv, null_opp) if o > 0]
        nr = sum(null_rates) / len(null_rates) if null_rates else float("nan")
        p_rate = ((sum(1 for v in null_rates if v >= obs_rate) + 1) / (len(null_rates) + 1)
                  if null_rates and obs_opp else float("nan"))

        for c, (aa, nev, evs) in site_hit.items():
            gap = residue_gap(seqs, c, aa, diel)
            hits.append(dict(gene=g, site=c + 1, residue=aa, n_events=nev,
                             events=";".join(map(str, evs)),
                             diurnal_gap=round(gap, 3),
                             phenotype_specific=gap >= MIN_DIURNAL_GAP,
                             p_site=(site_null[c] + 1) / (done + 1)))

        rows.append(dict(gene=g, n_sites=ncol, n_events=len(events),
                         n_perm=done, n_resolutions=N_RESOLUTIONS,
                         null_draws_needing_relaxed_bl=relaxed_total,
                         sites_2plus_events_changed=obs_opp,
                         null_mean_opportunity=round(no, 1),
                         convergent_sites=obs_conv,
                         null_mean_convergent=round(nc, 2),
                         obs_same_residue_rate=(round(obs_rate, 4) if obs_opp else None),
                         null_same_residue_rate=round(nr, 4),
                         p_rate=(round(p_rate, 4) if obs_opp else None)))
        rate_s = f"{obs_rate:.3f}" if obs_opp else "  n/a"
        p_s = f"{p_rate:.3f}" if obs_opp else "  n/a"
        print(f"{g:10s} opp={obs_opp:4d} (null {no:6.1f})  conv={obs_conv:3d} "
              f"(null {nc:5.2f})  rate={rate_s} (null {nr:.3f}, p={p_s})  "
              f"perms={done}", flush=True)

    # ---- per-site FDR across every gene ------------------------------------
    if hits:
        qs = bh([h["p_site"] for h in hits])
        for h, q in zip(hits, qs):
            h["q_site"] = round(q, 4)
            h["p_site"] = round(h["p_site"], 5)
        hits.sort(key=lambda h: h["p_site"])

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

    tot_opp = sum(r["sites_2plus_events_changed"] for r in rows)
    tot_conv = sum(r["convergent_sites"] for r in rows)
    no_t = sum(r["null_mean_opportunity"] for r in rows)
    nc_t = sum(r["null_mean_convergent"] for r in rows)
    print(f"\nPOOLED across {len(rows)} genes")
    print(f"  opportunity (sites changing in >=2 independent diurnal lineages): "
          f"{tot_opp}  (null {no_t:.1f})")
    print(f"  of those, same residue in >=2 lineages: {tot_conv}  (null {nc_t:.1f})")
    if tot_opp and no_t:
        o_r, n_r = tot_conv / tot_opp, nc_t / no_t
        print(f"  same-residue RATE: {o_r:.4f} observed vs {n_r:.4f} null "
              f"({'EXCESS' if o_r > n_r else 'no excess'})")
    if hits:
        sig = [h for h in hits if h["q_site"] <= FDR]
        both = [h for h in sig if h["phenotype_specific"]]
        print(f"\n  per-site test: {len(hits)} candidate sites")
        print(f"    statistically unusual (q <= {FDR}):            {len(sig)}")
        print(f"    ALSO phenotype-specific (gap >= {MIN_DIURNAL_GAP}): {len(both)}")
        print(f"\n  {'gene':9s} {'site':>5s} {'res':>4s} {'events':>7s} {'gap':>7s} "
              f"{'q':>7s} {'phenotype-specific':>19s}")
        for h in sig[:15]:
            print(f"  {h['gene']:9s} {h['site']:5d} {h['residue']:>4s} "
                  f"{h['n_events']:7d} {h['diurnal_gap']:+7.3f} {h['q_site']:7.3f} "
                  f"{str(h['phenotype_specific']):>19s}")
        if not both:
            print("\n  No site is both statistically unusual and phenotype-specific.")
            print("  Sites passing the permutation test alone are homoplasy at tolerant")
            print("  positions, not convergence: their residue is as common in nocturnal")
            print("  species as in diurnal ones.")


if __name__ == "__main__":
    main()
