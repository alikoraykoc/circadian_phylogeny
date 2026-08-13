#!/usr/bin/env python3
"""How many diurnal lineages must share a change before our analysis can see it?

The 18-gene power calibration showed PCOC detects COMPLETE convergence with
certainty: when all 10 gains of diurnality carry the derived profile, every
simulated site is recovered. But the real sweep DECLARES all 10 events convergent
regardless of the truth, and a first 5-point test found that when only 2 of the
10 lineages actually converged, detection collapsed to zero, while the same data
detected with the correct 2-event scenario scored a perfect 1.000. Blindness, not
weakness, and caused by the declaration rather than by lack of signal.

This sweep maps the whole curve: k = 2..10 converging lineages, several
independent draws of WHICH lineages at each k, always detected with the full
declared 10-event scenario, exactly as the real analysis did.

The trick that makes it affordable
----------------------------------
pcoc_det's runtime is dominated by fixed per-tree setup, not by column count, and
the declared scenario is IDENTICAL for every k (that is the whole point). So all
the simulated alignments are concatenated into one and detected in a single run,
then the posteriors are split back apart by block. That turns ~27 detection runs
into 1. PCOC treats sites independently, so concatenation changes nothing.

Group count is not the only variable
------------------------------------
A draw of 3 small clades can carry less evolutionary material than a draw of 2
large ones (an earlier run drew 3 events spanning 9 branches against 2 events
spanning 25). Every replicate therefore records its branch count and total branch
length, so power can be read against both the number of lineages and the amount
of sequence change available.

Usage: python scripts/pcoc_partial_convergence_sweep.py [GENE]   (default CLOCK)
Output: results/pcoc_sim/sweep/<gene>/power_by_k.csv     per-replicate results
        results/pcoc_sim/sweep/<gene>/blocks.csv         the block index
        results/figures/partial_convergence_<gene>.pdf
"""
import csv
import glob
import itertools
import os
import random
import shutil
import subprocess
import sys

from ete3 import Tree

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IMG = "carinerey/pcoc"
GENE = sys.argv[1] if len(sys.argv) > 1 else "CLOCK"

K_VALUES = list(range(2, 11))
N_REPLICATES = 3        # independent draws of WHICH lineages converge, per k
N_SITES = 100
N_COUPLES = 2           # ancestral/convergent profile pairs per replicate
CPU = 4
SEED = 20260813

OUT = os.path.join(PROJ, "results", "pcoc_sim", "sweep", GENE)
SCEN_FILE = os.path.join(PROJ, "results", "pcoc", "scenarios", "ER_gain", f"{GENE}.scenario")
TREE_DIR = os.path.join(PROJ, "results", "pcoc_sim", "trees", GENE)
TREE_FILE = os.path.join(PROJ, "results", "branchlengths", "pcoc", f"{GENE}.treefile")


def docker(args):
    return subprocess.run(["docker", "run", "--rm", "-v", f"{PROJ}:/proj", IMG] + args,
                          capture_output=True, text=True)


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


def numbered_tree():
    """PCOC numbering: postorder from 0 (events_placing.py:init_tree)."""
    t = Tree(TREE_FILE, format=1)
    by_id = {}
    for i, n in enumerate(t.traverse("postorder")):
        n.add_features(ND=i)
        by_id[i] = n
    return by_id


def is_null_alignment(path):
    """A<x>_C<x> is the no-change control; only A<x>_C<y> carries the signal."""
    parts = os.path.basename(path).replace(".fa", "").split("_")
    a = [p for p in parts if p.startswith("A")][-1][1:]
    c = [p for p in parts if p.startswith("C")][-1][1:]
    return a == c


def posteriors(tsv):
    """PCOC_V1 per site; blank cells underflowed to 0, i.e. strongly rejected."""
    out = []
    with open(tsv) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        col = next(c for c in rd.fieldnames if c.startswith("PCOC"))
        for r in rd:
            v = (r.get(col) or "").strip()
            out.append(float(v) if v else 0.0)
    return out


def main():
    full = open(SCEN_FILE).read().strip()
    events = full.split("/")
    n_full = len(events)
    by_id = numbered_tree()
    print(f"{GENE}: {n_full} events, "
          f"{sum(len(e.split(',')) for e in events)} convergent branches")

    os.makedirs(OUT, exist_ok=True)
    os.makedirs(TREE_DIR, exist_ok=True)
    shutil.copy(TREE_FILE, TREE_DIR)

    # ---- simulate every (k, replicate) -------------------------------------
    blocks, merged, offset = [], {}, 0
    for k in K_VALUES:
        # At k = n_full there is only one possible subset, so replicates would be
        # identical draws; one is enough.
        combos = list(itertools.combinations(range(n_full), k))
        reps = 1 if k == n_full else min(N_REPLICATES, len(combos))
        rng = random.Random(SEED + k)
        draws = rng.sample(combos, reps)

        for rep, idx in enumerate(draws):
            idx = sorted(idx)
            sub = "/".join(events[i] for i in idx)
            ids = [int(x) for e in (events[i] for i in idx) for x in e.split(",")]
            nbr = len(ids)
            bl = sum(by_id[i].dist for i in ids)

            tag = f"k{k}_r{rep}"
            simdir = os.path.join(OUT, tag)
            shutil.rmtree(simdir, ignore_errors=True)
            print(f"  sim {tag}: events {idx}, {nbr} branches, bl={bl:.4f}", flush=True)

            r = docker(["pcoc_sim.py",
                        "-td", f"/proj/results/pcoc_sim/trees/{GENE}",
                        "-o", f"/proj/results/pcoc_sim/sweep/{GENE}/{tag}",
                        "-m", sub, "-c", str(k),
                        "-n_sites", str(N_SITES),
                        "-nb_sampled_couple", str(N_COUPLES),
                        "--no_clean_seqs", "--no_cleanup", "-cpu", str(CPU)])
            if r.returncode != 0:
                print(f"    FAILED: {r.stderr[-300:]}")
                continue

            fas = sorted(f for f in glob.glob(os.path.join(
                simdir, "RUN_*", "Tree_1", "sequences", "Scenario_1", "*.fa"))
                if not is_null_alignment(f))
            if not fas:
                print("    no convergent alignment produced")
                continue

            n_here = 0
            for f in fas:
                s = read_fasta(f)
                width = len(next(iter(s.values())))
                for name, seq in s.items():
                    merged.setdefault(name, []).append(seq)
                n_here += width
            blocks.append(dict(k=k, replicate=rep, events=";".join(map(str, idx)),
                               n_branches=nbr, branch_length=round(bl, 5),
                               start=offset, end=offset + n_here))
            offset += n_here

    if not blocks:
        sys.exit("no simulated data produced")

    cat = os.path.join(OUT, "all_blocks.fa")
    with open(cat, "w") as fh:
        for name, parts in merged.items():
            fh.write(f">{name}\n{''.join(parts)}\n")
    print(f"\nconcatenated {len(blocks)} blocks -> {offset} sites; "
          f"detecting once with the full {n_full}-event declared scenario", flush=True)

    # ---- one detection, with the scenario the real analysis declared --------
    det = os.path.join(OUT, "detect_declared_full")
    shutil.rmtree(det, ignore_errors=True)
    r = docker(["pcoc_det.py",
                "-t", f"/proj/results/branchlengths/pcoc/{GENE}.treefile",
                "-aa", f"/proj/results/pcoc_sim/sweep/{GENE}/all_blocks.fa",
                "-m", full,
                "-o", f"/proj/results/pcoc_sim/sweep/{GENE}/detect_declared_full",
                "-f", "0.0", "-cpu", str(CPU), "--gamma"])
    res = glob.glob(os.path.join(det, "RUN_*", "*.results.tsv"))
    if not res:
        sys.exit(f"detection failed: {r.stderr[-600:]}")
    post = posteriors(res[0])
    print(f"got {len(post)} posteriors for {offset} sites")

    # ---- split back apart --------------------------------------------------
    rows = []
    for b in blocks:
        p = post[b["start"]:b["end"]]
        if not p:
            continue
        row = dict(gene=GENE, k=b["k"], replicate=b["replicate"], events=b["events"],
                   n_branches=b["n_branches"], branch_length=b["branch_length"],
                   n_sites=len(p))
        for thr in (0.8, 0.9, 0.99):
            row[f"power_{thr}"] = round(sum(1 for x in p if x >= thr) / len(p), 4)
        row["median_posterior"] = round(sorted(p)[len(p) // 2], 4)
        row["max_posterior"] = round(max(p), 4)
        rows.append(row)

    with open(os.path.join(OUT, "power_by_k.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(OUT, "blocks.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(blocks[0].keys()))
        w.writeheader()
        w.writerows(blocks)

    print(f"\n{'k':>3s} {'rep':>4s} {'branches':>9s} {'bl':>8s} "
          f"{'power@0.9':>10s} {'medianPost':>11s}")
    for r_ in rows:
        print(f"{r_['k']:3d} {r_['replicate']:4d} {r_['n_branches']:9d} "
              f"{r_['branch_length']:8.4f} {r_['power_0.9']:10.3f} "
              f"{r_['median_posterior']:11.4f}")

    by_k = {}
    for r_ in rows:
        by_k.setdefault(r_["k"], []).append(r_["power_0.9"])
    print("\nmean power at 0.9 by number of converging lineages:")
    recovered = None
    for k in sorted(by_k):
        m = sum(by_k[k]) / len(by_k[k])
        print(f"  k={k:2d}  {m:.3f}")
        if recovered is None and m >= 0.5:
            recovered = k
    print(f"\nDetection recovers at k = {recovered}" if recovered
          else "\nDetection never recovers across the tested range")

    _figure(rows)


def _figure(rows):
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.12,
                        subplot_titles=("Power vs number of converging lineages",
                                        "Power vs convergent branch length"))
    fig.add_trace(go.Scatter(x=[r["k"] for r in rows], y=[r["power_0.9"] for r in rows],
                             mode="markers", name="replicate",
                             marker=dict(size=9, opacity=0.75)), row=1, col=1)
    by_k = {}
    for r in rows:
        by_k.setdefault(r["k"], []).append(r["power_0.9"])
    ks = sorted(by_k)
    fig.add_trace(go.Scatter(x=ks, y=[sum(by_k[k]) / len(by_k[k]) for k in ks],
                             mode="lines+markers", name="mean",
                             line=dict(width=3)), row=1, col=1)
    fig.add_trace(go.Scatter(x=[r["branch_length"] for r in rows],
                             y=[r["power_0.9"] for r in rows],
                             mode="markers", name="replicate",
                             marker=dict(size=9, opacity=0.75), showlegend=False),
                  row=1, col=2)
    fig.update_xaxes(title_text="converging lineages (of 10 declared)", row=1, col=1)
    fig.update_xaxes(title_text="convergent branch length (subs/site)", row=1, col=2)
    fig.update_yaxes(title_text="power at posterior 0.9", range=[-0.03, 1.03], row=1, col=1)
    fig.update_yaxes(title_text="power at posterior 0.9", range=[-0.03, 1.03], row=1, col=2)
    fig.update_layout(template="simple_white", width=1200, height=520,
                      title=(f"{GENE}: detection when only some declared lineages "
                             f"actually converged"))
    out = os.path.join(PROJ, "results", "figures",
                       f"partial_convergence_{GENE}.pdf")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.write_image(out)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
