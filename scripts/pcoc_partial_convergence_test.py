#!/usr/bin/env python3
"""Measure PCOC power when only SOME of the declared transitions truly converged.

Why this is not covered by run_pcoc_sim.sh
------------------------------------------
run_pcoc_sim.sh simulates and detects with the same 10-event scenario, so it
measures power against COMPLETE convergence: every one of the 52 convergent
branches carries the derived profile. That is the easy case, and PCOC scores a
perfect 1.0 on it.

The realistic failure mode is different. Our detection sweep DECLARES all 10
gains of diurnality as one convergent class. If a given site actually converged
in only 3 of those 10 lineages, PCOC is fitting a convergent profile to 7 clades
that never had it, and the likelihood of the convergent model is dragged down by
the mismatch. A null could then mean "convergence happened, but only in a
subset" rather than "no convergence".

pcoc_sim's -p_conv flag looks like it would test this, but it is declared in the
argument parser and never read anywhere in the codebase, so it does nothing.
This script does it in two stages instead:

  1. simulate with a sub-scenario of k events (only k clades converge)
  2. detect with the FULL 10-event scenario, exactly as the real sweep did

k = 10 is the positive control and should reproduce power 1.0. Each k is also
detected with the CORRECT k-event scenario, which is the upper bound: the gap
between the two curves is the cost of declaring events that did not converge.

Profile couples are concatenated into a single alignment before detection.
pcoc_det's runtime is dominated by fixed per-tree setup rather than by the
number of columns, so one 600-site run costs far less than three 200-site runs,
and PCOC treats sites independently so concatenation changes nothing.

Usage: python scripts/pcoc_partial_convergence_test.py [GENE]   (default CLOCK)
Output: results/pcoc_sim/partial/<gene>/partial_power.csv
"""
import csv
import glob
import os
import random
import shutil
import subprocess
import sys

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IMG = "carinerey/pcoc"
GENE = sys.argv[1] if len(sys.argv) > 1 else "CLOCK"

K_VALUES = [2, 3, 5, 7, 10]
N_SITES = 200
N_COUPLES = 3
CPU = 4
SEED = 20260813

OUT = os.path.join(PROJ, "results", "pcoc_sim", "partial", GENE)
SCEN_FILE = os.path.join(PROJ, "results", "pcoc", "scenarios", "ER_gain", f"{GENE}.scenario")
TREE_DIR = os.path.join(PROJ, "results", "pcoc_sim", "trees", GENE)


def docker(args):
    cmd = ["docker", "run", "--rm", "-v", f"{PROJ}:/proj", IMG] + args
    return subprocess.run(cmd, capture_output=True, text=True)


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


def main():
    full = open(SCEN_FILE).read().strip()
    events = full.split("/")
    n_full = len(events)
    print(f"{GENE}: full scenario has {n_full} events, "
          f"{sum(len(e.split(',')) for e in events)} convergent branches")

    os.makedirs(OUT, exist_ok=True)
    os.makedirs(TREE_DIR, exist_ok=True)
    shutil.copy(os.path.join(PROJ, "results", "branchlengths", "pcoc", f"{GENE}.treefile"),
                TREE_DIR)

    rows = []
    for k in K_VALUES:
        rng = random.Random(SEED + k)
        # Sample which events truly converge. Order is preserved so the
        # transition node stays first within each event group, which pcoc
        # requires.
        idx = sorted(rng.sample(range(n_full), k)) if k < n_full else list(range(n_full))
        sub = "/".join(events[i] for i in idx)
        sub_branches = sum(len(events[i].split(",")) for i in idx)
        simdir = os.path.join(OUT, f"k{k}")
        shutil.rmtree(simdir, ignore_errors=True)

        print(f"\n=== k={k}: simulating {k} converging events "
              f"({sub_branches} branches), events {idx}", flush=True)
        r = docker(["pcoc_sim.py",
                    "-td", f"/proj/results/pcoc_sim/trees/{GENE}",
                    "-o", f"/proj/results/pcoc_sim/partial/{GENE}/k{k}",
                    "-m", sub, "-c", str(k),
                    "-n_sites", str(N_SITES),
                    "-nb_sampled_couple", str(N_COUPLES),
                    "--no_clean_seqs", "--no_cleanup", "-cpu", str(CPU)])
        if r.returncode != 0:
            print(f"  sim FAILED: {r.stderr[-500:]}")
            continue

        # Keep only the alignments whose ancestral and convergent profiles
        # differ; the A<x>_C<x> file is the null and is not needed here, the
        # main calibration already pins FPR at 0.
        fas = [f for f in glob.glob(os.path.join(simdir, "RUN_*", "Tree_1",
                                                 "sequences", "Scenario_1", "*.fa"))
               if not _same_profile(f)]
        if not fas:
            print("  no convergent alignments produced")
            continue

        merged, per_file = {}, []
        for f in sorted(fas):
            s = read_fasta(f)
            n = len(next(iter(s.values())))
            per_file.append((os.path.basename(f), n))
            for name, seq in s.items():
                merged.setdefault(name, []).append(seq)
        cat = os.path.join(OUT, f"k{k}_convergent.fa")
        with open(cat, "w") as fh:
            for name, parts in merged.items():
                fh.write(f">{name}\n{''.join(parts)}\n")
        total_sites = sum(n for _, n in per_file)
        print(f"  concatenated {len(fas)} couples -> {total_sites} convergent sites")

        # Detect twice: with the full declared scenario (what the real sweep
        # did) and with the correct k-event scenario (the achievable ceiling).
        for tag, scen in (("declared_full", full), ("correct_k", sub)):
            det = os.path.join(OUT, f"k{k}_{tag}")
            shutil.rmtree(det, ignore_errors=True)
            print(f"  detecting [{tag}] ...", flush=True)
            r = docker(["pcoc_det.py",
                        "-t", f"/proj/results/branchlengths/pcoc/{GENE}.treefile",
                        "-aa", f"/proj/results/pcoc_sim/partial/{GENE}/k{k}_convergent.fa",
                        "-m", scen,
                        "-o", f"/proj/results/pcoc_sim/partial/{GENE}/k{k}_{tag}",
                        "-f", "0.0", "-cpu", str(CPU), "--gamma"])
            res = glob.glob(os.path.join(det, "RUN_*", "*.results.tsv"))
            if not res:
                print(f"    det FAILED: {r.stderr[-400:]}")
                continue
            post = _posteriors(res[0])
            row = dict(gene=GENE, k=k, n_events_declared=n_full if tag == "declared_full" else k,
                       converging_events=k, converging_branches=sub_branches,
                       detection=tag, n_sites=len(post))
            for thr in (0.8, 0.9, 0.95, 0.99):
                row[f"power_{thr}"] = round(sum(1 for p in post if p >= thr) / len(post), 4) if post else 0.0
            row["median_posterior"] = round(sorted(post)[len(post) // 2], 4) if post else 0.0
            row["max_posterior"] = round(max(post), 4) if post else 0.0
            rows.append(row)
            print(f"    power@0.9={row['power_0.9']:.3f}  median post={row['median_posterior']:.4f}")

    if rows:
        out = os.path.join(OUT, "partial_power.csv")
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print(f"\nwrote {out}")


def _same_profile(path):
    """True for the A<x>_C<x> null alignment, where both profiles are identical."""
    base = os.path.basename(path).replace(".fa", "")
    parts = base.split("_")
    a = [p for p in parts if p.startswith("A")][-1][1:]
    c = [p for p in parts if p.startswith("C")][-1][1:]
    return a == c


def _posteriors(tsv):
    """PCOC_V1 column; blank cells are posteriors that underflowed to 0."""
    out = []
    with open(tsv) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        col = next((c for c in rd.fieldnames if c.startswith("PCOC")), None)
        for r in rd:
            v = (r.get(col) or "").strip()
            out.append(float(v) if v else 0.0)
    return out


if __name__ == "__main__":
    main()
