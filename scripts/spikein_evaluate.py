#!/usr/bin/env python3
"""Score the pipeline against the planted convergent sites.

This is the ONLY script that should read ANSWER_KEY.csv. Everything upstream
runs blind: the spiked alignments are ordinary alignments as far as PCOC, TDG09,
the parsimony test and the residue screen are concerned.

What it measures, per method and per difficulty level:

  recall     of the planted sites, how many were flagged
  false pos  of the flagged sites, how many were not planted

Recall answers "does signal survive the plumbing", which is the question nine
silent interface defects made necessary. False positives answer a second
question that has been hanging over the results: TDG09 currently calls 23 percent
of variable sites significant while PCOC and Contrast-FEL call none, and nothing
in the real data can adjudicate that because the truth is unknown. Here it is
known.

Expected, if the pipeline is sound:
  level 'all'  (10 of 10 events)  PCOC recovers most; model-free methods too
  level 'most' (7 of 10)          PCOC borderline, per the k-curve in 0A.5b
  level 'few'  (3 of 10)          PCOC near zero BY DESIGN, model-free should hold

A failure at level 'all' would mean signal is being lost in the plumbing, and
the entire null result would be uninterpretable until it is found.

Usage: python scripts/spikein_evaluate.py
"""
import csv
import glob
import json
import os
from collections import defaultdict

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SPIKE = os.path.join(PROJ, "results", "spikein")

PCOC_THRESHOLD = 0.99
FDR_THRESHOLD = 0.05
LEVELS = ("all", "most", "few")


def load_key():
    p = os.path.join(SPIKE, "ANSWER_KEY.csv")
    if not os.path.exists(p):
        raise SystemExit(f"no answer key at {p}; run spikein_generate.py first")
    with open(p) as fh:
        return [dict(r, trimmed_col=int(r["trimmed_col"])) for r in csv.DictReader(fh)]


def load_pcoc_hits():
    """{(gene, col)} above the calibrated posterior threshold."""
    hits, tested = set(), defaultdict(int)
    for f in glob.glob(os.path.join(SPIKE, "pcoc", "*", "RUN_*", "*.trim.results.tsv")):
        gene = os.path.basename(f).split(".")[0]
        with open(f) as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            col = next((c for c in rd.fieldnames if c.startswith("PCOC")), None)
            site_col = next((c for c in rd.fieldnames if c.lower().startswith("site")), None)
            for i, r in enumerate(rd, start=1):
                site = int(r[site_col]) if site_col and r.get(site_col) else i
                tested[gene] += 1
                v = (r.get(col) or "").strip()
                if v and float(v) >= PCOC_THRESHOLD:
                    hits.add((gene, site))
    return hits, tested


def load_tdg09_hits():
    hits, tested = set(), defaultdict(int)
    for f in glob.glob(os.path.join(SPIKE, "tdg09", "*.tdg09.out")):
        gene = os.path.basename(f).split(".")[0]
        inblock = False
        for line in open(f):
            if line.startswith("FullResults:"):
                inblock = True
                continue
            if not inblock:
                continue
            if line.startswith("- ["):
                p = [x.strip() for x in line.strip()[3:-1].split(",")]
                if len(p) < 2 or p[-1] == "NA":
                    continue
                tested[gene] += 1
                try:
                    if float(p[-1]) <= FDR_THRESHOLD:
                        hits.add((gene, int(p[0])))
                except ValueError:
                    pass
            elif line.strip() and not line.startswith("#"):
                break
    return hits, tested


def load_parsimony_hits():
    """Sites the hardened parsimony test calls convergent AND phenotype-specific."""
    p = os.path.join(SPIKE, "observed_convergent_sites.csv")
    if not os.path.exists(p):
        return set(), set()
    flagged, specific = set(), set()
    with open(p) as fh:
        for r in csv.DictReader(fh):
            k = (r["gene"], int(r["site"]))
            flagged.add(k)
            if r.get("phenotype_specific", "").lower() == "true":
                specific.add(k)
    return flagged, specific


def score(name, hits, key_by_site, planted, tested_total):
    """Recall per level, plus how many flagged sites were not planted."""
    rows = []
    for lv in LEVELS:
        want = {k for k, v in key_by_site.items() if v["level"] == lv}
        got = want & hits
        rows.append((lv, len(got), len(want)))
    fp = len(hits - planted)
    return rows, fp


def main():
    key = load_key()
    key_by_site = {(k["gene"], k["trimmed_col"]): k for k in key}
    planted = set(key_by_site)

    pcoc_hits, pcoc_tested = load_pcoc_hits()
    tdg_hits, tdg_tested = load_tdg09_hits()
    pars_flagged, pars_specific = load_parsimony_hits()

    print(f"planted: {len(planted)} sites across {len({g for g, _ in planted})} genes")
    print(f"  by level: " + ", ".join(
        f"{lv}={sum(1 for v in key_by_site.values() if v['level']==lv)}" for lv in LEVELS))
    print()

    methods = [
        ("PCOC", pcoc_hits, sum(pcoc_tested.values())),
        ("TDG09", tdg_hits, sum(tdg_tested.values())),
        ("parsimony (flagged)", pars_flagged, 0),
        ("parsimony (phenotype-specific)", pars_specific, 0),
    ]

    print(f"{'method':32s} {'all':>9s} {'most':>9s} {'few':>9s} {'false pos':>10s}")
    print("-" * 74)
    for name, hits, tested in methods:
        if not hits and tested == 0:
            print(f"{name:32s} {'(not run)':>9s}")
            continue
        rows, fp = score(name, hits, key_by_site, planted, tested)
        cells = " ".join(f"{g:>4d}/{w:<4d}" for _, g, w in rows)
        print(f"{name:32s} {cells} {fp:10d}")

    print()
    for name, hits, tested in methods:
        if tested:
            print(f"  {name}: flagged {len(hits)} of {tested} tested sites "
                  f"({100*len(hits)/tested:.1f} percent)")

    # ---- verdict -----------------------------------------------------------
    rows, _ = score("PCOC", pcoc_hits, key_by_site, planted, 0)
    all_recall = rows[0][1] / rows[0][2] if rows[0][2] else 0
    print()
    if all_recall >= 0.8:
        print(f"PLUMBING OK: PCOC recovered {rows[0][1]}/{rows[0][2]} of the fully")
        print("convergent sites, so signal survives scenario building, tree handling")
        print("and detection end to end. The null on the real data is interpretable.")
    elif all_recall > 0:
        print(f"PARTIAL: PCOC recovered only {rows[0][1]}/{rows[0][2]} fully convergent")
        print("sites. Signal reaches the detector but is being attenuated somewhere.")
    else:
        print("PLUMBING BROKEN: PCOC recovered NONE of the fully convergent sites.")
        print("Signal is lost between the alignment and the detector. The null result")
        print("on the real data means nothing until this is located.")

    out = os.path.join(SPIKE, "spikein_scorecard.csv")
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["method", "level", "recovered", "planted", "false_positives"])
        for name, hits, tested in methods:
            if not hits and tested == 0:
                continue
            rows, fp = score(name, hits, key_by_site, planted, tested)
            for lv, got, want in rows:
                w.writerow([name, lv, got, want, fp])
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
