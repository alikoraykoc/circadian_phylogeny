#!/usr/bin/env python3
"""Manuscript Tables 1 and 2, written as CSV and as Markdown for MANUSCRIPT.md.

Generated rather than hand-typed, so the tables cannot drift from the data.

Outputs: shared_results/tables/table1_genes.csv
         shared_results/tables/table2_scenarios.csv
         shared_results/tables/tables.md
"""
import csv
import glob
import os
import re
from collections import Counter, defaultdict

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(PROJ, "shared_results", "tables")
os.makedirs(OUT, exist_ok=True)

GENES = ("CLOCK NPAS2 ARNTL PER1 PER2 PER3 CRY1 CRY2 NR1D1 NR1D2 "
         "RORA RORB RORC CSNK1D CSNK1E FBXL3 BHLHE40 BHLHE41").split()

MODULE = {g: m for m, gs in {
    "Positive arm": "CLOCK NPAS2 ARNTL",
    "Negative arm": "PER1 PER2 PER3 CRY1 CRY2",
    "Auxiliary loop": "NR1D1 NR1D2 RORA RORB RORC",
    "Post-translational": "CSNK1D CSNK1E FBXL3",
    "Output repressors": "BHLHE40 BHLHE41",
}.items() for g in gs.split()}


def fasta_dims(path):
    n, first, cur = 0, None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            n += 1
            if n == 2:
                first = "".join(cur)
        elif line and n == 1:
            cur.append(line)
    if first is None:
        first = "".join(cur)
    return n, len(first)


AA = set("ACDEFGHIKLMNPQRSTVWY")


def site_classes(path):
    """Constant, variable and parsimony-informative column counts.

    A column is parsimony-informative when at least two distinct amino acids
    each occur in at least two sequences. Gaps and ambiguity codes are ignored,
    matching the convention IQ-TREE reports under the same name.
    """
    rows, cur = [], []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if cur:
                rows.append("".join(cur)); cur = []
        elif line:
            cur.append(line)
    if cur:
        rows.append("".join(cur))
    L = len(rows[0])
    const = pinf = 0
    for j in range(L):
        cnt = Counter(r[j] for r in rows if r[j] in AA)
        if len(cnt) <= 1:
            const += 1
        if sum(1 for v in cnt.values() if v >= 2) >= 2:
            pinf += 1
    return const, L - const, pinf


def tdg09_testable(gene, root):
    p = os.path.join(root, f"{gene}.tdg09.out")
    if not os.path.exists(p):
        return None
    txt = open(p).read()
    if "FullResults:" not in txt:
        return None
    n = 0
    for line in txt.split("FullResults:")[1].splitlines():
        if line.startswith("- ["):
            parts = [x.strip() for x in line.strip()[3:-1].split(",")]
            if len(parts) >= 2 and parts[-1] != "NA":
                n += 1
        elif line.strip() and not line.startswith("#"):
            break
    return n


def table1():
    qc = {}
    p = os.path.join(PROJ, "shared_results", "selection", "foreground_qc.csv")
    if os.path.exists(p):
        qc = {r["gene"]: r for r in csv.DictReader(open(p))}
    rows = []
    for g in GENES:
        un = os.path.join(PROJ, "data", "alignments", f"{g}_aligned.fa")
        tr = os.path.join(PROJ, "results", "trim", f"{g}.trim.fa")
        if not (os.path.exists(un) and os.path.exists(tr)):
            continue
        n_un, l_un = fasta_dims(un)
        n_tr, l_tr = fasta_dims(tr)
        _const, n_var, n_pinf = site_classes(tr)
        # TDG09 returns NA at every column that is not parsimony-informative, so
        # its testable set must equal the parsimony-informative set exactly.
        # Verified site by site across all 18 genes; asserted here so a change in
        # trimming, taxon sampling or TDG09 version cannot silently break it.
        n_td = tdg09_testable(g, os.path.join(PROJ, "results", "tdg09"))
        if n_td is not None and n_td != n_pinf:
            raise SystemExit(
                f"{g}: TDG09 testable sites ({n_td}) != parsimony-informative "
                f"columns ({n_pinf}); the two are expected to be identical")
        rows.append(dict(
            gene=g, module=MODULE.get(g, ""), taxa=n_tr,
            untrimmed_columns=l_un, trimmed_columns=l_tr,
            variable_columns=n_var, parsimony_informative_sites=n_pinf,
            codons_analysed=qc.get(g, {}).get("codons_used", "")))
    with open(os.path.join(OUT, "table1_genes.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    return rows


def scenario_counts(model):
    p = os.path.join(PROJ, "results", "scenario", f"convergent_events_{model}.csv")
    if not os.path.exists(p):
        return {}
    ev, br, tips = defaultdict(set), Counter(), defaultdict(set)
    for r in csv.DictReader(open(p)):
        d = r["direction"].strip('"')
        ev[d].add(r["event_id"]); br[d] += 1
        tips[d] |= set(r["tips"].strip('"').split(";"))
    return {d: (len(ev[d]), br[d], len(tips[d])) for d in ev}


def calibration(setname):
    p = os.path.join(PROJ, "results", "pcoc_sim", setname, "threshold_choice.txt")
    if not os.path.exists(p):
        p = os.path.join(PROJ, "shared_results", f"pcoc_sim_calibration_{setname}",
                         "threshold_choice.txt")
    if not os.path.exists(p):
        return None, None
    pw, fpr = [], []
    for line in open(p):
        m = re.search(r"power=([\d.]+).*FPR=([\d.]+)", line)
        if m:
            pw.append(float(m.group(1))); fpr.append(float(m.group(2)))
    return (min(pw) if pw else None), (max(fpr) if fpr else None)


def sites_above(setname):
    pat = os.path.join(PROJ, "results", "pcoc", setname, "*", "RUN_*", "*.trim.results.tsv")
    files = glob.glob(pat)
    if not files:
        return None
    n = 0
    for f in files:
        rd = list(csv.DictReader(open(f), delimiter="\t"))
        if not rd:
            continue
        cols = [c for c in rd[0] if c.startswith("PCOC")]
        n += sum(1 for r in rd for c in cols if float(r[c] or 0) >= 0.99)
    return n


def table2():
    total_branches = 118
    rows = []
    for model, setname, direction in (("ER", "ER_gain", "gain"),
                                      ("ER", "ER_reversal", "reversal"),
                                      ("ARD", "ARD_gain", "gain"),
                                      ("ARD", "ARD_reversal", "reversal")):
        sc = scenario_counts(model).get(direction)
        if not sc:
            continue
        n_ev, n_br, n_tips = sc
        pw, fpr = calibration(setname)
        sa = sites_above(setname)
        rows.append(dict(
            reconstruction=model,
            direction="Gains of diurnality" if direction == "gain" else "Reversals to nocturnality",
            events=n_ev, convergent_branches=n_br,
            pct_of_branches=round(100 * n_br / total_branches, 1),
            convergent_leaves=n_tips,
            min_gene_power="" if pw is None else f"{pw:.3f}",
            worst_gene_fpr="" if fpr is None else f"{fpr:.4f}",
            sites_above_threshold="not run" if sa is None else sa))
    with open(os.path.join(OUT, "table2_scenarios.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    return rows


def md(rows, header, title, note):
    keys = list(rows[0].keys())
    out = [f"**{title}**", "", note, "",
           "| " + " | ".join(header) + " |",
           "|" + "|".join("---" for _ in keys) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(r[k]) for k in keys) + " |")
    return "\n".join(out)


if __name__ == "__main__":
    t1 = table1()
    t2 = table2()
    blocks = [
        md(t1, ["Gene", "Module", "Taxa", "Untrimmed columns", "Trimmed columns",
                "Variable columns", "Parsimony-informative sites",
                "Codons analysed"],
           "Table 1. Gene set, taxon occupancy and alignment dimensions.",
           "Trimmed columns are those retained by trimAl `-automated1`. A column "
           "is variable when it holds more than one amino acid, and "
           "parsimony-informative when at least two amino acids each occur in at "
           "least two sequences. The parsimony-informative columns are exactly "
           "the sites TDG09 could test; it returns NA at every other column. "
           "Codons analysed are those remaining after removal of all-gap "
           "columns."),
        "",
        md(t2, ["Reconstruction", "Direction", "Events", "Convergent branches",
                "Percent of branches", "Convergent leaves", "Minimum gene power",
                "Worst gene FPR", "Sites above threshold"],
           "Table 2. Convergent scenarios, calibrated detection power and outcome.",
           "Branch and leaf counts are on the 60-taxon species tree (118 branches). "
           "Power and false positive rate are the worst value across the 18 genes at "
           "the calibrated posterior threshold of 0.99. All sites above threshold "
           "were subsequently rejected as alignment artefacts."),
    ]
    with open(os.path.join(OUT, "tables.md"), "w") as fh:
        fh.write("\n".join(blocks) + "\n")
    print("wrote", OUT)
    print("\n".join(blocks))
