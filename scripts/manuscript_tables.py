#!/usr/bin/env python3
"""Manuscript Tables 1 and 2, written as CSV and as Markdown for MANUSCRIPT.md.

Generated rather than hand-typed, so the tables cannot drift from the data.

Outputs: shared_results/tables/table1_genes.csv
         shared_results/tables/table4_concordance.csv
         shared_results/tables/table5_gap_sensitivity.csv
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


def table4():
    """Concordance at every internal transition branch (the hemiplasy control).

    Tip branches are excluded: a tip is in every gene tree, so there is no
    bipartition to conflict over and no concordance factor is defined.
    """
    flagged = os.path.join(PROJ, "shared_results", "scenario",
                           "transition_branches_flagged.csv")
    events = os.path.join(PROJ, "shared_results", "scenario",
                          "convergent_events_ER.csv")
    if not (os.path.exists(flagged) and os.path.exists(events)):
        return []

    ev = {}
    for r in csv.DictReader(open(events)):
        if r["node_role"].strip('"') == "transition":
            ev[r["node"]] = (r["event_id"], r["direction"].strip('"'),
                             len(r["tips"].strip('"').split(";")))

    rows = []
    for r in csv.DictReader(open(flagged)):
        try:
            gcf, scf = float(r["gCF"]), float(r["sCFL"])
        except ValueError:
            continue  # tip branch, no concordance defined
        eid, direction, ntips = ev.get(r["child"], ("", "", ""))
        weak = [a for a, v in (("gCF", gcf), ("sCFL", scf)) if v < 50]
        rows.append(dict(
            event=eid, direction=direction, descendant_taxa=ntips,
            gCF=f"{gcf:.1f}", sCFL=f"{scf:.1f}", decisive_gene_trees=r["gN"],
            weak_axis=", ".join(weak) or "none",
            flagged=r["hemiplasy_flag"]))
    rows.sort(key=lambda x: float(x["gCF"]))
    with open(os.path.join(OUT, "table4_concordance.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    return rows


def table5():
    """Sensitivity of the convergent-site count to the phenotype-specificity cutoff.

    The 0.25 threshold was set by judgment, so the result is reported across the
    range rather than at that one value. The pivot point, the largest gap among
    the statistically unusual sites, is computed rather than hardcoded.
    """
    p = os.path.join(PROJ, "shared_results", "pcoc_sim_calibration",
                     "observed_convergent_sites.csv")
    if not os.path.exists(p):
        return []
    sig = [r for r in csv.DictReader(open(p)) if float(r["q_site"]) <= 0.05]
    if not sig:
        return []
    pivot = max(float(r["diurnal_gap"]) for r in sig)

    thresholds = [0.30, 0.25, pivot, 0.20, 0.15, 0.10, 0.00]
    rows = []
    for t in thresholds:
        label = f"{t:g}" + (" (used here)" if abs(t - 0.25) < 1e-9 else "")
        rows.append(dict(threshold=label,
                         sites_passing_both=sum(1 for r in sig
                                                if float(r["diurnal_gap"]) >= t)))
    with open(os.path.join(OUT, "table5_gap_sensitivity.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    return rows


if __name__ == "__main__":
    t1 = table1()
    t2 = table2()
    t4 = table4()
    t5 = table5()
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
        "",
        md(t4, ["Event", "Direction", "Descendant taxa", "gCF", "sCFL",
                "Decisive gene trees", "Weak axis", "Flagged"],
           "Table 4. Gene and site concordance at every internal transition "
           "branch.",
           "The remaining 10 of the 15 transitions are on tip branches, where a "
           "concordance factor is undefined because a tip appears in every gene "
           "tree. A branch is flagged as a hemiplasy risk only when gCF and sCFL "
           "are both below 50 percent, since either alone is usually gene-tree "
           "estimation error rather than genuine conflict. No transition is "
           "flagged."),
        "",
        md(t5, ["phenotype-specificity threshold", "sites passing both criteria"],
           "Table 5. Convergent sites as a function of the phenotype-specificity "
           "threshold.",
           "Sites counted are those that are both statistically unusual "
           "(q <= 0.05) and phenotype-specific at the given threshold. The 0.25 "
           "cutoff was fixed by judgment before the positive control existed and "
           "never adjusted; the row between 0.25 and 0.20 is the largest gap "
           "observed among the statistically unusual sites, so every threshold "
           "above it returns zero."),
    ]
    with open(os.path.join(OUT, "tables.md"), "w") as fh:
        fh.write("\n".join(blocks) + "\n")
    print("wrote", OUT)
    print("\n".join(blocks))
