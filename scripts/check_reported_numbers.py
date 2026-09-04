#!/usr/bin/env python3
"""Verify that the numbers quoted in the write-ups still match the data.

Three defects motivated this script, all of the same shape: a number that was
correct when it was written, and stale by the time someone read it.

  1. PROGRESS_REPORT.md quoted 645 / 213 parsimony sites after the hardened run
     had replaced them with 648 / 201.
  2. shared_results/README.md quoted the pre-rooting-fix k-curve as current.
  3. shared_results/figures/ held figures from before three rendering fixes.

Part A recomputes each headline number from the data and fails if the write-ups
no longer contain it. Part B scans for known superseded numbers appearing
without a nearby marker saying they are superseded.

Usage: python scripts/check_reported_numbers.py
Exit status is non-zero if anything is stale.
"""
import csv
import os
import re
import sys
from collections import Counter, defaultdict

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SHARED = os.path.join(PROJ, "shared_results")
AA = set("ACDEFGHIKLMNPQRSTVWY")

GENES = ("CLOCK NPAS2 ARNTL PER1 PER2 PER3 CRY1 CRY2 NR1D1 NR1D2 "
         "RORA RORB RORC CSNK1D CSNK1E FBXL3 BHLHE40 BHLHE41").split()

# Narrative documents whose numbers must track the data. The archive README is
# excluded on purpose: quoting superseded numbers is its entire job.
DOCS = ["MANUSCRIPT.md", "RESULTS.md", "PROGRESS_REPORT.md", "METHODS.md",
        "RESUME.md", "shared_results/README.md",
        "shared_results/rerconverge/README.md"]

problems = []


def doc_text(name):
    p = os.path.join(PROJ, name)
    return open(p).read() if os.path.exists(p) else ""


def claims_check(label, pattern, allowed, *, required=True):
    """Every occurrence of `pattern` must capture a value in `allowed`.

    Anchoring on surrounding wording, rather than asking whether a number
    appears somewhere in the file, is what makes this catch a single stale
    sentence in a document that states the same quantity correctly elsewhere.
    """
    allowed = {str(a) for a in allowed}
    allowed |= {f"{int(a):,}" for a in allowed if str(a).isdigit()}
    seen = 0
    for name in DOCS:
        txt = doc_text(name)
        if not txt:
            continue
        # These documents are hard-wrapped, so a claim regularly straddles a line
        # break. Swap newlines for spaces rather than collapsing whitespace, so
        # character offsets survive and line numbers stay exact.
        flat = txt.replace("\n", " ")
        for m in re.finditer(pattern, flat):
            seen += 1
            got = m.group(1)
            if got.replace(",", "") not in {a.replace(",", "") for a in allowed}:
                line = txt[:m.start()].count("\n") + 1
                problems.append(
                    f"{name}:{line} {label}: found {got}, expected "
                    f"{sorted(allowed)[0] if len(allowed) == 1 else sorted(allowed)}")
    if required and seen == 0:
        problems.append(f"{label}: no document states this any more "
                        f"(pattern {pattern!r} matched nothing)")


def read_fa_rows(path):
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
    return rows


# ---------------------------------------------------------------- Part A
def check_alignment_dimensions():
    total = const = pinf = 0
    for g in GENES:
        p = os.path.join(PROJ, "results", "trim", f"{g}.trim.fa")
        if not os.path.exists(p):
            return
        rows = read_fa_rows(p)
        L = len(rows[0])
        total += L
        for j in range(L):
            cnt = Counter(r[j] for r in rows if r[j] in AA)
            if len(cnt) <= 1:
                const += 1
            if sum(1 for v in cnt.values() if v >= 2) >= 2:
                pinf += 1

    claims_check("trimmed columns",
                 r"trimmed dataset comprises ([\d,]+) columns", [total])
    claims_check("trimmed columns",
                 r"protein alignment: ([\d,]+) trimmed columns", [total])
    claims_check("constant columns", r"([\d,]+) are constant", [const])
    claims_check("variable columns", r"constant, ([\d,]+) variable", [const and var_of(total, const)])
    claims_check("parsimony-informative columns",
                 r"variable and ([\d,]+) parsimony-informative", [pinf])
    # A large "of N sites/columns" denominator must be one the data supports:
    # the whole trimmed alignment, the parsimony-informative subset TDG09 could
    # test, or the codon count. Anything else is a typo or a stale figure.
    claims_check("alignment denominator",
                 r"(?:of|in|across) ([\d,]{5,}) (?:sites|columns)\b",
                 [total, pinf, codon_count()])
    # Two testable denominators are legitimate: the real data and the spike-in,
    # whose planted substitutions make a few more columns informative.
    spike = spikein_pinf()
    claims_check("testable-site denominator",
                 r"of ([\d,]+) testable sites",
                 [pinf] + ([spike] if spike else []))


def var_of(total, const):
    return total - const


def spikein_pinf():
    total = 0
    for g in GENES:
        p = os.path.join(PROJ, "results", "spikein", "trim", f"{g}.trim.fa")
        if not os.path.exists(p):
            return None
        rows = read_fa_rows(p)
        for j in range(len(rows[0])):
            cnt = Counter(r[j] for r in rows if r[j] in AA)
            if sum(1 for v in cnt.values() if v >= 2) >= 2:
                total += 1
    return total


def check_parsimony():
    p = os.path.join(SHARED, "pcoc_sim_calibration",
                     "observed_convergent_substitutions.csv")
    if not os.path.exists(p):
        return
    rows = list(csv.DictReader(open(p)))
    sites = sum(int(r["sites_2plus_events_changed"]) for r in rows)
    same = sum(int(r["convergent_sites"]) for r in rows)
    nul = sum(float(r["null_mean_convergent"]) for r in rows)
    opp = sum(float(r["null_mean_opportunity"]) for r in rows)

    claims_check("parsimony same-residue count",
                 r"([\d,]+) of those changed to the same residue", [same])
    claims_check("parsimony same-residue count",
                 r"(?:of the|of) ([\d,]+) same-residue sites", [same])
    claims_check("parsimony same-residue count",
                 r"same residue in 2\+ independent diurnal lineages \| ([\d,]+)",
                 [same])
    claims_check("parsimony same-residue count",
                 r"Pooled, ([\d,]+) of [\d,]+ opportunity sites", [same])
    claims_check("parsimony opportunity count",
                 r"Pooled, [\d,]+ of ([\d,]+) opportunity sites", [sites])
    claims_check("parsimony opportunity count",
                 r"([\d,]+) sites changed in 2 or more independent diurnal", [sites],
                 required=False)
    claims_check("parsimony observed rate",
                 r"(?:observed rate of|Rate) ([\d.]+)", [f"{same / sites:.3f}"],
                 required=False)

    q = os.path.join(SHARED, "pcoc_sim_calibration", "observed_convergent_sites.csv")
    if os.path.exists(q):
        n = len(list(csv.DictReader(open(q))))
        if n != same:
            problems.append(
                f"per-site table has {n} rows but the summary reports {same} "
                f"same-residue sites")


def check_scenario():
    p = os.path.join(SHARED, "scenario", "convergent_events_ER.csv")
    if not os.path.exists(p):
        return
    ev, br = defaultdict(set), Counter()
    for r in csv.DictReader(open(p)):
        d = r["direction"].strip('"')
        ev[d].add(r["event_id"]); br[d] += 1

    # The k-curve must be built on the same scenario the study reports. This is
    # the check that would have caught the swapped partial_convergence dirs.
    k = os.path.join(SHARED, "partial_convergence", "power_by_k.csv")
    if os.path.exists(k):
        rows = list(csv.DictReader(open(k)))
        top = max(int(r["n_branches"]) for r in rows)
        if top != br["gain"]:
            problems.append(
                f"partial_convergence/ tops out at {top} convergent branches but "
                f"the ER gain scenario has {br['gain']}; the k-curve directory is "
                f"probably the superseded pre-rooting-fix run")


def codon_count():
    p = os.path.join(SHARED, "selection", "contrastfel_sites.csv")
    if not os.path.exists(p):
        return 0
    return sum(1 for _ in csv.DictReader(open(p)))


def check_contrastfel():
    claims_check("Contrast-FEL codon denominator",
                 r"of ([\d,]+) codons", [codon_count()], required=False)


def check_spikein():
    p = os.path.join(SHARED, "spikein", "spikein_scorecard.csv")
    if not os.path.exists(p):
        return
    for r in csv.DictReader(open(p)):
        if r["method"] == "TDG09":
            claims_check("spike-in TDG09 false positives",
                         r"\*\*([\d,]+) false positives\*\* among",
                         [r["false_positives"]], required=False)
            break


def check_permulations():
    """The permulation null must still say what the write-ups claim it says."""
    mins = {}
    for mode in ("ssm", "cc"):
        p = os.path.join(SHARED, "rerconverge", f"rer_permulation_pvalues_{mode}.csv")
        if not os.path.exists(p):
            return
        rows = list(csv.DictReader(open(p)))
        n_sig = sum(1 for r in rows if float(r["permulation_P"]) < 0.05)
        if n_sig:
            problems.append(
                f"{mode} permulations: {n_sig} gene(s) now below p = 0.05, but the "
                f"write-ups state that none reach it")
        mins[mode] = (min(float(r["permulation_P"]) for r in rows),
                      min(float(r["permulation_p_adj"]) for r in rows),
                      min(rows, key=lambda r: float(r["permulation_P"]))["gene"])

    ssm_p, ssm_adj, ssm_gene = mins["ssm"]
    cc_p, cc_adj, _ = mins["cc"]
    claims_check("smallest ssm permulation p",
                 r"(?:permulation p-value being|permulation p is|"
                 r"permulation p-value of|becomes \*\*)(\d+\.\d+)", [f"{ssm_p:.3f}"])
    claims_check("smallest cc permulation p",
                 r"[Cc]omplete-case permulations (?:agree|give)[^(]*\((\d+\.\d+)",
                 [f"{cc_p:.3f}"], required=False)
    claims_check("smallest ssm adjusted permulation p",
                 r"smallest adjusted value(?:s)? (?:is|are) (\d+\.\d+)",
                 [f"{ssm_adj:.3f}"], required=False)

    # The anticonservatism claim is a count, so recompute it rather than trust it.
    par = {r["gene"]: (float(r["parametric_P"]), float(r["permulation_P"]))
           for r in csv.DictReader(open(os.path.join(
               SHARED, "rerconverge", "rer_permulation_pvalues_ssm.csv")))}
    n_anti = sum(1 for a, b in par.values() if a < b)
    claims_check("genes where the parametric p is smaller",
                 r"(?:smaller than the permulation p-value in|"
                 r"anticonservative in) ([\d]+) of (?:the )?18", [n_anti])


# Surnames that are two words, so the "X et al." form ends on the second word
# while the reference entry begins on the first.
COMPOUND_SURNAMES = {"Pond": "Kosakovsky"}


def check_references():
    """Every in-text citation must be listed, and every entry must be cited."""
    txt = doc_text("MANUSCRIPT.md")
    if "# References" not in txt:
        problems.append("MANUSCRIPT.md has no References section")
        return
    body, refs = txt.split("# References", 1)
    # Collapse wrapping so "Boyko and Beaulieu\n2021" reads as one citation.
    flat = re.sub(r"\s+", " ", body)

    cited = set()
    # "Surname et al. 2018" and "Surname and Surname (2015)".
    for m in re.finditer(
            r"\b([A-Z][a-zA-Z-]+) (?:et al\.?|and [A-Z][a-zA-Z-]+) \(?(\d{4})\)?", flat):
        cited.add((COMPOUND_SURNAMES.get(m.group(1), m.group(1)), m.group(2)))
    # Single-author and corporate forms, always parenthesised as "(Name 2023)".
    # The lookbehind keeps this off the second surname of "X and Y 2015".
    for m in re.finditer(
            r"\((?![^()]*\band\b)([A-Z][a-zA-Z-]+)(?: [A-Z][A-Za-z.]+)* (\d{4})\)", flat):
        cited.add((COMPOUND_SURNAMES.get(m.group(1), m.group(1)), m.group(2)))

    listed = set()
    for line in refs.split("\n"):
        m = re.match(r"^([A-Z][a-zA-Z-]+)[^(]*\((\d{4})\)", line)
        if m:
            listed.add((m.group(1), m.group(2)))

    for name, year in sorted(cited - listed):
        problems.append(f"MANUSCRIPT.md cites {name} {year} but the References "
                        f"section does not list it")
    for name, year in sorted(listed - cited):
        problems.append(f"References lists {name} {year} but nothing in the text "
                        f"cites it")


def check_gap_sensitivity():
    """The threshold-sensitivity table must match a recount from the site table.

    The 0.25 phenotype-specificity cutoff was a judgment call, so the manuscript
    reports the result as a function of it. That table is the part a reviewer
    will check, and it must not drift.
    """
    p = os.path.join(SHARED, "pcoc_sim_calibration", "observed_convergent_sites.csv")
    if not os.path.exists(p):
        return
    rows = list(csv.DictReader(open(p)))
    sig = [r for r in rows if float(r["q_site"]) <= 0.05]
    if not sig:
        problems.append("no sites at q <= 0.05; the sensitivity table assumes 12")
        return

    txt = doc_text("MANUSCRIPT.md")
    # Scope to the sensitivity table itself; other tables also hold decimals.
    head = "| phenotype-specificity threshold | sites passing both criteria |"
    if head not in txt:
        problems.append("MANUSCRIPT.md no longer carries the gap-sensitivity table")
        return
    block = txt.split(head, 1)[1].split("\n\n", 1)[0]
    for m in re.finditer(r"\| (0\.\d+)(?: \(used here\))? \| (\d+) \|", block):
        thr, claimed = float(m.group(1)), int(m.group(2))
        actual = sum(1 for r in sig if float(r["diurnal_gap"]) >= thr)
        if actual != claimed:
            problems.append(
                f"gap-sensitivity table: at threshold {thr} the manuscript says "
                f"{claimed} sites, the data give {actual}")

    top = max(float(r["diurnal_gap"]) for r in sig)
    claims_check("largest gap among significant sites",
                 r"largest gap among the 12\s+statistically unusual sites is (\d+\.\d+)",
                 [f"{top:g}"], required=False)

    # The serine enrichment is load-bearing: it is presented as the observable
    # signature of Fitch's unordered counting, so it must still hold.
    n_ser = sum(1 for r in sig if r["residue"] == "S")
    claims_check("serine count among significant sites",
                 r"(\d+)(?: of \d+)? are serine", [n_ser], required=False)
    if n_ser * 2 != len(sig):
        problems.append(
            f"{n_ser} of {len(sig)} significant sites are serine, but the "
            f"write-ups describe them as half")

    # Named borderline sites must still be the ones a relaxed threshold admits.
    for gene, site in (("PER1", "876"), ("BHLHE40", "359")):
        hit = [r for r in sig if r["gene"] == gene and r["site"] == site]
        if not hit:
            problems.append(f"{gene} site {site} is named in the manuscript as a "
                            f"borderline candidate but is no longer significant")


def check_clock675():
    """The CLOCK 675 worked example must still describe the actual column.

    Three write-ups quote it as the clearest case of homoplasy, with specific
    counts. Those counts are recomputed here rather than trusted.
    """
    aln = os.path.join(PROJ, "results", "trim", "CLOCK.trim.fa")
    diel = os.path.join(PROJ, "data", "diel_activity.csv")
    if not (os.path.exists(aln) and os.path.exists(diel)):
        return
    lab = {r["species"]: r["activity"] for r in csv.DictReader(open(diel))}

    names, rows, cur = [], [], []
    for line in open(aln):
        line = line.rstrip()
        if line.startswith(">"):
            if cur:
                rows.append("".join(cur)); cur = []
            names.append(line[1:].split()[0])
        elif line:
            cur.append(line)
    if cur:
        rows.append("".join(cur))

    col = 675 - 1
    states = Counter(r[col] for r in rows)
    n_met, n_tax = states.get("M", 0), len(rows)
    claims_check("CLOCK 675 methionine count",
                 r"fixed for\s+methionine in (\d+) of \d+ species", [n_met])
    claims_check("CLOCK 675 taxon count",
                 r"methionine in \d+ of (\d+)\s+species", [n_tax])

    carriers = [names[i] for i, r in enumerate(rows) if r[col] == "V"]
    by_diel = Counter(lab.get(t, "?") for t in carriers)
    if by_diel.get("diurnal") != 2 or by_diel.get("nocturnal") != 2:
        problems.append(
            f"CLOCK 675 valine carriers are now {dict(by_diel)}, but the "
            f"write-ups say two diurnal and two nocturnal")


def check_hemiplasy():
    """The concordance control must still say what the write-ups claim."""
    p = os.path.join(SHARED, "scenario", "transition_branches_flagged.csv")
    if not os.path.exists(p):
        return
    rows = list(csv.DictReader(open(p)))

    def num(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return None

    internal = [r for r in rows if num(r["gCF"]) is not None]
    terminal = [r for r in rows if num(r["gCF"]) is None]
    both_weak = [r for r in internal
                 if num(r["gCF"]) < 50 and num(r["sCFL"]) < 50]
    one_weak = [r for r in internal
                if (num(r["gCF"]) < 50) != (num(r["sCFL"]) < 50)]
    flagged = [r for r in rows if r["hemiplasy_flag"] == "True"]

    if both_weak or flagged:
        problems.append(
            f"hemiplasy control: {len(both_weak)} transition(s) weak on both axes "
            f"and {len(flagged)} flagged, but the write-ups state none")

    # The write-ups spell these counts as words, so match either form.
    words = {"one": 1, "two": 2, "three": 3, "four": 4, "five": 5, "six": 6,
             "seven": 7, "eight": 8, "nine": 9, "ten": 10, "eleven": 11,
             "twelve": 12, "fifteen": 15}

    def written(pattern, expected, label):
        for name in DOCS:
            txt = doc_text(name)
            if not txt:
                continue
            for m in re.finditer(pattern, txt.replace("\n", " "), re.I):
                tok = m.group(1)
                got = words.get(tok.lower(), None) if not tok.isdigit() else int(tok)
                if got is not None and got != expected:
                    line = txt[:m.start()].count("\n") + 1
                    problems.append(f"{name}:{line} {label}: says {tok}, "
                                    f"data give {expected}")

    written(r"(\w+) of the 15 transitions (?:sit|are) on tip",
            len(terminal), "transitions on tip branches")
    written(r"(\w+) are weak on one axis", len(one_weak),
            "internal transitions weak on one axis")

    # The branch the write-ups single out as the lowest gCF must still be it.
    worst = min(internal, key=lambda r: num(r["gCF"]))
    claims_check("lowest transition gCF",
                 r"gCF is (\d+\.\d+), the lowest of any transition",
                 [f"{num(worst['gCF']):g}"], required=False)
    claims_check("lowest-gCF branch sCFL",
                 r"the lowest of any\s+transition, against\s+sCFL of (\d+\.\d+)",
                 [f"{num(worst['sCFL']):g}"], required=False)


# ---------------------------------------------------------------- Part B
# Numbers that were correct in an earlier generation of the analysis. They may
# appear only next to a marker saying so.
SUPERSEDED = {
    "213": "pre-hardening parsimony same-residue count (current: 201)",
    "206": "pre-rooting-fix parsimony same-residue count (current: 201)",
    "645": "pre-hardening parsimony opportunity count (current: 648)",
    "639": "pre-rooting-fix parsimony opportunity count (current: 648)",
}
MARKERS = ("supersede", "superseded", "pre-rooting", "prereroot", "pre-fix",
           "pre-hardening", "archive", "historical", "do not quote", "v1",
           "replaced", "earlier")


def marker_window(lines, i):
    """The text a marker may legitimately live in, for the number on line i.

    Deliberately narrow: a marker four paragraphs away does not make a number
    on this line labelled. The one widening is markdown tables, where the
    "superseded" label belongs in the header row rather than on each data row,
    so a table row carries its whole table's header with it.
    """
    lo, hi = max(0, i - 2), min(len(lines), i + 3)
    window = lines[lo:hi]
    if lines[i].lstrip().startswith("|"):
        j = i
        while j > 0 and lines[j - 1].lstrip().startswith("|"):
            j -= 1
        window = window + lines[max(0, j - 1):i]
    return " ".join(window).lower()


def check_superseded_numbers():
    for name in DOCS:
        txt = doc_text(name)
        if not txt:
            continue
        lines = txt.split("\n")
        for i, line in enumerate(lines):
            for num, meaning in SUPERSEDED.items():
                if not re.search(rf"(?<![\d.,]){num}(?![\d.,])", line):
                    continue
                if not any(m in marker_window(lines, i) for m in MARKERS):
                    problems.append(
                        f"{name}:{i + 1} quotes {num} ({meaning}) with no nearby "
                        f"marker saying it is superseded")


if __name__ == "__main__":
    check_alignment_dimensions()
    check_parsimony()
    check_scenario()
    check_contrastfel()
    check_spikein()
    check_permulations()
    check_references()
    check_gap_sensitivity()
    check_clock675()
    check_hemiplasy()
    check_superseded_numbers()

    if problems:
        print("STALE NUMBERS FOUND\n")
        for p in problems:
            print("  -", p)
        sys.exit(1)
    print("All reported numbers match the data.")
