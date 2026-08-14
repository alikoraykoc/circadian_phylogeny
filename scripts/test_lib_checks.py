#!/usr/bin/env python3
"""Prove each check FIRES on the historical bug it was written for.

A check that has never failed is not known to work. Every case below
reconstructs an actual defect from this project, not a hypothetical one, and
asserts that the corresponding check raises. Each case also asserts that the
check PASSES on the corrected input, so a function that raises unconditionally
cannot masquerade as working.

Run: python scripts/test_lib_checks.py
"""
import os
import sys

from ete3 import Tree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_checks import (PipelineCheckError, check_coordinate_map,  # noqa: E402
                        check_no_fallbacks, check_parsed_rows,
                        check_same_rooting, check_same_taxa,
                        check_scenario_against_tree)

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PASS, FAIL = [], []


def expect_raise(name, fn, *a, **kw):
    try:
        fn(*a, **kw)
    except PipelineCheckError as e:
        PASS.append(f"{name}: fired -> {str(e).splitlines()[0][:88]}")
        return
    FAIL.append(f"{name}: DID NOT FIRE (the check is asleep)")


def expect_pass(name, fn, *a, **kw):
    try:
        fn(*a, **kw)
        PASS.append(f"{name}: passes on corrected input")
    except PipelineCheckError as e:
        FAIL.append(f"{name}: FALSE ALARM on good input -> {e}")


def main():
    # ---- 1. Gorilla_gorilla_gorilla vs Gorilla_gorilla ---------------------
    real = {"Homo_sapiens", "Gorilla_gorilla", "Mus_musculus"}
    codon_export = {"Homo_sapiens", "Gorilla_gorilla_gorilla", "Mus_musculus"}
    expect_raise("taxa mismatch (Gorilla trinomial)",
                 check_same_taxa, codon_export, real, "codon_export vs species tree")
    expect_pass("taxa mismatch (after rename)",
                check_same_taxa, real, real, "renamed")

    # ---- 2. the rooting bug ------------------------------------------------
    sp = Tree(os.path.join(PROJ, "data", "species_tree.nwk"), format=1)
    unrooted = os.path.join(PROJ, "results", "branchlengths", "pcoc",
                            "ARNTL.treefile.unrooted")
    rooted = os.path.join(PROJ, "results", "branchlengths", "pcoc", "ARNTL.treefile")
    if os.path.exists(unrooted):
        expect_raise("unrooted gene tree (THE rooting bug)",
                     check_same_rooting, Tree(unrooted, format=1), sp,
                     "ARNTL: pre-fix tree")
    else:
        FAIL.append("rooting: no .unrooted backup to test against")
    if os.path.exists(rooted):
        expect_pass("re-rooted gene tree",
                    check_same_rooting, Tree(rooted, format=1), sp, "ARNTL: fixed")

    # ---- 3. scenario with a descendant first, not the transition ----------
    ids = set(range(120))
    transitions = {94, 105, 0, 40, 22}
    good = "94,93,90/105/0"
    bad_order = "93,94,90/105/0"          # transition no longer first
    bad_missing = "94,93,90/105/0/999"    # node absent from the tree
    expect_pass("scenario well formed",
                check_scenario_against_tree, good, ids, transitions, 5, "ok")
    expect_raise("scenario begins with a descendant (dropped-transition bug)",
                 check_scenario_against_tree, bad_order, ids, transitions, 5, "bad order")
    expect_raise("scenario references a missing node",
                 check_scenario_against_tree, bad_missing, ids, transitions, 6, "missing")
    expect_raise("scenario lost branches (52 of 55)",
                 check_scenario_against_tree, good, ids, transitions, 55, "short count")

    # ---- 4. incomplete coordinate map --------------------------------------
    total = {i: i for i in range(1, 11)}
    partial = {i: i for i in range(1, 8)}
    out_of_range = dict(total)
    out_of_range[5] = 999
    expect_pass("coordinate map total", check_coordinate_map, total, 10, 10, "ok")
    expect_raise("coordinate map drops positions silently",
                 check_coordinate_map, partial, 10, 10, "incomplete")
    expect_raise("coordinate map lands out of range",
                 check_coordinate_map, out_of_range, 10, 10, "out of range")

    # ---- 5. the TDG09 progress-log trap ------------------------------------
    # CLOCK declares SiteCount 841. The FullResults block has 841 rows; the
    # '# Result{...}' comments are a progress log with lrt=0.0 everywhere.
    expect_pass("parser row count matches", check_parsed_rows, 841, 841, "CLOCK FullResults")
    expect_raise("parser read the wrong block",
                 check_parsed_rows, 177, 841, "CLOCK: parsed only variable sites")

    # ---- 6. a fallback that fires silently ---------------------------------
    expect_pass("no fallbacks used", check_no_fallbacks, 0, "prep_tdg09_inputs")
    expect_raise("majority-rule fallback fired (hid the rooting bug)",
                 check_no_fallbacks, 7, "prep_tdg09_inputs: 7 nodes")

    # ---- report ------------------------------------------------------------
    for line in PASS:
        print(f"  ok    {line}")
    for line in FAIL:
        print(f"  FAIL  {line}")
    print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
    if FAIL:
        sys.exit(1)
    print("Every check fires on its historical bug and passes on corrected input.")


if __name__ == "__main__":
    main()
