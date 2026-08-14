#!/usr/bin/env python3
"""Re-root every per-gene tree to the Upham rooting.

CLAUDE.md states the invariant plainly: "Re-root each per-gene tree to the Upham
rooting before scenario building." It was never implemented.
`05_branchlengths_fixed.sh:74` only echoes a reminder to do it, so the trees have
been used exactly as IQ-TREE emitted them, which is unrooted: a trifurcating
basal node with children of sizes 1, 1, 58 against the species tree's 6, 54.

What that cost
--------------
Scenario building matches species-tree clades onto gene trees by descendant tip
set. Under a different rooting a clade's descendant set is not preserved, so
three nodes of the largest gain event failed to match in EVERY gene. Worse, one
of the three is that event's TRANSITION node, so the event group in each
`.scenario` file began with a descendant instead of the transition branch. PCOC
requires the transition node first, and `prep_pcoc_scenarios.py` assumed the
ancestor could never be dropped while its descendants survived. That assumption
holds only if the two trees share a rooting.

Effect per gene: 52 of 55 gain branches declared, and event 1 (24 of the 55, the
largest gain) given the wrong transition branch.

The marsupial outgroup is monophyletic in every gene tree, so re-rooting is
well defined and lossless. Originals are preserved as <gene>.treefile.unrooted.

Re-rooting changes PCOC's postorder node numbering, so scenarios must be rebuilt
and any analysis keyed to scenario node ids must be re-run.

Usage: python scripts/reroot_gene_trees.py [--check]
  --check  report what would change without writing anything
"""
import os
import shutil
import sys

from ete3 import Tree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_checks import PipelineCheckError, check_same_rooting, check_same_taxa

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TREES = os.path.join(PROJ, "results", "branchlengths", "pcoc")
SPTREE = os.path.join(PROJ, "data", "species_tree.nwk")

GENES = ("CLOCK NPAS2 ARNTL PER1 PER2 PER3 CRY1 CRY2 NR1D1 NR1D2 "
         "RORA RORB RORC CSNK1D CSNK1E FBXL3 BHLHE40 BHLHE41").split()

CHECK = "--check" in sys.argv


def outgroup_tips():
    """The smaller of the species tree root's two clades, i.e. the marsupials."""
    sp = Tree(SPTREE, format=1)
    return min((set(c.get_leaf_names()) for c in sp.children), key=len)


def main():
    og = outgroup_tips()
    print(f"Upham outgroup: {len(og)} tips ({', '.join(sorted(og)[:3])}, ...)")
    print(f"mode: {'CHECK ONLY, nothing written' if CHECK else 'REWRITING trees in place'}\n")

    ok = skipped = 0
    for g in GENES:
        path = os.path.join(TREES, f"{g}.treefile")
        if not os.path.exists(path):
            print(f"{g:9s} missing treefile")
            continue
        t = Tree(path, format=1)
        tips = set(t.get_leaf_names())
        here = og & tips
        before = sorted(len(c) for c in t.children)

        if not here:
            print(f"{g:9s} SKIP: no outgroup taxa present")
            skipped += 1
            continue

        # ete3 needs a single node to root on. With more than one outgroup tip
        # that node is their common ancestor, which must be monophyletic or the
        # rooting would silently move unrelated taxa across the root.
        if len(here) == 1:
            target = t & next(iter(here))
        else:
            mono, why, _ = t.check_monophyly(values=list(here), target_attr="name")
            if not mono:
                print(f"{g:9s} SKIP: outgroup not monophyletic ({why})")
                skipped += 1
                continue
            target = t.get_common_ancestor(list(here))

        tips_before = set(t.get_leaf_names())
        t.set_outgroup(target)
        after = sorted(len(c) for c in t.children)

        # Re-rooting must move the root and nothing else. ete3's set_outgroup
        # rearranges the topology around the new root, so verify the taxon set
        # survived and that the result now agrees with the species tree, which
        # is the whole point of doing this.
        check_same_taxa(t.get_leaf_names(), tips_before,
                        f"{g}: after re-rooting", "re-rooted", "original")
        check_same_rooting(t, Tree(SPTREE, format=1), f"{g}: after re-rooting")

        if not CHECK:
            backup = path + ".unrooted"
            if not os.path.exists(backup):
                shutil.copy(path, backup)
            t.write(outfile=path, format=1)

        print(f"{g:9s} root children {before} -> {after}   outgroup {len(here)} tips")
        ok += 1

    print(f"\n{ok} trees re-rooted, {skipped} skipped")
    if not CHECK and ok:
        print("\nNEXT, in order:")
        print("  1. python scripts/prep_pcoc_scenarios.py     (node numbering changed)")
        print("  2. verify every event group starts with its transition node")
        print("  3. re-run PCOC, its calibration, and anything keyed to scenario ids")


if __name__ == "__main__":
    main()
