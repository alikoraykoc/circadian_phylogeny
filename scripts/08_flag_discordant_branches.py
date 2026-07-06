#!/usr/bin/env python3
"""
08_flag_discordant_branches.py

Cross-reference the diel transition branches (from 07) with the gene-concordance
factors (from 06). A diel transition sitting on a low-gCF branch is where
hemiplasy is most likely to manufacture false convergence, so those branches get
flagged for extra scrutiny or exclusion from the scenario.

This is the honest control most convergence papers skip. The one fiddly part is
matching branch identities between corHMM's node numbering and IQ-TREE's
concordance output; both reference the SAME Upham topology, so match by the set
of descendant tips of each branch (a bipartition), which is numbering-agnostic.
"""
import os
import sys
import pandas as pd

PROJ = os.environ.get("PROJ", os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
RES = os.path.join(PROJ, "results")

trans_csv = os.path.join(RES, "scenario", "transition_branches.csv")
cf_stat   = os.path.join(RES, "concordance", "concord.cf.stat")

GCF_FLAG_THRESHOLD = 50.0   # EDIT: flag transition branches with gCF below this (%)

def main():
    if not (os.path.exists(trans_csv) and os.path.exists(cf_stat)):
        sys.exit("Run steps 06 and 07 first (need transition_branches.csv and concord.cf.stat).")

    trans = pd.read_csv(trans_csv)

    # concord.cf.stat is whitespace-delimited with a comment header block.
    cf = pd.read_csv(cf_stat, sep=r"\s+", comment="#")
    # Typical columns include: ID, gCF, gDF1, ..., sCF, ..., Label / branch descriptors.

    # MATCHING STRATEGY (implement to your data):
    #   1. From the labelled Upham tree, compute the descendant-tip set of each
    #      transition branch's child node.
    #   2. From concord.cf.tree, compute the descendant-tip set of each annotated
    #      branch. Match branches by identical tip sets (bipartitions).
    #   3. Join gCF/sCF onto the transition table.
    #
    # Left as an explicit TODO because it depends on how you exported node labels
    # in step 07. The ete3 snippet below is the intended shape.

    print("TODO: implement bipartition matching (see docstring).")
    print(f"Transition branches: {len(trans)} | gCF threshold: {GCF_FLAG_THRESHOLD}%")

    # ---- intended output once matched ----
    # merged = trans.join(gcf_by_branch)
    # merged["hemiplasy_flag"] = merged["gCF"] < GCF_FLAG_THRESHOLD
    # merged.to_csv(os.path.join(RES, "scenario", "transition_branches_flagged.csv"), index=False)
    # n = int(merged["hemiplasy_flag"].sum())
    # print(f"Flagged {n} transition branch(es) on low-concordance branches.")

if __name__ == "__main__":
    main()
