#!/usr/bin/env python3
"""
Prepare TDG09 inputs: prefix-labelled PHYLIP alignments and trees with
pre-labelled internal nodes. TDG09 determines groups by a two-letter prefix
on each sequence/tip name, and also reads internal node labels if present.
We use 'Di' for diurnal and 'No' for nocturnal.

Internal node states come straight from corHMM's marginal reconstruction, via
the tip-set table `results/scenario/node_tipsets_ER.csv` written by step 07.

The previous version re-derived them instead, by starting at an assumed
nocturnal root and flipping state at every branch listed in
transition_branches.csv. That was fragile in two ways. It rebuilt a hand-rolled
ape node numbering (assuming ete3's leaf order matches ape's tip order, which is
not guaranteed) and it propagated any single error in the transition list
through the whole subtree below it. Reading the states directly removes both,
and matching by DESCENDANT TIP SET rather than node id makes the lookup
invariant to numbering and to per-gene pruning. That is the same defect class as
the rooting bug (PROGRESS_REPORT.md 0A.9), where tip-set matching silently
failed because two trees disagreed about what a clade was.
"""
import os
import csv
from ete3 import Tree

PROJ = os.environ.get("PROJ", os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
RES = os.path.join(PROJ, "results")
GENES = "CLOCK NPAS2 ARNTL PER1 PER2 PER3 CRY1 CRY2 NR1D1 NR1D2 RORA RORB RORC CSNK1D CSNK1E FBXL3 BHLHE40 BHLHE41".split()

# Load diel activity
diel = {}
with open(os.path.join(PROJ, "data", "diel_activity.csv")) as f:
    for row in csv.DictReader(f):
        prefix = "Di" if row["activity"] == "diurnal" else "No"
        diel[row["species"]] = prefix

# corHMM's marginal ancestral states, keyed by descendant tip set.
# Written by 07_scenario.R; see this module's docstring for why the states are
# read rather than re-derived by flipping at transitions.
tipset_state = {}
with open(os.path.join(RES, "scenario", "node_tipsets_ER.csv")) as f:
    for row in csv.DictReader(f):
        tips = frozenset(t for t in row["tips"].split(";") if t)
        tipset_state[tips] = "Di" if row["state"] == "diurnal" else "No"

sptree = Tree(os.path.join(PROJ, "data", "species_tree.nwk"), format=1)

out_dir = os.path.join(RES, "tdg09", "inputs")
os.makedirs(out_dir, exist_ok=True)

for gene in GENES:
    # Read trimmed FASTA
    aln_path = os.path.join(RES, "trim", f"{gene}.trim.fa")
    seqs = {}
    cur = None
    with open(aln_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                cur = line[1:]
                seqs[cur] = ""
            elif cur:
                seqs[cur] += line

    ntaxa = len(seqs)
    nsites = len(next(iter(seqs.values())))

    # Write PHYLIP with prefixed names
    phy_path = os.path.join(out_dir, f"{gene}.phy")
    with open(phy_path, "w") as f:
        f.write(f" {ntaxa} {nsites}\n")
        for sp, seq in seqs.items():
            prefix = diel.get(sp, "No")
            prefixed = f"{prefix}_{sp}"
            f.write(f"{prefixed}  {seq}\n")

    # Read PCOC tree and label ALL nodes (tips + internal)
    tree_path = os.path.join(RES, "branchlengths", "pcoc", f"{gene}.treefile")
    t = Tree(tree_path, format=1)

    # Prefix tip labels
    for leaf in t.get_leaves():
        prefix = diel.get(leaf.name, "No")
        leaf.name = f"{prefix}_{leaf.name}"

    gene_tips_orig = {l.name.split("_", 1)[1] for l in t.get_leaves()}

    # Label internal nodes from the corHMM states, matched by descendant tip set.
    # The gene trees are now re-rooted to the Upham rooting (0A.9), so a clade
    # here is the same clade there and the exact lookup succeeds; before that fix
    # it silently fell through to the majority-rule fallback below.
    fallback = 0
    for node in t.traverse():
        if node.is_leaf():
            continue
        desc = frozenset(l.name.split("_", 1)[1] for l in node.get_leaves())
        state = tipset_state.get(desc)
        if state is None:
            # The node is not a species-tree clade in its own right, which
            # happens when a gene is missing taxa. Fall back to the smallest
            # species-tree clade containing exactly these tips.
            cand = [(s, st) for s, st in tipset_state.items()
                    if (s & set(gene_tips_orig)) == desc]
            state = min(cand, key=lambda x: len(x[0]))[1] if cand else None
        if state is None:
            di = sum(1 for l in node.get_leaves() if l.name.startswith("Di_"))
            state = "Di" if di > len(node) / 2 else "No"
            fallback += 1
        node.name = state

    # Write tree with internal node labels (format=8: all names, all distances)
    tree_out = os.path.join(out_dir, f"{gene}.tree")
    # format_root_node=True so the ROOT's label is actually serialised; ete3 omits
    # it otherwise. TDG09 turns out not to read it (see 10_tdg09.sh, where the
    # root group is set by the ORDER of -groups), but writing a tree whose root
    # is unlabelled while every other node is labelled is a trap for the next
    # reader, and for any other tool pointed at these files.
    t.write(outfile=tree_out, format=8, format_root_node=True)

    # A non-zero fallback count means some internal node could not be matched to
    # any corHMM node and was labelled by majority rule instead. That is a
    # warning, not a routine outcome: it is how the rooting bug hid for so long.
    flag = f"  WARNING: {fallback} nodes labelled by majority-rule fallback" if fallback else ""
    print(f"{gene}: {ntaxa} taxa, {nsites} sites{flag}")

# Write groups file
with open(os.path.join(out_dir, "groups.txt"), "w") as f:
    f.write("No Di\n")

print(f"\nTDG09 inputs written to {out_dir}")
