#!/usr/bin/env python3
"""Build TDG09 inputs from the spiked alignments.

A thin wrapper that reuses prep_tdg09_inputs.py's logic against a different
alignment directory and output directory, so the control exercises the same node
labelling code as the real run. The labelling is the part that carried two of the
nine defects (states re-derived by flipping from an assumed root, and a root
group taken from a CLI flag's argument order), which is exactly why it must not
be reimplemented here.

Env: TRIM_DIR alignments to read, RES_OVERRIDE where results/tdg09/inputs goes.
"""
import csv
import os

from ete3 import Tree

PROJ = os.environ.get("PROJ", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RES = os.path.join(PROJ, "results")
TRIM = os.environ["TRIM_DIR"]
OUTROOT = os.environ["RES_OVERRIDE"]

GENES = ("CLOCK NPAS2 ARNTL PER1 PER2 PER3 CRY1 CRY2 NR1D1 NR1D2 "
         "RORA RORB RORC CSNK1D CSNK1E FBXL3 BHLHE40 BHLHE41").split()

diel = {}
with open(os.path.join(PROJ, "data", "diel_activity.csv")) as f:
    for row in csv.DictReader(f):
        diel[row["species"]] = "Di" if row["activity"] == "diurnal" else "No"

tipset_state = {}
with open(os.path.join(RES, "scenario", "node_tipsets_ER.csv")) as f:
    for row in csv.DictReader(f):
        tips = frozenset(t for t in row["tips"].split(";") if t)
        tipset_state[tips] = "Di" if row["state"] == "diurnal" else "No"

out_dir = os.path.join(OUTROOT, "tdg09", "inputs")
os.makedirs(out_dir, exist_ok=True)

for gene in GENES:
    aln_path = os.path.join(TRIM, f"{gene}.trim.fa")
    if not os.path.exists(aln_path):
        continue
    seqs, cur = {}, None
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
    with open(os.path.join(out_dir, f"{gene}.phy"), "w") as f:
        f.write(f" {ntaxa} {nsites}\n")
        for sp, seq in seqs.items():
            f.write(f"{diel.get(sp, 'No')}_{sp}  {seq}\n")

    t = Tree(os.path.join(RES, "branchlengths", "pcoc", f"{gene}.treefile"), format=1)
    for leaf in t.get_leaves():
        leaf.name = f"{diel.get(leaf.name, 'No')}_{leaf.name}"
    gene_tips = {l.name.split("_", 1)[1] for l in t.get_leaves()}

    fallback = 0
    for node in t.traverse():
        if node.is_leaf():
            continue
        desc = frozenset(l.name.split("_", 1)[1] for l in node.get_leaves())
        state = tipset_state.get(desc)
        if state is None:
            cand = [(s, st) for s, st in tipset_state.items()
                    if (s & gene_tips) == desc]
            state = min(cand, key=lambda x: len(x[0]))[1] if cand else None
        if state is None:
            di = sum(1 for l in node.get_leaves() if l.name.startswith("Di_"))
            state = "Di" if di > len(node) / 2 else "No"
            fallback += 1
        node.name = state

    t.write(outfile=os.path.join(out_dir, f"{gene}.tree"), format=8,
            format_root_node=True)
    flag = f"  WARNING: {fallback} majority-rule fallbacks" if fallback else ""
    print(f"{gene}: {ntaxa} taxa, {nsites} sites{flag}")

with open(os.path.join(out_dir, "groups.txt"), "w") as f:
    f.write("No Di\n")
print(f"\nspiked TDG09 inputs in {out_dir}")
