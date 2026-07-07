#!/usr/bin/env python3
"""
Prepare TDG09 inputs: prefix-labelled PHYLIP alignments and trees with
pre-labelled internal nodes. TDG09 determines groups by a two-letter prefix
on each sequence/tip name, and also reads internal node labels if present.
We use 'Di' for diurnal and 'No' for nocturnal, and label internal nodes
using the MAP ancestral states from corHMM (step 07).
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

# Load transition branches to reconstruct internal node states
# We rebuild ancestral states: start from root (nocturnal = ancestral),
# then flip at each transition branch
trans_csv = os.path.join(RES, "scenario", "transition_branches.csv")
transitions = set()
with open(trans_csv) as f:
    for row in csv.DictReader(f):
        transitions.add(int(row["child"]))

# Build ape-style node states on the species tree
sptree = Tree(os.path.join(PROJ, "data", "species_tree.nwk"), format=1)
tip_order = [l.name for l in sptree.get_leaves()]
ntip = len(tip_order)

# Assign ape-style IDs
ape_id_map = {}  # ete3 node id -> ape id
for i, name in enumerate(tip_order):
    for leaf in sptree.get_leaves():
        if leaf.name == name:
            ape_id_map[id(leaf)] = i + 1
            break

counter = ntip + 1
ape_id_map[id(sptree)] = counter
counter += 1
for node in sptree.traverse("preorder"):
    if node.is_leaf():
        continue
    if id(node) not in ape_id_map:
        ape_id_map[id(node)] = counter
        counter += 1

# Propagate states: root is nocturnal (0), transitions flip
# state 1 = nocturnal (corHMM coding: 0+1=1), state 2 = diurnal (1+1=2)
# In our binary: 0 = nocturnal, 1 = diurnal
node_state = {}  # ape_id -> "Di" or "No"

# First set tip states
for leaf in sptree.get_leaves():
    ape_id = ape_id_map[id(leaf)]
    node_state[ape_id] = diel.get(leaf.name, "No")

# For internal nodes, reconstruct from transitions
# Root is nocturnal (ancestral state)
def assign_states(node, parent_state):
    ape_id = ape_id_map[id(node)]
    if ape_id in transitions:
        # Flip state
        current = "Di" if parent_state == "No" else "No"
    else:
        current = parent_state
    if not node.is_leaf():
        node_state[ape_id] = current
        for child in node.children:
            assign_states(child, current)

assign_states(sptree, "No")  # root = nocturnal

# Build reverse map: ete3 node -> state
ete_node_state = {}
for node in sptree.traverse():
    ape_id = ape_id_map[id(node)]
    if ape_id in node_state:
        ete_node_state[id(node)] = node_state[ape_id]

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

    # Map species tree internal node states onto gene tree by descendant tip matching
    sp_bip_state = {}
    for node in sptree.traverse():
        if node.is_leaf():
            continue
        tips = frozenset(l.name for l in node.get_leaves())
        ape_id = ape_id_map[id(node)]
        if ape_id in node_state:
            sp_bip_state[tips] = node_state[ape_id]

    # Prefix tip labels
    for leaf in t.get_leaves():
        prefix = diel.get(leaf.name, "No")
        leaf.name = f"{prefix}_{leaf.name}"

    # Label internal nodes by matching bipartitions
    gene_tips_orig = {l.name.split("_", 1)[1] for l in t.get_leaves()}
    for node in t.traverse():
        if node.is_leaf():
            continue
        desc_tips = frozenset(l.name.split("_", 1)[1] for l in node.get_leaves())
        state = sp_bip_state.get(desc_tips)
        if state:
            node.name = state
        else:
            # Fallback: majority rule from descendant tips
            di_count = sum(1 for l in node.get_leaves() if l.name.startswith("Di_"))
            node.name = "Di" if di_count > len(list(node.get_leaves())) / 2 else "No"

    # Write tree with internal node labels (format=8: all names, all distances)
    tree_out = os.path.join(out_dir, f"{gene}.tree")
    t.write(outfile=tree_out, format=8)

    print(f"{gene}: {ntaxa} taxa, {nsites} sites")

# Write groups file
with open(os.path.join(out_dir, "groups.txt"), "w") as f:
    f.write("Di No\n")

print(f"\nTDG09 inputs written to {out_dir}")
