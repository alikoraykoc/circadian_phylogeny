#!/usr/bin/env python3
"""
Build per-gene PCOC scenario strings (-m flag) by matching species-tree
transition branches (from corHMM step 07) to PCOC's own node numbering
in each gene tree. Matching is by descendant-tip bipartitions.
"""
import os
import csv
import re
from ete3 import Tree

PROJ = os.environ.get("PROJ", os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
RES = os.path.join(PROJ, "results")
GENES = "CLOCK NPAS2 ARNTL PER1 PER2 PER3 CRY1 CRY2 NR1D1 NR1D2 RORA RORB RORC CSNK1D CSNK1E FBXL3 BHLHE40 BHLHE41".split()

# 1. Build species tree bipartitions for transition branches
sptree = Tree(os.path.join(PROJ, "data", "species_tree.nwk"), format=1)
tip_order = [l.name for l in sptree.get_leaves()]
ntip = len(tip_order)

# Assign ape-style IDs to species tree
ape_id = {}
for i, name in enumerate(tip_order):
    for leaf in sptree.get_leaves():
        if leaf.name == name:
            ape_id[id(leaf)] = i + 1
            break
counter = ntip + 1
ape_id[id(sptree)] = counter
counter += 1
for node in sptree.traverse("preorder"):
    if node.is_leaf():
        continue
    if id(node) not in ape_id:
        ape_id[id(node)] = counter
        counter += 1

# Map ape node ID -> descendant tips
ape_to_tips = {}
for node in sptree.traverse():
    aid = ape_id[id(node)]
    ape_to_tips[aid] = frozenset(l.name for l in node.get_leaves())

# Read transition branches
trans_csv = os.path.join(RES, "scenario", "transition_branches.csv")
transitions = []
with open(trans_csv) as f:
    for row in csv.DictReader(f):
        child_id = int(row["child"])
        tips = ape_to_tips.get(child_id, frozenset())
        transitions.append(tips)

print(f"Loaded {len(transitions)} transition branches from species tree")

# 2. For each gene, parse PCOC numbered tree and map transitions
def parse_pcoc_numbered_tree(nwk_str):
    """Parse a PCOC-numbered newick. PCOC appends _N to each node name.
    Returns dict: node_number -> frozenset of original tip names."""
    t = Tree(nwk_str, format=1)
    node_map = {}
    # Rename tips: strip the _N suffix to get original name
    tip_rename = {}
    for leaf in t.get_leaves():
        # PCOC format: OriginalName_N
        parts = leaf.name.rsplit("_", 1)
        if len(parts) == 2 and parts[1].isdigit():
            pcoc_id = int(parts[1])
            orig_name = parts[0]
            # But species names can have underscores (Genus_species),
            # so the last _N is the PCOC number
            tip_rename[leaf.name] = (orig_name, pcoc_id)
        else:
            # Fallback
            tip_rename[leaf.name] = (leaf.name, -1)

    # For internal nodes, PCOC numbers them too
    for node in t.traverse():
        if node.is_leaf():
            orig, pcoc_id = tip_rename[node.name]
            tips = frozenset([orig])
            node_map[pcoc_id] = tips
        else:
            # Internal node name is the PCOC number
            try:
                pcoc_id = int(node.name)
            except (ValueError, TypeError):
                continue
            tips = frozenset()
            for leaf in node.get_leaves():
                orig, _ = tip_rename[leaf.name]
                tips = tips | {orig}
            node_map[pcoc_id] = tips
    return node_map

out_dir = os.path.join(RES, "pcoc", "scenarios")
os.makedirs(out_dir, exist_ok=True)

for gene in GENES:
    nwk_path = os.path.join(RES, "pcoc", f"{gene}_num.nwk")
    tree_path = os.path.join(RES, "branchlengths", "pcoc", f"{gene}.treefile")

    if not os.path.exists(nwk_path):
        # Generate numbered tree
        print(f"  Generating numbered tree for {gene}...")
        os.system(
            f'docker run --rm -v {PROJ}:/proj carinerey/pcoc '
            f'pcoc_num_tree.py -t /proj/results/branchlengths/pcoc/{gene}.treefile '
            f'-o /proj/results/pcoc/{gene}_num.pdf -n '
            f'-u /proj/results/pcoc/{gene}_num.nwk 2>/dev/null'
        )

    if not os.path.exists(nwk_path):
        print(f"  WARNING: could not generate numbered tree for {gene}, skipping")
        continue

    with open(nwk_path) as f:
        nwk_str = f.read().strip()

    pcoc_map = parse_pcoc_numbered_tree(nwk_str)
    gene_tips = set()
    for tips in pcoc_map.values():
        gene_tips |= tips

    # Invert: tips -> pcoc_id
    tips_to_pcoc = {v: k for k, v in pcoc_map.items()}

    # Match each transition to a PCOC node
    matched = []
    skipped = 0
    for trans_tips in transitions:
        # If this gene doesn't have all species in this transition clade,
        # intersect with available tips
        avail = trans_tips & gene_tips
        if len(avail) == 0:
            skipped += 1
            continue
        if len(avail) == 1:
            # Single tip transition: use the tip's PCOC number
            tip_name = next(iter(avail))
            for pcoc_id, ptips in pcoc_map.items():
                if ptips == frozenset([tip_name]):
                    matched.append(str(pcoc_id))
                    break
        else:
            # Internal node: find the MRCA of available tips in the gene tree
            pcoc_id = tips_to_pcoc.get(frozenset(avail))
            if pcoc_id is not None:
                matched.append(str(pcoc_id))
            else:
                # Try the exact transition tips that are in the gene
                # The MRCA may include extra tips not in the transition
                skipped += 1

    scenario = "/".join(matched)

    # Write scenario file
    with open(os.path.join(out_dir, f"{gene}.scenario"), "w") as f:
        f.write(scenario + "\n")

    print(f"  {gene}: {len(matched)} transitions matched, {skipped} skipped -> {scenario}")

print(f"\nScenarios written to {out_dir}")
