#!/usr/bin/env python3
"""
08_flag_discordant_branches.py

Cross-reference the diel transition branches (from 07) with the gene-concordance
factors (from 06). A diel transition sitting on a low-gCF branch is where
hemiplasy is most likely to manufacture false convergence, so those branches get
flagged for extra scrutiny or exclusion from the scenario.

Matching strategy: both corHMM (step 07) and IQ-TREE concordance (step 06)
reference the same Upham topology, but use different internal node numbering.
We match branches by their descendant-tip sets (bipartitions), which is
numbering-agnostic.
"""
import os
import sys
import csv
from ete3 import Tree

PROJ = os.environ.get("PROJ", os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
RES = os.path.join(PROJ, "results")

trans_csv = os.path.join(RES, "scenario", "transition_branches.csv")
cf_stat = os.path.join(RES, "concordance", "concord.cf.stat")
cf_branch = os.path.join(RES, "concordance", "concord.cf.tree")
sptree_path = os.path.join(PROJ, "data", "species_tree.nwk")

GCF_FLAG_THRESHOLD = 50.0


def get_descendant_tips(node):
    """Return frozenset of leaf names under a node."""
    return frozenset(l.name for l in node.get_leaves())


def build_ape_bipartitions(tree_path):
    """Build a map from ape-style node ID to descendant tip set.

    ape numbers tips 1..Ntip (alphabetical in the Newick), then internal
    nodes (Ntip+1)..(Ntip+Nnode). The root is Ntip+1. We replicate this
    by reading the tree with ete3 and assigning IDs in ape's traversal order.
    """
    t = Tree(tree_path, format=1)
    tips = sorted(t.get_leaf_names())
    ntip = len(tips)

    # ape tip IDs: 1-based index in the order tips appear in tip.label
    # tip.label order = the order tips are encountered in the Newick (left to right)
    tip_order = [l.name for l in t.get_leaves()]
    tip_id = {name: i + 1 for i, name in enumerate(tip_order)}

    # ape internal node IDs: root = Ntip+1, then postorder
    # We assign by traversal to match ape's read.tree ordering
    internal_id = {}
    counter = ntip + 1
    # ape assigns root first, then in the order edges appear
    internal_id[id(t)] = counter
    counter += 1
    for node in t.traverse("preorder"):
        if node.is_leaf():
            continue
        if id(node) not in internal_id:
            internal_id[id(node)] = counter
            counter += 1

    # Build bipartition map: ape_node_id -> frozenset of tip names
    bip_map = {}
    for node in t.traverse():
        if node.is_leaf():
            ape_id = tip_id[node.name]
        else:
            ape_id = internal_id[id(node)]
        bip_map[ape_id] = get_descendant_tips(node)

    return bip_map, t


def parse_cf_stat(cf_stat_path):
    """Parse concord.cf.stat into a dict of branch_id -> gCF, sCF."""
    cf_data = {}
    with open(cf_stat_path) as f:
        header = None
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if header is None:
                header = parts
                continue
            row = dict(zip(header, parts))
            try:
                bid = int(row["ID"])
                gcf = float(row["gCF"])
                scf = float(row["sCF"])
                gn = int(row["gN"])
                cf_data[bid] = {"gCF": gcf, "sCF": scf, "gN": gn}
            except (ValueError, KeyError):
                continue
    return cf_data


def build_iqtree_bipartitions(cf_branch):
    """Parse the cf.branch tree (which has IQ-TREE node IDs) and build
    a map from IQ-TREE branch ID to descendant tip set."""
    t = Tree(cf_branch, format=1)
    bip_map = {}
    for node in t.traverse():
        if node.is_leaf():
            continue
        try:
            bid = int(node.name)
        except (ValueError, TypeError):
            continue
        bip_map[bid] = get_descendant_tips(node)
    return bip_map


def main():
    if not all(os.path.exists(p) for p in [trans_csv, cf_stat, cf_branch, sptree_path]):
        sys.exit("Run steps 06 and 07 first.")

    # 1. Read transition branches (ape node IDs)
    with open(trans_csv) as f:
        reader = csv.DictReader(f)
        transitions = list(reader)

    # 2. Build bipartition maps for both numbering systems
    ape_bip, ape_tree = build_ape_bipartitions(sptree_path)
    iqt_bip = build_iqtree_bipartitions(cf_branch)
    cf_data = parse_cf_stat(cf_stat)

    # 3. Build IQ-TREE bipartition -> (branch_id, gCF, sCF) lookup
    bip_to_cf = {}
    for bid, tips in iqt_bip.items():
        if bid in cf_data:
            bip_to_cf[tips] = {"iqtree_id": bid, **cf_data[bid]}

    # 4. Match each transition branch to its concordance values
    results = []
    for tr in transitions:
        child_node = int(tr["child"])
        child_tips = ape_bip.get(child_node, frozenset())
        child_label = tr["child_label"] if tr["child_label"] != "" else None

        matched = bip_to_cf.get(child_tips)
        if matched:
            gcf = matched["gCF"]
            scf = matched["sCF"]
            gn = matched["gN"]
            flag = gcf < GCF_FLAG_THRESHOLD
        else:
            # Tip branches don't have concordance factors
            gcf = None
            scf = None
            gn = None
            flag = False

        results.append({
            "edge_id": tr["edge_id"],
            "parent": tr["parent"],
            "child": tr["child"],
            "child_label": child_label or "",
            "child_tips": ",".join(sorted(child_tips)) if len(child_tips) <= 10 else f"{len(child_tips)} taxa",
            "gCF": f"{gcf:.1f}" if gcf is not None else "tip",
            "sCF": f"{scf:.1f}" if scf is not None else "tip",
            "gN": str(gn) if gn is not None else "tip",
            "hemiplasy_flag": flag,
        })

    # 5. Write output
    out_path = os.path.join(RES, "scenario", "transition_branches_flagged.csv")
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    # 6. Report
    n_internal = sum(1 for r in results if r["gCF"] != "tip")
    n_flagged = sum(1 for r in results if r["hemiplasy_flag"])
    n_tip = sum(1 for r in results if r["gCF"] == "tip")

    print(f"Transition branches: {len(results)} total ({n_internal} internal, {n_tip} tip)")
    print(f"gCF threshold: {GCF_FLAG_THRESHOLD}%")
    print(f"Flagged for hemiplasy risk: {n_flagged}")
    print()
    for r in results:
        marker = " ** FLAGGED" if r["hemiplasy_flag"] else ""
        label = r["child_label"] or r["child_tips"]
        print(f"  edge {r['edge_id']:>3}: gCF={r['gCF']:>5}  sCF={r['sCF']:>5}  {label}{marker}")

    print(f"\nOutput: {out_path}")


if __name__ == "__main__":
    main()
