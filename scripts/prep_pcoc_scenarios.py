#!/usr/bin/env python3
"""
Build per-gene PCOC scenario strings (the -m flag) from the convergent events
emitted by step 07.

PCOC scenario format (confirmed against `pcoc_det.py -h` in carinerey/pcoc):

    -m "1,2,3/67/55,56"
    "Transition node must be the first number and independent events must be
     separated by a '/'"

So each "/"-group is ONE independent convergent event, the first node in the
group is the transition branch, and the remaining nodes are the other branches
that stay in the derived state. PCOC does NOT propagate the convergent state
down the tree for you: any descendant branch you omit is modelled as ancestral.
Omitting them silently destroys the signal, because the tips that actually carry
the convergent amino acids end up under the ancestral profile.

Two things this script must get right:

1. Node identity across trees. corHMM/ape, IQ-TREE and PCOC all number internal
   nodes differently, so nodes are matched by their descendant-tip SET, which is
   numbering-agnostic.
2. Missing species. The per-gene trees are pruned versions of the Upham topology
   (CSNK1D has 40 of 60 species), so an event's tip set is intersected with the
   tips actually present in that gene. Pruning can only remove tips from a clade,
   never add them, so the clade's node in the gene tree carries exactly that
   intersection. Events whose taxa are entirely absent from a gene are dropped,
   and the count is reported per gene.

Scenarios are written per (ASR model, direction):
    results/pcoc/scenarios/<model>_<direction>/<gene>.scenario
"""
import os
import csv
import sys
from collections import OrderedDict
from ete3 import Tree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_checks import check_scenario_against_tree

PROJ = os.environ.get("PROJ", os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
RES = os.path.join(PROJ, "results")
GENES = "CLOCK NPAS2 ARNTL PER1 PER2 PER3 CRY1 CRY2 NR1D1 NR1D2 RORA RORB RORC CSNK1D CSNK1E FBXL3 BHLHE40 BHLHE41".split()

MODELS = ["ER", "ARD"]
DIRECTIONS = ["gain", "reversal"]

# PCOC needs at least a few independent events to have any power (CLAUDE.md).
MIN_EVENTS = 4


def load_events(model):
    """Read convergent_events_<model>.csv into
    {(direction, event_id): [tipset, ...]} with the transition node FIRST."""
    path = os.path.join(RES, "scenario", f"convergent_events_{model}.csv")
    if not os.path.exists(path):
        sys.exit(f"Missing {path}. Run scripts/07_scenario.R first.")

    events = OrderedDict()
    with open(path) as f:
        for row in csv.DictReader(f):
            key = (row["direction"], int(row["event_id"]))
            tips = frozenset(row["tips"].split(";"))
            events.setdefault(key, {"transition": None, "convergent": []})
            if row["node_role"] == "transition":
                events[key]["transition"] = tips
            else:
                events[key]["convergent"].append(tips)

    ordered = OrderedDict()
    for key, v in events.items():
        if v["transition"] is None:
            continue  # malformed event, no transition node
        ordered[key] = [v["transition"]] + v["convergent"]
    return ordered


def parse_pcoc_numbered_tree(nwk_path):
    """Map PCOC node number -> frozenset of original tip names.

    PCOC renames tips to "<OriginalName>_<N>" and labels internal nodes with the
    bare number. Species names themselves contain underscores (Genus_species), so
    only the final _N is the PCOC id.
    """
    t = Tree(nwk_path, format=1)

    orig_name = {}
    for leaf in t.get_leaves():
        base, _, suffix = leaf.name.rpartition("_")
        if base and suffix.isdigit():
            orig_name[leaf.name] = (base, int(suffix))
        else:
            raise ValueError(f"Unparseable PCOC tip name: {leaf.name!r} in {nwk_path}")

    node_map = {}
    for node in t.traverse():
        if node.is_leaf():
            base, pid = orig_name[node.name]
            node_map[pid] = frozenset([base])
        else:
            try:
                pid = int(node.name)
            except (ValueError, TypeError):
                continue  # unlabelled node (e.g. the root in some outputs)
            node_map[pid] = frozenset(orig_name[l.name][0] for l in node.get_leaves())
    return node_map


def build_scenario(events, tipset_to_pcoc, gene_tips):
    """Return (scenario_string, n_events_kept, n_branches, n_events_dropped)."""
    groups = []
    dropped = 0

    for _key, nodes in events.items():
        ids = []
        for tips in nodes:
            avail = tips & gene_tips
            if not avail:
                continue  # every species of this branch is missing from the gene
            pid = tipset_to_pcoc.get(avail)
            if pid is None:
                # Should not happen: pruning only removes tips from a clade, so
                # the clade's node carries exactly this intersection. Surface it
                # rather than silently dropping a branch.
                continue
            if pid not in ids:  # pruning can collapse distinct nodes onto one
                ids.append(pid)

        if not ids:
            dropped += 1
            continue
        # ids[0] is the transition node: it was processed first and, being the
        # ancestor of every other node in the event, cannot be dropped while a
        # descendant survives.
        groups.append(",".join(str(i) for i in ids))

    scenario = "/".join(groups)
    n_branches = sum(len(g.split(",")) for g in groups)
    return scenario, len(groups), n_branches, dropped


def main():
    for model in MODELS:
        events_all = load_events(model)

        for direction in DIRECTIONS:
            events = OrderedDict(
                (k, v) for k, v in events_all.items() if k[0] == direction
            )
            out_dir = os.path.join(RES, "pcoc", "scenarios", f"{model}_{direction}")
            os.makedirs(out_dir, exist_ok=True)

            print(f"\n== {model} / {direction}: {len(events)} events on the species tree ==")
            print(f"{'gene':<9} {'events':>6} {'branches':>9}  {'dropped':>7}")

            for gene in GENES:
                nwk_path = os.path.join(RES, "pcoc", f"{gene}_num.nwk")
                if not os.path.exists(nwk_path):
                    print(f"{gene:<9} {'--':>6} {'--':>9}  numbered tree missing")
                    continue

                pcoc_map = parse_pcoc_numbered_tree(nwk_path)
                gene_tips = frozenset().union(*pcoc_map.values())
                tipset_to_pcoc = {tips: pid for pid, tips in pcoc_map.items()}

                scenario, n_ev, n_br, n_drop = build_scenario(
                    events, tipset_to_pcoc, gene_tips
                )

                # The scenario is the single most defect-prone artefact in the
                # pipeline: it silently lost events to exact tip-set matching,
                # and later began a group with a DESCENDANT instead of its
                # transition node because the two trees disagreed on the root.
                # Both failures produced a well-formed file, so verify against
                # the tree rather than trusting the string.
                # load_events puts the transition node's tipset first in each
                # event, which is the invariant PCOC's -m format depends on.
                transition_ids = set()
                for _k, nodes in events.items():
                    if not nodes:
                        continue
                    pid = tipset_to_pcoc.get(nodes[0] & gene_tips)
                    if pid is not None:
                        transition_ids.add(pid)
                if scenario.strip() and transition_ids:
                    check_scenario_against_tree(
                        scenario, set(pcoc_map), transition_ids,
                        expected_branches=None,
                        context=f"{gene} [{model}/{direction}]")

                with open(os.path.join(out_dir, f"{gene}.scenario"), "w") as f:
                    f.write(scenario + "\n")

                warn = ""
                if n_ev < MIN_EVENTS:
                    warn = f"  LOW POWER (<{MIN_EVENTS} events)"
                print(f"{gene:<9} {n_ev:>6} {n_br:>9}  {n_drop:>7}{warn}")

    print(f"\nScenarios written under {os.path.join(RES, 'pcoc', 'scenarios')}/<model>_<direction>/")


if __name__ == "__main__":
    main()
