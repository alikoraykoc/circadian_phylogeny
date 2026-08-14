"""Assertions for the boundaries between pipeline components.

Nine defects have been found in this pipeline. Not one raised an error, and every
one was the same kind of mistake: an unchecked assumption that two
representations of the same object agreed.

  ape node numbering        against  ete3 node numbering
  species-tree clades       against  gene-tree clades under a different rooting
  node labels in memory     against  node labels written to disk
  a CLI flag's argument order against a tree's root state
  alignment columns before pruning against after
  a tool's progress log     against  its results table

Two of them had the pipeline assuming a diurnal ancestral mammal, which is the
exact inverse of the hypothesis under test, and they survived multiple sessions
and reviews.

Every function here RAISES rather than warns. That is deliberate: eight of the
nine defects produced states that looked entirely survivable, and a warning is
precisely what let the rooting bug hide for weeks (prep_tdg09_inputs.py had a
majority-rule fallback that fired silently whenever the exact clade lookup
failed).

Usage:
    from lib_checks import check_same_taxa
    check_same_taxa(aln_names, tree_tips, context=f"{gene}: PCOC input")
"""


class PipelineCheckError(AssertionError):
    """Raised when two representations of the same object disagree."""


def _fail(context, msg, detail=""):
    text = f"{context}: {msg}"
    if detail:
        text += f"\n  {detail}"
    raise PipelineCheckError(text)


def check_same_taxa(names_a, names_b, context, label_a="alignment", label_b="tree"):
    """Two taxon sets must be identical, not merely overlapping.

    CLAUDE.md names label mismatch as the top cause of downstream failure. The
    real instance here was Gorilla_gorilla_gorilla against Gorilla_gorilla, which
    silently dropped one species from every gene in codon_export.
    """
    a, b = set(names_a), set(names_b)
    if a == b:
        return
    only_a = sorted(a - b)[:5]
    only_b = sorted(b - a)[:5]
    _fail(context, f"{label_a} and {label_b} taxa differ "
                   f"({len(a)} vs {len(b)})",
          f"only in {label_a}: {only_a}\n  only in {label_b}: {only_b}")


def outgroup_clade(tree):
    """The smaller of the root's two clades. None if the tree is unrooted."""
    kids = tree.children
    if len(kids) != 2:
        return None            # trifurcating root: IQ-TREE's unrooted output
    return min((frozenset(c.get_leaf_names()) for c in kids), key=len)


def check_same_rooting(tree, reference_tree, context):
    """Two trees must agree on where the root is.

    THE rooting bug. Scenario building matches species-tree clades onto gene
    trees by descendant tip set, and a tip set is not preserved across a
    different rooting. IQ-TREE emits unrooted trees (a trifurcating basal node)
    while the species tree is rooted, so three nodes of the largest gain event
    failed to match in every gene, one of them that event's transition branch.
    Checking at the point of use would have caught it immediately.
    """
    og_ref = outgroup_clade(reference_tree)
    if og_ref is None:
        _fail(context, "reference tree is unrooted (trifurcating root)")
    og = outgroup_clade(tree)
    if og is None:
        _fail(context, "tree is UNROOTED (trifurcating root); "
                       "re-root it before matching clades",
              f"root children sizes: {sorted(len(c) for c in tree.children)}")
    tips = set(tree.get_leaf_names())
    expected = og_ref & tips
    if not expected:
        return                 # no outgroup taxa retained; nothing to compare
    if og != expected and (tips - og) != expected:
        _fail(context, "tree and reference disagree on the root",
              f"tree outgroup has {len(og)} tips, reference implies {len(expected)}")


def check_scenario_against_tree(scenario, node_ids_in_tree, transition_ids,
                                expected_branches, context):
    """A PCOC scenario must reference real nodes, transition-first, none lost.

    Catches two distinct historical defects at once:
      - the dropped-transition bug, where a group began with a DESCENDANT
        because the true transition node failed to match across rootings, and
        PCOC requires the transition node first;
      - the silently-dropped-events bug, where exact tip-set matching lost
        branches and nobody noticed because the scenario still looked well formed.
    """
    groups = [g for g in scenario.strip().split("/") if g]
    if not groups:
        _fail(context, "scenario is empty")

    seen = 0
    for i, grp in enumerate(groups):
        try:
            ids = [int(x) for x in grp.split(",")]
        except ValueError:
            _fail(context, f"group {i} is not a comma-separated list of ints", grp[:60])
        missing = [n for n in ids if n not in node_ids_in_tree]
        if missing:
            _fail(context, f"group {i} references nodes absent from the tree", str(missing[:5]))
        if ids[0] not in transition_ids:
            _fail(context, f"group {i} does not begin with a transition node",
                  f"first id {ids[0]} is not in the transition set; "
                  f"PCOC requires the transition branch first")
        seen += len(ids)

    if expected_branches is not None and seen != expected_branches:
        _fail(context, f"scenario has {seen} branches, expected {expected_branches}",
              "branches were dropped or duplicated during matching")


def check_coordinate_map(mapping, domain_size, target_max, context,
                         require_total=True):
    """A coordinate map must be total over its domain and land in range.

    13_consensus.py chains trimmed column -> untrimmed column -> codon site ->
    reference residue, and had no checks at all. A map that silently loses keys
    produces a results table that is simply missing sites, which looks like a
    negative result.
    """
    if not mapping:
        _fail(context, "coordinate map is empty")
    bad = [k for k, v in mapping.items() if not (1 <= v <= target_max)]
    if bad:
        _fail(context, f"{len(bad)} mapped values fall outside 1..{target_max}",
              f"e.g. {[(k, mapping[k]) for k in bad[:5]]}")
    if require_total and domain_size is not None and len(mapping) != domain_size:
        _fail(context, f"map covers {len(mapping)} of {domain_size} positions",
              "an incomplete map drops sites silently")


def check_parsed_rows(n_parsed, n_expected, context, tolerance=0):
    """A parser must return as many rows as the input claims to contain.

    THE TDG09 trap. Its per-site '# Result{...}' comment lines report lrt=0.0 at
    every site, including ones whose two model likelihoods differ by 9.7 log
    units, because they are a progress log written before the test is evaluated.
    The real numbers live in the FullResults block. Parsing the comments yields a
    clean all-zero null that looks completely plausible, and nothing about the
    output announces the mistake.
    """
    if n_expected is None:
        return
    if abs(n_parsed - n_expected) > tolerance:
        _fail(context, f"parsed {n_parsed} rows, input declares {n_expected}",
              "check whether the parser is reading a progress log rather than "
              "the results table")


def check_no_fallbacks(n_fallbacks, context):
    """A fallback path must not fire in normal operation.

    prep_tdg09_inputs.py falls back to majority-rule labelling when a gene-tree
    clade has no species-tree counterpart. Under the rooting bug that fired
    routinely and silently, mislabelling 7 internal nodes per gene and
    overcounting diurnal ancestors by a third. A fallback that fires quietly is
    indistinguishable from one that never fires.
    """
    if n_fallbacks:
        _fail(context, f"{n_fallbacks} positions used the fallback path",
              "the primary lookup failed; investigate rather than accept the "
              "fallback's answer")
