#!/usr/bin/env python3
"""Build per-gene codon alignments and foreground-labelled trees for HyPhy.

The selection track asks a different question from the convergence track: not
"did the same residue appear in diurnal lineages" but "is the selective regime
different there". Those routinely dissociate, so a convergence null says nothing
either way about selection.

Foreground = every DIURNAL branch, not just the transition branches. The
hypothesis is a sustained regime in diurnal lineages, so the transition branch
plus the descendant branches that stay diurnal all belong to the test set. That
is the same branch set the PCOC convergent class uses.

Why this needs its own script
-----------------------------
Each gene has CDS for a different subset of species, 32 to 58 of 60, so a single
shared labelled tree cannot work. Per gene:

  1. the codon alignment is rewritten with the project's species naming
     (codon_export uses Gorilla_gorilla_gorilla where everything else uses
     Gorilla_gorilla). codon_export/ is left untouched so its provenance stands;
  2. that gene's tree is pruned to its CDS taxa. The topology is the Upham
     scaffold (CLAUDE.md design decision 1), but the branch lengths must come
     from the per-gene tree, NOT from the species tree. The Upham tree is dated
     in millions of years, and handing HyPhy a branch length of 46 as if it were
     substitutions per site means total saturation; in practice it crashed the
     GTR pre-fit with a reversible-model numerical error. The per-gene trees are
     in substitutions per site, which is the scale CLAUDE.md requires;
  3. branches are labelled {Foreground} where the child node's ER state is
     diurnal, matched by DESCENDANT TIP SET rather than by node id.

Tip-set matching is the load-bearing part. ape, ete3, IQ-TREE and PCOC all number
nodes differently, and a node id means nothing once a tree is pruned. It is also
what the rooting bug came down to: matching by tip set is only valid when both
trees share a rooting, which is why the pruned tree is derived from the species
tree itself rather than from a gene tree.

Inputs : codon_export/codon/<gene>.codon.fasta
         results/scenario/node_tipsets_ER.csv   (from 07_scenario.R)
         results/branchlengths/pcoc/<gene>.treefile   (re-rooted, subs/site)
         data/species_tree.nwk                        (state lookup only)
Outputs: results/codon/<gene>.codon.fasta          renamed, HyPhy-ready
         results/selection/trees/<gene>.labelled.nwk
         results/selection/foreground_qc.csv
"""
import csv
import os

from ete3 import Tree

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODON_SRC = os.path.join(PROJ, "codon_export", "codon")
CODON_OUT = os.path.join(PROJ, "results", "codon")
TREE_OUT = os.path.join(PROJ, "results", "selection", "trees")
SPTREE = os.path.join(PROJ, "data", "species_tree.nwk")
GENE_TREES = os.path.join(PROJ, "results", "branchlengths", "pcoc")
TIPSETS = os.path.join(PROJ, "results", "scenario", "node_tipsets_ER.csv")

GENES = ("CLOCK NPAS2 ARNTL PER1 PER2 PER3 CRY1 CRY2 NR1D1 NR1D2 "
         "RORA RORB RORC CSNK1D CSNK1E FBXL3 BHLHE40 BHLHE41").split()

# codon_export uses the subspecies trinomial; the rest of the project uses the
# binomial. Rename at the point of consumption, never in codon_export itself.
RENAME = {"Gorilla_gorilla_gorilla": "Gorilla_gorilla"}

# HyPhy needs enough branches in each class for the test to mean anything.
# Contrast-FEL compares two rate classes, so a handful of foreground branches
# gives an uninterpretable result rather than a negative one.
MIN_FOREGROUND = 5
MIN_BACKGROUND = 5


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        elif line:
            buf.append(line)
    if name:
        seqs[name] = "".join(buf)
    return seqs


def load_states():
    """{frozenset(tips): state} for every node of the species tree."""
    out = {}
    with open(TIPSETS) as fh:
        for r in csv.DictReader(fh):
            tips = frozenset(t for t in r["tips"].split(";") if t)
            out[tips] = r["state"]
    return out


def main():
    os.makedirs(CODON_OUT, exist_ok=True)
    os.makedirs(TREE_OUT, exist_ok=True)
    states = load_states()
    sp = Tree(SPTREE, format=1)
    all_tips = set(sp.get_leaf_names())

    qc = []
    print(f"{'gene':9s} {'taxa':>5s} {'codons':>7s} {'allgap':>8s} {'fg':>4s} {'bg':>4s} "
          f"{'fg_tips':>8s} {'flag':>10s}")
    for g in GENES:
        src = os.path.join(CODON_SRC, f"{g}.codon.fasta")
        if not os.path.exists(src):
            print(f"{g:9s} missing codon alignment")
            continue

        seqs = {RENAME.get(k, k): v for k, v in read_fasta(src).items()}
        unknown = set(seqs) - all_tips
        if unknown:
            raise SystemExit(f"{g}: names absent from the species tree: {sorted(unknown)}")

        width = len(next(iter(seqs.values())))
        if width % 3:
            raise SystemExit(f"{g}: codon alignment width {width} is not a multiple of 3")

        # Prune the GENE tree, which carries substitutions-per-site branch
        # lengths on the Upham topology and has been re-rooted to the Upham
        # rooting. Rooting matters: the tip-set lookup below is only valid when
        # this tree and the species tree agree on what a clade is.
        gtree = os.path.join(GENE_TREES, f"{g}.treefile")
        if not os.path.exists(gtree):
            print(f"{g:9s} missing gene tree")
            continue
        t = Tree(gtree, format=1)
        t.prune(sorted(seqs), preserve_branch_length=True)

        fg = bg = 0
        fg_tips = 0
        for node in t.traverse():
            if node.is_root():
                continue
            tips = frozenset(node.get_leaf_names())
            state = states.get(tips)
            if state is None:
                # The node is not a species-tree clade in its own right, which
                # happens when pruning collapses a clade whose other members
                # lack CDS. Fall back to the smallest species-tree clade that
                # contains exactly these tips among the retained taxa.
                cand = [(s, st) for s, st in states.items()
                        if (s & set(seqs)) == set(tips)]
                state = min(cand, key=lambda x: len(x[0]))[1] if cand else None
            if state == "diurnal":
                node.add_features(fg=True)
                fg += 1
                if node.is_leaf():
                    fg_tips += 1
            else:
                bg += 1

        # HyPhy reads {Foreground} appended to the node label.
        for node in t.traverse():
            if node.is_root():
                continue
            if getattr(node, "fg", False):
                node.name = f"{node.name}{{Foreground}}" if node.name else "{Foreground}"

        out_tree = os.path.join(TREE_OUT, f"{g}.labelled.nwk")
        t.write(outfile=out_tree, format=1, format_root_node=True)

        # Drop codon columns that are all gaps in the retained taxa. The codon
        # alignments were built across all 60 species, so a column occupied only
        # by species without CDS becomes empty once pruned: 13 of the 18 genes
        # have some, PER3 has 482 of 2122 and CSNK1D 239 of 792. HyPhy does not
        # merely ignore them, it dies on them with
        #   ASSERTION FAILED: Non-constant site passed to
        #   ComputeCompressedSubstitutionConstantSite
        # inside its constant-site optimisation. An all-gap column carries no
        # information, so removing it changes no result, but it DOES renumber
        # the sites, hence the map written alongside.
        ncod = width // 3
        keep = [c for c in range(ncod)
                if not all(s[c * 3:c * 3 + 3] == "---" for s in seqs.values())]
        dropped = ncod - len(keep)

        out_aln = os.path.join(CODON_OUT, f"{g}.codon.fasta")
        with open(out_aln, "w") as fh:
            for name in sorted(seqs):
                s = seqs[name]
                fh.write(f">{name}\n{''.join(s[c*3:c*3+3] for c in keep)}\n")

        # HyPhy site index (1-based, after cleaning) -> original codon index,
        # which is also the residue index in the untrimmed protein alignment.
        # Step 13 needs this to reach human residue numbering.
        with open(os.path.join(CODON_OUT, f"{g}.codonmap.csv"), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["hyphy_site", "original_codon"])
            for i, c in enumerate(keep, start=1):
                w.writerow([i, c + 1])

        flag = ""
        if fg < MIN_FOREGROUND or bg < MIN_BACKGROUND:
            flag = "LOW_POWER"
        qc.append(dict(gene=g, taxa=len(seqs), codons_in=ncod,
                       codons_used=len(keep), allgap_dropped=dropped,
                       foreground_branches=fg, background_branches=bg,
                       foreground_tips=fg_tips, flag=flag))
        print(f"{g:9s} {len(seqs):5d} {len(keep):7d} {dropped:8d} {fg:4d} {bg:4d} "
              f"{fg_tips:8d} {flag:>10s}")

    with open(os.path.join(PROJ, "results", "selection", "foreground_qc.csv"),
              "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(qc[0].keys()))
        w.writeheader()
        w.writerows(qc)
    low = [r["gene"] for r in qc if r["flag"]]
    print(f"\n{len(qc)} genes prepared; LOW_POWER: {', '.join(low) if low else 'none'}")


if __name__ == "__main__":
    main()
