#!/usr/bin/env bash
# 04_branchlengths_fixed.sh
# THE REQUIRED TREE STEP. For each gene, fix the Upham topology (-te) and
# re-estimate ONLY branch lengths + model parameters under one shared model.
# Output = species topology carrying gene-specific branch lengths (subs/site).
# This is what PCOC, TDG09, and (collected) RERconverge consume.
set -euo pipefail
source "$(dirname "$0")/../config.sh"

: > "$RES/branchlengths/rer_input.trees"   # RERconverge collection: "<gene>\t<newick>"

for g in $GENES; do
  aln="$RES/trim/$g.trim.$AAEXT"

  # If a gene is missing taxa, prune the species tree to that gene's tips first.
  # (See prune_tree.py helper; here we assume full taxon overlap.)
  sp="$SPTREE"
  # EDIT: sp="$RES/branchlengths/$g.sp.nwk"   # if you pruned per gene

  echo "== fixed-topology BL: $g =="
  iqtree2 -s "$aln" -te "$sp" -m "$FIXED_MODEL" \
          -pre "$RES/branchlengths/$g" -T "$THREADS" -redo
  # No -B: with a fixed topology, branch supports are meaningless.

  # collect for RERconverge (tab-separated: gene name, then the newick)
  printf '%s\t%s\n' "$g" "$(cat "$RES/branchlengths/$g.treefile")" \
        >> "$RES/branchlengths/rer_input.trees"
done

echo
echo "Per-gene fixed-topology trees in $RES/branchlengths/"
echo "SANITY: branch lengths should be small decimals (subs/site), NOT the"
echo "million-year scale of the Upham tree, and no branch should be 0."
echo "RE-ROOT each *.treefile to the Upham rooting before scenario building (07)."
