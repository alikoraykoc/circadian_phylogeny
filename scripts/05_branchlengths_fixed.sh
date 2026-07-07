#!/usr/bin/env bash
# 05_branchlengths_fixed.sh
# THE REQUIRED TREE STEP. Two tree sets on the fixed Upham topology:
#
#   (A) PCOC / TDG09 trees: per-gene MFP. Independent per-gene analyses with no
#       cross-gene comparison, so each gene gets its best-fit branch lengths.
#
#   (B) RERconverge trees: one uniform model (from supermatrix MFP, step 05).
#       RER normalizes scale but not model-driven shape differences (long-vs-short
#       branch ratios from different rate-heterogeneity choices). A shared model
#       prevents that confound.
#
# Step 05 must run BEFORE this script so that best_model.txt exists.
set -euo pipefail
source "$(dirname "$0")/../config.sh"

# ---- resolve RER model ----
if [[ -f "$RER_MODEL_FILE" ]]; then
  RER_MODEL=$(cat "$RER_MODEL_FILE")
  echo "RER uniform model (from supermatrix): $RER_MODEL"
else
  echo "WARNING: $RER_MODEL_FILE not found. Run step 05 first."
  echo "Falling back to $RER_MODEL_FALLBACK"
  RER_MODEL="$RER_MODEL_FALLBACK"
fi

# ---- prepare output dirs ----
mkdir -p "$RES/branchlengths/pcoc" "$RES/branchlengths/rer"
: > "$RES/branchlengths/rer_input.trees"

for g in $GENES; do
  aln="$RES/trim/$g.trim.$AAEXT"

  # Prune species tree to this gene's taxa if needed
  aln_taxa=$(grep '^>' "$aln" | sed 's/>//' | sort)
  tree_taxa=$(sed 's/[,();]/ /g' "$SPTREE" | tr ' ' '\n' | grep -oE '^[A-Za-z][A-Za-z_]+' | sort -u)
  if [[ "$aln_taxa" == "$tree_taxa" ]]; then
    sp="$SPTREE"
  else
    sp="$RES/branchlengths/$g.pruned.nwk"
    python3 -c "
from ete3 import Tree
t = Tree('$SPTREE', format=1)
aln_tips = set(open('$aln').read().split('>')[1:])
aln_tips = {s.split('\n')[0].strip() for s in aln_tips}
t.prune(list(aln_tips & {l.name for l in t.get_leaves()}), preserve_branch_length=True)
t.write(outfile='$sp', format=5)
"
    echo "  pruned species tree for $g ($(echo "$aln_taxa" | wc -l | tr -d ' ') taxa)"
  fi

  # (A) PCOC / TDG09 tree: per-gene MFP
  echo "== PCOC/TDG09 tree (MFP): $g =="
  iqtree -s "$aln" -te "$sp" -m "$PCOC_MODEL" --keep-ident \
          -pre "$RES/branchlengths/pcoc/$g" -T "$THREADS" -redo

  # (B) RERconverge tree: uniform model
  echo "== RER tree ($RER_MODEL): $g =="
  iqtree -s "$aln" -te "$sp" -m "$RER_MODEL" --keep-ident \
          -pre "$RES/branchlengths/rer/$g" -T "$THREADS" -redo

  # collect RER trees (tab-separated: gene name, then the newick)
  printf '%s\t%s\n' "$g" "$(cat "$RES/branchlengths/rer/$g.treefile")" \
        >> "$RES/branchlengths/rer_input.trees"
done

echo
echo "PCOC/TDG09 trees: $RES/branchlengths/pcoc/"
echo "RER trees:        $RES/branchlengths/rer/"
echo "RER collection:   $RES/branchlengths/rer_input.trees"
echo
echo "SANITY: branch lengths should be small decimals (subs/site), NOT the"
echo "million-year scale of the Upham tree, and no branch should be 0."
echo "RE-ROOT each *.treefile to the Upham rooting before scenario building (07)."
