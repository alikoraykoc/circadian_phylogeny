#!/usr/bin/env bash
# 00_check_inputs.sh
# Verify inputs exist and that alignment tips are consistent with the species tree.
set -euo pipefail
source "$(dirname "$0")/../config.sh"

echo "== Checking inputs =="
[[ -f "$SPTREE" ]] || { echo "MISSING species tree: $SPTREE"; exit 1; }
[[ -f "$DIEL"   ]] || { echo "MISSING diel activity table: $DIEL"; exit 1; }

missing=0
for g in $GENES; do
  f="$ALN/$g.$AAEXT"
  if [[ ! -f "$f" ]]; then echo "MISSING alignment: $f"; missing=1; fi
done
[[ $missing -eq 0 ]] || { echo "Fix missing alignments before continuing."; exit 1; }

echo "== Tips present in species tree =="
# crude tip extraction from a Newick file
grep -oE '[A-Za-z0-9_\.]+' "$SPTREE" | sort -u > "$RES/_tree_tips.txt" || true
echo "Tree tip tokens written to $RES/_tree_tips.txt (inspect manually)."

echo
echo "NEXT: for each gene, confirm that alignment headers match tree tips exactly."
echo "Mismatched labels are the most common cause of IQ-TREE -te failures."
echo "All inputs present. OK."
