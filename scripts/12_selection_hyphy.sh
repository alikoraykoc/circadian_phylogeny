#!/usr/bin/env bash
# 12_selection_hyphy.sh   (needs codon alignments from step 02)
# Contrast-FEL: site-level test for DIFFERENT dN/dS between diurnal (foreground)
#               and the rest.
# RELAX:        gene-level test for RELAXED vs INTENSIFIED selection on foreground.
# Both need the codon alignment and a tree whose foreground branches are labelled.
set -euo pipefail
source "$(dirname "$0")/../config.sh"

LABELLED_TREE="$RES/scenario/species_labelled_foreground.nwk"   # EDIT: produced by labelling step

if [[ ! -f "$LABELLED_TREE" ]]; then
  echo "Foreground-labelled tree missing. Create it from step 07 outputs:"
  echo "  hyphy label-tree ... --label Foreground  (tag the transition branches)"
  exit 1
fi

for g in $GENES; do
  codon="$RES/codon/$g.codon.fasta"
  [[ -f "$codon" ]] || { echo "No codon alignment for $g, skipping."; continue; }

  echo "== Contrast-FEL: $g =="
  hyphy contrast-fel --alignment "$codon" --tree "$LABELLED_TREE" \
        --branch-set Foreground --output "$RES/selection/$g.contrastfel.json"

  echo "== RELAX: $g =="
  hyphy relax --alignment "$codon" --tree "$LABELLED_TREE" \
        --test Foreground --output "$RES/selection/$g.relax.json"
done

echo "Selection results (JSON) in results/selection/."
