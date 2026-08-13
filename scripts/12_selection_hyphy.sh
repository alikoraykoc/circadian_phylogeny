#!/usr/bin/env bash
# 12_selection_hyphy.sh
#
# Contrast-FEL: site-level test for DIFFERENT dN/dS between diurnal branches
#               (foreground) and the rest.
# RELAX:        gene-level test for RELAXED vs INTENSIFIED selection on the
#               foreground, reported as the selection intensity parameter K.
#
# This asks a different question from the convergence track. Convergence asks
# whether the SAME residue appeared repeatedly; selection asks whether the
# regime differs at all. A gene can sit under sharply different selection in
# diurnal lineages while showing no convergent residue whatsoever, so the
# convergence null constrains this result not at all.
#
# Per-gene trees are required, not one shared tree. Each gene has CDS for a
# different subset of species (32 to 58 of 60), so a single labelled tree cannot
# match every alignment. scripts/prep_selection_inputs.py builds them.
#
# Foreground is every DIURNAL branch (transition plus descendants that stay
# diurnal), because the hypothesis is a sustained regime rather than an event at
# the moment of transition.
set -uo pipefail
cd "$(dirname "$0")/.."
source config.sh

TREES="$RES/selection/trees"
CODON="$RES/codon"
OUT="$RES/selection"

if [[ ! -d "$TREES" ]]; then
  echo "No labelled trees. Run: python scripts/prep_selection_inputs.py" >&2
  exit 1
fi
mkdir -p "$OUT"

# HyPhy saturates the cores it is given and both analyses are per gene, so genes
# run sequentially rather than in parallel.
CPU="${HYPHY_CPU:-4}"

echo "== selection track == $(date)"
for g in $GENES; do
  aln="$CODON/$g.codon.fasta"
  tre="$TREES/$g.labelled.nwk"
  if [[ ! -f "$aln" || ! -f "$tre" ]]; then
    echo "SKIP $g: missing alignment or labelled tree"
    continue
  fi

  # Resumable: a finished JSON means the gene is done. HyPhy writes the JSON
  # only on success, so a partial run leaves nothing and is cleanly redone.
  cf="$OUT/$g.contrastfel.json"
  if [[ -s "$cf" ]]; then
    echo "SKIP $g contrast-fel (done)"
  else
    echo "== Contrast-FEL: $g == $(date)"
    hyphy CPU="$CPU" contrast-fel --alignment "$aln" --tree "$tre" \
          --branch-set Foreground --output "$cf" 2>&1 | tail -3
  fi

  rx="$OUT/$g.relax.json"
  if [[ -s "$rx" ]]; then
    echo "SKIP $g relax (done)"
  else
    echo "== RELAX: $g == $(date)"
    hyphy CPU="$CPU" relax --alignment "$aln" --tree "$tre" \
          --test Foreground --output "$rx" 2>&1 | tail -3
  fi
done

echo "== SELECTION COMPLETE == $(date)"
echo "JSON in $OUT; summarize with scripts/summarize_selection.py"
