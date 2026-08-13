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
mkdir -p "$OUT" "$OUT/logs"

# HyPhy saturates the cores it is given and both analyses are per gene, so genes
# run sequentially rather than in parallel.
CPU="${HYPHY_CPU:-4}"

# Two passes rather than both analyses per gene.
#
# Contrast-FEL is reliable here (every gene attempted has converged) and it is
# the analysis that matters most: it works site by site and feeds the step 13
# consensus. RELAX is unreliable on this data, roughly half its attempts fail on
# flat likelihood surfaces, and each attempt costs 30 to 45 minutes. Interleaving
# them lets a stuck RELAX block the site-level results for every later gene.
#
#   bash scripts/12_selection_hyphy.sh contrastfel   run the reliable pass first
#   bash scripts/12_selection_hyphy.sh relax         then best-effort RELAX
#   bash scripts/12_selection_hyphy.sh both          original interleaved behaviour
MODE="${1:-both}"
case "$MODE" in
  contrastfel|relax|both) ;;
  *) echo "usage: $0 [contrastfel|relax|both]" >&2; exit 1 ;;
esac

echo "== selection track ($MODE) == $(date)"
for g in $GENES; do
  aln="$CODON/$g.codon.fasta"
  tre="$TREES/$g.labelled.nwk"
  if [[ ! -f "$aln" || ! -f "$tre" ]]; then
    echo "SKIP $g: missing alignment or labelled tree"
    continue
  fi

  # Resumable: a finished JSON means the gene is done. HyPhy writes the JSON
  # only on success, so a partial run leaves nothing and is cleanly redone.
  # Full per-gene logs. An earlier version piped through `tail -3`, which threw
  # away the diagnostic for a RELAX failure and left only "Check errors.log",
  # by which time the next gene had already truncated errors.log. Keep
  # everything; these logs are small and are the only record of why a gene died.
  cf="$OUT/$g.contrastfel.json"
  if [[ "$MODE" == "relax" ]]; then
    :
  elif [[ -s "$cf" ]]; then
    echo "SKIP $g contrast-fel (done)"
  else
    echo "== Contrast-FEL: $g == $(date)"
    hyphy CPU="$CPU" contrast-fel --alignment "$aln" --tree "$tre" \
          --branch-set Foreground --output "$cf" > "$OUT/logs/$g.contrastfel.log" 2>&1 \
      || echo "  FAILED contrast-fel $g, see $OUT/logs/$g.contrastfel.log"
    tail -2 "$OUT/logs/$g.contrastfel.log"
  fi

  # RELAX fails intermittently on this data, and the SAME inputs can succeed on a
  # retry (NPAS2 failed twice then completed). The cause is its own diagnostic:
  # "Potential convergence issues due to flat likelihood surfaces". RELAX fits a
  # separate rate parameter per branch, so on trees carrying many near-zero
  # branches the surface for K is close to unidentifiable and the optimizer can
  # wander off. Retry rather than lose the gene, and treat any K that only
  # appears on some attempts as unreliable (see results/selection/replicates/).
  rx="$OUT/$g.relax.json"
  if [[ "$MODE" == "contrastfel" ]]; then
    :
  elif [[ -s "$rx" ]]; then
    echo "SKIP $g relax (done)"
  else
    for attempt in 1 2 3; do
      echo "== RELAX: $g (attempt $attempt) == $(date)"
      hyphy CPU="$CPU" relax --alignment "$aln" --tree "$tre" \
            --test Foreground --output "$rx" \
            > "$OUT/logs/$g.relax.attempt$attempt.log" 2>&1
      if [[ -s "$rx" ]]; then
        cp "$OUT/logs/$g.relax.attempt$attempt.log" "$OUT/logs/$g.relax.log"
        break
      fi
      echo "  attempt $attempt failed"
    done
    [[ -s "$rx" ]] || echo "  RELAX GAVE UP on $g after 3 attempts"
  fi
done

echo "== SELECTION COMPLETE == $(date)"
echo "JSON in $OUT; summarize with scripts/summarize_selection.py"
