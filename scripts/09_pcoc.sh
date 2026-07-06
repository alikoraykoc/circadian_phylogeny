#!/usr/bin/env bash
# 09_pcoc.sh
# Convergence detection with PCOC (Docker image), per gene, plus a power/FPR
# simulation on your actual tree + scenario. Runs on the fixed-topology gene
# trees (04) and the trimmed protein alignments (01), using the scenario (07).
set -euo pipefail
source "$(dirname "$0")/../config.sh"

PCOC_IMG="carinerey/pcoc"           # docker pull carinerey/pcoc
SCEN="$(cat "$RES/scenario/pcoc_scenario.txt")"

run_pcoc () {  # $1 = docker mount root ; runs pcoc inside container
  docker run --rm -v "$PROJ":/proj "$PCOC_IMG" "$@"
}

# ---- power / FPR simulation FIRST (calibrate the posterior threshold) ----
# Use one representative gene tree; the scenario drives power more than the gene.
echo "== PCOC power simulation =="
run_pcoc pcoc_sim.py \
  -td /proj/results/branchlengths \
  -o  /proj/results/pcoc/sim \
  -m  "$SCEN"  || echo "Check pcoc_sim.py flags for your image version (-h)."

# ---- detection per gene ----
for g in $GENES; do
  echo "== PCOC detect: $g =="
  run_pcoc pcoc_det.py \
    -t  "/proj/results/branchlengths/$g.treefile" \
    -aa "/proj/results/trim/$g.trim.$AAEXT" \
    -m  "$SCEN" \
    -o  "/proj/results/pcoc/$g" \
    --plot  || echo "Check pcoc_det.py flags for your image version (-h)."
done

echo
echo "Per-site PCOC/PC/OC posteriors in results/pcoc/<gene>/"
echo "Apply the threshold chosen from the simulation, not a fixed 0.8 by default."
