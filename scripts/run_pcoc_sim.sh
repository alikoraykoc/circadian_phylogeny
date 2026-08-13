#!/usr/bin/env bash
# Power and false-positive calibration for the PCOC detection sweep.
#
# Usage: bash scripts/run_pcoc_sim.sh [SET]      (default: ER_gain)
#
# Why this exists
# ---------------
# CLAUDE.md requires the PCOC posterior threshold to come from pcoc_sim
# power/FPR calibration on the ACTUAL tree and scenario, not from a fixed 0.8.
# It is also the only way to tell a biological null from an underpowered one:
# the ER_gain sweep returned zero sites, and that number means nothing until we
# know the design could have detected a planted signal.
#
# Each gene is calibrated on its OWN tree, because mean branch length spans a
# 37-fold range across the 18 genes (ARNTL 0.0035 to CSNK1E 0.130 subs/site).
# Power depends on how much evolutionary time the convergent branches have, so a
# single pooled calibration would be meaningless.
#
# The -c flag is load-bearing
# ---------------------------
# pcoc_sim reads -m only to define the POOL of candidate events, then simulates
# n_events of them, where n_events = random.randrange(c_min, c_max+1) whenever
# -c is left at its default of 0 (pcoc_sim.py:496). A pilot run on RORB silently
# simulated 7 of the 10 real events, discarding 3 at random
# (events_placing.py:placeNTransitionsInTree deletes c_max - c events). Passing
# -c equal to the gene's real event count is what forces the simulated scenario
# to match the detection scenario branch for branch.
set -euo pipefail
cd /Users/koray/Desktop/circadian_analysis
source config.sh

SET="${1:-ER_gain}"
PCOC_IMG="carinerey/pcoc"
SCEN_DIR="$RES/pcoc/scenarios/$SET"
SIM_ROOT="$RES/pcoc_sim/$SET"
TREE_ROOT="$RES/pcoc_sim/trees"

# Sites per simulated alignment, and number of ancestral/convergent profile
# couples drawn from C60 per gene. The couples are the real variance here: power
# depends on how far apart the two profiles are, so a spread of couples gives a
# power curve rather than a single point estimate.
N_SITES=100
N_COUPLES=10

mkdir -p "$SIM_ROOT" "$TREE_ROOT"
echo "== PCOC SIM CALIBRATION: $SET == $(date)"

for g in $GENES; do
  if [[ -f "$SIM_ROOT/$g/DONE" ]]; then
    echo "SKIP $g (already calibrated)"
    continue
  fi

  SCEN="$(cat "$SCEN_DIR/$g.scenario")"
  if [[ -z "$SCEN" ]]; then
    echo "SKIP $g: empty scenario"
    continue
  fi

  # Event count is the number of "/"-separated groups in the scenario.
  NEV=$(awk -F'/' '{print NF}' <<< "$SCEN")

  # pcoc_sim takes a DIRECTORY of trees, so each gene needs its own single-tree
  # directory rather than a file path.
  mkdir -p "$TREE_ROOT/$g"
  cp "$RES/branchlengths/pcoc/$g.treefile" "$TREE_ROOT/$g/"

  echo "== sim: $g [$NEV events, $N_COUPLES couples] == $(date)"
  docker run --rm -v "$PROJ":/proj "$PCOC_IMG" pcoc_sim.py \
    -td "/proj/results/pcoc_sim/trees/$g" \
    -o  "/proj/results/pcoc_sim/$SET/$g" \
    -m  "$SCEN" \
    -c  "$NEV" \
    -n_sites "$N_SITES" \
    -nb_sampled_couple "$N_COUPLES" \
    --pcoc \
    -cpu 4 2>&1 | tail -3

  # Only mark DONE once the benchmark table actually exists, so an interrupted
  # gene is cleanly redone rather than skipped.
  if compgen -G "$SIM_ROOT/$g/RUN_*/Tree_1/BenchmarkResults.tsv" > /dev/null; then
    touch "$SIM_ROOT/$g/DONE"
    echo "  done: $g $(date)"
  else
    echo "  FAILED: $g produced no BenchmarkResults.tsv" >&2
  fi
  echo
done

echo "== PCOC SIM COMPLETE: $SET == $(date)"
