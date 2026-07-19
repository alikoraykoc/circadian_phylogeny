#!/usr/bin/env bash
# Run PCOC detection over all genes for ONE scenario set, sequentially,
# 4 CPU threads per gene.
#
# Usage: bash scripts/run_pcoc_all.sh <SET>
#   where <SET> is a directory under results/pcoc/scenarios/, i.e. one
#   combination of ASR model and transition direction:
#
#     ER_gain       primary analysis: derived diurnality under the ER model
#     ER_reversal   reversals to nocturnality under the ER model
#     ARD_gain      sensitivity: derived diurnality under the ARD model
#     ARD_reversal  sensitivity: reversals under the ARD model
#
# Gains and reversals are run as SEPARATE convergent classes. PCOC fits a single
# convergent amino-acid profile per run, so branches moving toward opposite
# phenotypes must never share one class.
#
# -f 0.0 keeps the posterior for every site rather than only those above a fixed
# cutoff, so the significance threshold stays a post-hoc decision that can be
# calibrated with pcoc_sim without re-running the sweep (see CLAUDE.md).
set -euo pipefail
cd /Users/koray/Desktop/circadian_analysis
source config.sh

SET="${1:-}"
if [[ -z "$SET" ]]; then
  echo "usage: $0 <ER_gain|ER_reversal|ARD_gain|ARD_reversal>" >&2
  exit 1
fi

# Cooperative stop. run_pcoc_sets.sh invokes this script as a fresh bash process
# per set, so dropping this flag lets an in-flight set finish its current gene
# sweep untouched while preventing any LATER set from starting. Remove the file
# to resume. Killing the loop instead would abandon a gene mid-computation.
STOP_FLAG="$RES/pcoc/STOP"
if [[ -f "$STOP_FLAG" ]]; then
  echo "STOP flag present ($STOP_FLAG): not starting set $SET"
  echo "Remove the file to resume: rm $STOP_FLAG"
  exit 0
fi

PCOC_IMG="carinerey/pcoc"
SCEN_DIR="$RES/pcoc/scenarios/$SET"
OUT_ROOT="$RES/pcoc/$SET"

if [[ ! -d "$SCEN_DIR" ]]; then
  echo "No scenario dir: $SCEN_DIR (run prep_pcoc_scenarios.py first)" >&2
  exit 1
fi

mkdir -p "$OUT_ROOT"
echo "== PCOC set: $SET == $(date)"

for g in $GENES; do
  # PCOC writes $OUT_ROOT/$g/RUN_<timestamp>/$g.trim.results.tsv only on success,
  # so a partial run leaves no results file and is cleanly redone.
  if compgen -G "$OUT_ROOT/$g/RUN_*/$g.trim.results.tsv" > /dev/null; then
    echo "SKIP $g (already done)"
    continue
  fi

  SCEN="$(cat "$SCEN_DIR/$g.scenario")"
  if [[ -z "$SCEN" ]]; then
    echo "SKIP $g: empty scenario (no events survive in this gene tree)"
    continue
  fi

  echo "== PCOC detect: $g [$SET] == $(date)"
  docker run --rm -v "$PROJ":/proj "$PCOC_IMG" pcoc_det.py \
    -t  "/proj/results/branchlengths/pcoc/$g.treefile" \
    -aa "/proj/results/trim/$g.trim.$AAEXT" \
    -m  "$SCEN" \
    -o  "/proj/results/pcoc/$SET/$g" \
    -f  0.0 \
    -cpu 4 \
    --gamma \
    --no_cleanup 2>&1
  echo "  done: $g $(date)"
  echo
done

echo "== PCOC SET COMPLETE: $SET == $(date)"
