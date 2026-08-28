#!/usr/bin/env bash
# Run the pipeline against the spiked alignments, blind.
#
# Everything here uses the SAME scripts, trees and scenario files as the real
# analysis. Only the alignment changes. That is the point: the nine defects
# found so far were all in the plumbing between components, so the control has
# to exercise the plumbing rather than re-test the components.
#
# Deliberately unchanged and shared with the real run:
#   results/branchlengths/pcoc/<gene>.treefile        re-rooted gene trees
#   results/pcoc/scenarios/ER_gain/<gene>.scenario    the convergent branch sets
#   data/diel_activity.csv                            tip phenotypes
#
# Usage: bash scripts/spikein_run.sh
# NOT set -u: conda's activate-gfortran script references an unbound GFORTRAN
# variable, which -u turns into a fatal error. That killed this script the first
# time, after PCOC finished but before TDG09 started.
set -o pipefail
cd "$(dirname "$0")/.."
source config.sh

SPIKE="$RES/spikein"
ALN="$SPIKE/trim"

if [[ ! -d "$ALN" ]]; then
  echo "No spiked alignments. Run: python scripts/spikein_generate.py" >&2
  exit 1
fi
mkdir -p "$SPIKE/pcoc" "$SPIKE/tdg09"

echo "########## SPIKE-IN RUN $(date) ##########"

# Preflight: PCOC runs in Docker and its output is sent to /dev/null, so a dead
# daemon makes every gene "succeed" in under a second while producing nothing.
# That happened on 2026-08-28 and silently skipped 12 of 18 genes. Fail here
# instead, and check each gene's output below.
if ! docker info > /dev/null 2>&1; then
  echo "FATAL: Docker daemon is not reachable. Start Docker Desktop, then re-run." >&2
  echo "  open -a Docker" >&2
  exit 1
fi

# ---- 1. PCOC, same tree and scenario as the real sweep ---------------------
echo "### PCOC $(date)"
for g in $GENES; do
  if compgen -G "$SPIKE/pcoc/$g/RUN_*/$g.trim.results.tsv" > /dev/null; then
    echo "SKIP $g (done)"; continue
  fi
  SCEN="$(cat "$RES/pcoc/scenarios/ER_gain/$g.scenario")"
  [[ -z "$SCEN" ]] && { echo "SKIP $g: empty scenario"; continue; }
  echo "== PCOC: $g == $(date)"
  docker run --rm -v "$PROJ":/proj carinerey/pcoc pcoc_det.py \
    -t  "/proj/results/branchlengths/pcoc/$g.treefile" \
    -aa "/proj/results/spikein/trim/$g.trim.fa" \
    -m  "$SCEN" \
    -o  "/proj/results/spikein/pcoc/$g" \
    -f  0.0 -cpu 4 --gamma --no_cleanup > /dev/null 2>&1
  # Verify PCOC actually wrote a results table. Without this, a container that
  # fails to start is indistinguishable from a gene that ran.
  if ! compgen -G "$SPIKE/pcoc/$g/RUN_*/$g.trim.results.tsv" > /dev/null; then
    echo "FATAL: PCOC produced no results table for $g" >&2
    exit 1
  fi
  echo "  done: $g $(date)"
done

# ---- 2. TDG09, inputs rebuilt from the spiked alignments -------------------
echo "### TDG09 $(date)"
source ~/miniforge3/etc/profile.d/conda.sh 2>/dev/null || source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate phylo
TRIM_DIR="$ALN" RES_OVERRIDE="$SPIKE" PROJ="$PROJ" \
  python scripts/spikein_prep_tdg09.py

for g in $GENES; do
  [[ -s "$SPIKE/tdg09/$g.tdg09.out" ]] && { echo "SKIP $g (done)"; continue; }
  echo "== TDG09: $g == $(date)"
  java -cp "$PROJ/tdg09-1.1.2/dist/tdg09.jar" tdg09.Analyse \
    -alignment "$SPIKE/tdg09/inputs/$g.phy" \
    -tree      "$SPIKE/tdg09/inputs/$g.tree" \
    -groups    No Di \
    -threads   4 \
    > "$SPIKE/tdg09/$g.tdg09.out" 2>&1
  echo "  done: $g"
done

# ---- 3. model-free methods, via the path overrides -------------------------
echo "### parsimony $(date)"
TRIM_DIR="$ALN" OUT_DIR="$SPIKE" python scripts/observed_convergent_substitutions.py \
  > "$SPIKE/parsimony.log" 2>&1
tail -6 "$SPIKE/parsimony.log"

echo "### residue screen $(date)"
TRIM_DIR="$ALN" OUT_DIR="$SPIKE" python scripts/model_free_convergence_screen.py \
  > "$SPIKE/residue_screen.log" 2>&1
tail -3 "$SPIKE/residue_screen.log"

echo "########## SPIKE-IN RUN COMPLETE $(date) ##########"
echo "score it with: python scripts/spikein_evaluate.py"
