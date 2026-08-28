#!/usr/bin/env bash
# Fast regression control: plant convergence in two small genes, run PCOC, score it.
#
# Purpose
# -------
# The full spike-in (spikein_run.sh) takes hours, so it gets run once at the end
# of a session, which is exactly when a defect introduced that morning is most
# expensive. This is the version you run after touching ANYTHING in the tree,
# scenario or alignment path. Target: about 5 minutes.
#
# Ten defects were found in this pipeline and not one raised an error. They were
# all interface errors between components, so a control has to exercise the
# assembled path rather than re-test the pieces. This does, at reduced scope.
#
# Scope reduction, and why it is safe
# -----------------------------------
# Two genes, 'all' difficulty only (every gain event converges). That is a
# pass/fail smoke test of the plumbing, NOT a measurement of sensitivity: it
# cannot tell you that PCOC misses partial convergence, because it never plants
# any. Use the full control for that.
#
# FBXL3 and BHLHE40 are chosen as the two smallest trimmed alignments (426 and
# 407 columns) that still have complete 60-taxon occupancy, so detection is fast
# and no gene-specific pruning path is exercised. If you change the pair, prefer
# full-occupancy genes for the same reason.
#
# Pass criterion: PCOC recovers the planted sites. Anything less means signal is
# dying between the alignment and the detector, and every null result in the
# project is void until it is located.
#
# Writes to results/spikein_quick/ so it can never touch the full control's
# alignments or answer key.
#
# Usage: bash scripts/spikein_quick.sh
# NOT set -u: conda's activate-gfortran script references an unbound GFORTRAN
# variable, which -u turns into a fatal error.
set -o pipefail
cd "$(dirname "$0")/.."
source config.sh

export SPIKE_OUT="$RES/spikein_quick"
export SPIKE_GENES="FBXL3 BHLHE40"
export SPIKE_LEVELS="all"
export SPIKE_SITES_PER_LEVEL="3"

if ! docker info > /dev/null 2>&1; then
  echo "FATAL: Docker daemon is not reachable. Start Docker Desktop, then re-run." >&2
  echo "  open -a Docker" >&2
  exit 1
fi

source ~/miniforge3/etc/profile.d/conda.sh 2>/dev/null || source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate phylo

# Always regenerate: the point is to test the CURRENT tree and scenario files.
rm -rf "$SPIKE_OUT"
mkdir -p "$SPIKE_OUT/pcoc"

echo "########## QUICK SPIKE-IN $(date) ##########"
python scripts/spikein_generate.py || exit 1

for g in $SPIKE_GENES; do
  SCEN="$(cat "$RES/pcoc/scenarios/ER_gain/$g.scenario")"
  [[ -z "$SCEN" ]] && { echo "FATAL: empty scenario for $g" >&2; exit 1; }
  echo "== PCOC: $g == $(date)"
  docker run --rm -v "$PROJ":/proj carinerey/pcoc pcoc_det.py \
    -t  "/proj/results/branchlengths/pcoc/$g.treefile" \
    -aa "/proj/results/spikein_quick/trim/$g.trim.fa" \
    -m  "$SCEN" \
    -o  "/proj/results/spikein_quick/pcoc/$g" \
    -f  0.0 -cpu 4 --gamma --no_cleanup > /dev/null 2>&1
  # A dead container is indistinguishable from a gene that ran unless you check.
  if ! compgen -G "$SPIKE_OUT/pcoc/$g/RUN_*/$g.trim.results.tsv" > /dev/null; then
    echo "FATAL: PCOC produced no results table for $g" >&2
    exit 1
  fi
done

echo
python scripts/spikein_evaluate.py
