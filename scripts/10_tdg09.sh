#!/usr/bin/env bash
# 10_tdg09.sh
# TDG09: site-specific test for a shift in amino-acid fitness between two sets of
# clades. It is the ML sibling of PCOC's PC submodel, so it cross-validates PCOC.
# Inputs: prefix-labelled PHYLIP alignments and trees from prep_tdg09_inputs.py.
set -euo pipefail
source "$(dirname "$0")/../config.sh"

TDG09_JAR="$PROJ/tdg09-1.1.2/dist/tdg09.jar"
INPUT_DIR="$RES/tdg09/inputs"

if [[ ! -f "$TDG09_JAR" ]]; then
  echo "tdg09.jar not found at $TDG09_JAR"
  exit 1
fi

for g in $GENES; do
  echo "== TDG09: $g =="
  java -cp "$TDG09_JAR" tdg09.Analyse \
    -alignment "$INPUT_DIR/$g.phy" \
    -tree      "$INPUT_DIR/$g.tree" \
    -groups    Di No \
    -threads   2 \
    > "$RES/tdg09/$g.tdg09.out" 2>&1
  echo "  done: $g"
done

echo "Per-site LRTs in results/tdg09/. Sites significant here AND in PCOC are strong."
