#!/usr/bin/env bash
# 10_tdg09.sh
# TDG09: site-specific test for a shift in amino-acid fitness between two sets of
# clades. It is the ML sibling of PCOC's PC submodel, so it cross-validates PCOC.
# Needs: trimmed protein alignment, fixed-topology gene tree, and a grouping of
# tips into the two sets being contrasted (e.g. diurnal vs the rest).
set -euo pipefail
source "$(dirname "$0")/../config.sh"

TDG09_JAR="tools/tdg09.jar"          # EDIT: path to tdg09.jar (from tdg09 release)
GROUPS="$RES/scenario/tdg09_groups.txt"   # EDIT: produce this from the scenario

if [[ ! -f "$TDG09_JAR" ]]; then
  echo "tdg09.jar not found at $TDG09_JAR. Download from the tdg09 release and set path."
  exit 1
fi

for g in $GENES; do
  echo "== TDG09: $g =="
  java -cp "$TDG09_JAR" tdg09.Analyse \
    -alignment "$RES/trim/$g.trim.$AAEXT" \
    -tree      "$RES/branchlengths/$g.treefile" \
    -groups    "$GROUPS" \
    -threads   4 \
    > "$RES/tdg09/$g.tdg09.out"
done

echo "Per-site LRTs in results/tdg09/. Sites significant here AND in PCOC are strong."
