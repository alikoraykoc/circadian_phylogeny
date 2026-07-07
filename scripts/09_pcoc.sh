#!/usr/bin/env bash
# 09_pcoc.sh
# Convergence detection with PCOC (Docker image), per gene.
# Uses per-gene scenario strings (from prep_pcoc_scenarios.py) that map
# species-tree transition branches onto each gene tree's PCOC node numbering.
set -euo pipefail
source "$(dirname "$0")/../config.sh"

PCOC_IMG="carinerey/pcoc"
SCEN_DIR="$RES/pcoc/scenarios"

if [[ ! -d "$SCEN_DIR" ]]; then
  echo "Run prep_pcoc_scenarios.py first to generate per-gene scenarios."
  exit 1
fi

# ---- detection per gene ----
for g in $GENES; do
  SCEN="$(cat "$SCEN_DIR/$g.scenario")"
  if [[ -z "$SCEN" ]]; then
    echo "SKIP $g: no scenario"
    continue
  fi
  echo "== PCOC detect: $g =="
  docker run --rm -v "$PROJ":/proj "$PCOC_IMG" pcoc_det.py \
    -t  "/proj/results/branchlengths/pcoc/$g.treefile" \
    -aa "/proj/results/trim/$g.trim.$AAEXT" \
    -m  "$SCEN" \
    -o  "/proj/results/pcoc/$g" \
    -f  0.8 \
    --gamma \
    --no_cleanup \
    2>&1 | tail -5
  echo "  done: $g"
done

echo
echo "Per-site PCOC/PC/OC posteriors in results/pcoc/<gene>/"
