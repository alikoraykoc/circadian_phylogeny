#!/usr/bin/env bash
# Run PCOC sequentially, one gene at a time, 4 CPU threads per gene.
set -euo pipefail
cd /Users/koray/Desktop/circadian_analysis
source config.sh

PCOC_IMG="carinerey/pcoc"
SCEN_DIR="$RES/pcoc/scenarios"

for g in $GENES; do
  # PCOC writes results to $RES/pcoc/$g/RUN_<timestamp>/$g.trim.results.tsv,
  # so check for any such file (glob may match multiple runs; -n on the list).
  if compgen -G "$RES/pcoc/$g/RUN_*/$g.trim.results.tsv" > /dev/null; then
    echo "SKIP $g (already done)"
    continue
  fi
  SCEN="$(cat "$SCEN_DIR/$g.scenario")"
  if [[ -z "$SCEN" ]]; then
    echo "SKIP $g: no scenario"
    continue
  fi
  echo "== PCOC detect: $g == $(date)"
  docker run --rm -v "$PROJ":/proj "$PCOC_IMG" pcoc_det.py \
    -t  "/proj/results/branchlengths/pcoc/$g.treefile" \
    -aa "/proj/results/trim/$g.trim.$AAEXT" \
    -m  "$SCEN" \
    -o  "/proj/results/pcoc/$g" \
    -f  0.8 \
    -cpu 4 \
    --gamma \
    --no_cleanup 2>&1
  echo "  done: $g $(date)"
  echo
done
echo "== ALL PCOC COMPLETE == $(date)"
