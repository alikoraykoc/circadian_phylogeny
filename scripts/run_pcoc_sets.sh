#!/usr/bin/env bash
# Run every PCOC scenario set in priority order.
#
# ER_gain is the primary analysis and runs first, so the headline result is
# available before the sensitivity sweeps finish. Each set is resumable: genes
# that already have a results.tsv are skipped, so an interrupted run picks up
# where it stopped.
set -uo pipefail
cd /Users/koray/Desktop/circadian_analysis

SETS=(ER_gain ER_reversal ARD_gain ARD_reversal)

for s in "${SETS[@]}"; do
  echo "############ START SET $s $(date) ############"
  if bash scripts/run_pcoc_all.sh "$s"; then
    echo "############ SET OK: $s $(date) ############"
  else
    # Do not let one bad set abort the rest; report and continue.
    echo "############ SET FAILED: $s (exit $?) $(date) ############"
  fi
  echo
done

echo "############ ALL PCOC SETS COMPLETE $(date) ############"
