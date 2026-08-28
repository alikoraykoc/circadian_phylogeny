#!/usr/bin/env bash
# ARD sensitivity sweeps, calibration first, run after finish_queue.sh completes.
#
# Why this is a separate script rather than extra stages in finish_queue.sh:
# bash re-reads a script file as it executes, so editing a running script shifts
# byte offsets and corrupts the running instance. That killed the TDG09 leg on
# 2026-08-28. finish_queue.sh is running, so ARD gets its own file that waits.
#
# What this tests, and what it does not
# ------------------------------------
# ARD is the sensitivity analysis for the ancestral-state model choice, NOT a
# second hypothesis test. It reconstructs a DIURNAL ancestral placental mammal,
# rooting diurnality at a 54-tip clade, and reframes the history as 6 gains and
# 10 reversals against ER's 10 and 5. That contradicts the nocturnal-bottleneck
# literature the study rests on, and is the artefact the balanced 30/30 taxon
# sample was expected to produce.
#
# So a null under ARD adds robustness. A POSITIVE under ARD would need to be
# reported as conditional on a reconstruction the study otherwise argues against,
# not as a finding. Say so in any write-up.
#
# Calibration comes first for each direction because ARD's event counts differ
# enough from ER's to change power, and CLAUDE.md forbids assuming a threshold.
#
# NOT set -u: conda's activate-gfortran script references an unbound GFORTRAN.
set -o pipefail
cd "$(dirname "$0")/.."
source config.sh
source ~/miniforge3/etc/profile.d/conda.sh 2>/dev/null || source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate phylo

LOG="$RES/ard_queue.log"
say() { echo "##### $* $(date)" | tee -a "$LOG"; }

say "ARD QUEUE START"

if ! docker info > /dev/null 2>&1; then
  say "FATAL: Docker daemon not reachable. Start Docker Desktop, then re-run."
  exit 1
fi

# ---- wait for the ER queue, so only one analysis runs at a time -------------
if pgrep -f "finish_queue.sh" > /dev/null; then
  say "waiting for finish_queue.sh to finish"
  while pgrep -f "finish_queue.sh" > /dev/null; do sleep 60; done
fi
# It may have exited on a FATAL. Do not pile ARD on top of a failed queue.
if ! grep -q "QUEUE COMPLETE" "$RES/finish_queue.log" 2>/dev/null; then
  say "FATAL: finish_queue.sh did not reach QUEUE COMPLETE. Fix that first."
  exit 1
fi

for SET in ARD_gain ARD_reversal; do
  say "CALIBRATE $SET"
  bash scripts/run_pcoc_sim.sh "$SET" >> "$LOG" 2>&1 \
    || { say "FATAL: calibration failed for $SET"; exit 1; }
  python scripts/summarize_pcoc_sim.py "$SET" 2>&1 | tee -a "$LOG" \
    || { say "FATAL: calibration summary failed for $SET"; exit 1; }

  say "SWEEP $SET"
  bash scripts/run_pcoc_all.sh "$SET" >> "$LOG" 2>&1 \
    || { say "FATAL: sweep failed for $SET"; exit 1; }
  say "$SET done"
done

say "ARD QUEUE COMPLETE"
