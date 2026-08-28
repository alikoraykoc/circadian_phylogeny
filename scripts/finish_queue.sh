#!/usr/bin/env bash
# Run the remaining analyses to completion, strictly one at a time.
#
# Order is chosen by what can invalidate already-written text, per the project
# plan. Each stage must succeed before the next starts: scoring a broken control
# or building on a superseded curve is how this project accumulated ten silent
# defects.
#
#   1. score the end-to-end spike-in control        (gates everything else)
#   2. re-run the partial-convergence k-curve       (published one is superseded)
#   3. ER_reversal sweep                            (never tested, in the premise)
#
# Deliberately NOT queued, and why:
#   ARD sensitivity sweeps   ~12 h, and ARD places a diurnal ancestral mammal,
#                            which contradicts the literature the study rests on.
#                            It is a robustness check, not a hypothesis test.
#   RELAX remainder          the test does not run reliably on this data; adding
#                            more non-converging genes does not fix that.
#   RELAX taxon-dropping     user deferred; and dropping taxa to make a test
#                            converge must be declared, not done quietly.
#
# NOT set -u: conda's activate-gfortran script references an unbound GFORTRAN.
set -o pipefail
cd "$(dirname "$0")/.."
source config.sh
source ~/miniforge3/etc/profile.d/conda.sh 2>/dev/null || source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate phylo

LOG="$RES/finish_queue.log"
say() { echo "##### $* $(date)" | tee -a "$LOG"; }

say "QUEUE START"

# ---- wait for any spike-in still running -----------------------------------
if pgrep -f "spikein_run.sh" > /dev/null; then
  say "waiting for spikein_run.sh to finish"
  while pgrep -f "spikein_run.sh" > /dev/null; do sleep 60; done
fi

# ---- gate: the control must actually be complete ---------------------------
NP=$(find "$RES/spikein/pcoc" -name "*.trim.results.tsv" 2>/dev/null | wc -l | tr -d ' ')
NT=$(ls -1 "$RES/spikein/tdg09"/*.tdg09.out 2>/dev/null | wc -l | tr -d ' ')
say "spike-in state: PCOC $NP/18, TDG09 $NT/18"
if [[ "$NP" -lt 18 ]]; then
  say "FATAL: spike-in PCOC incomplete ($NP/18). Not scoring a partial control."
  exit 1
fi

# ---- 1. score the control --------------------------------------------------
say "STAGE 1: score spike-in"
python scripts/spikein_evaluate.py 2>&1 | tee -a "$LOG" || { say "FATAL: scoring failed"; exit 1; }

# ---- 2. k-curve on corrected scenarios -------------------------------------
say "STAGE 2: partial-convergence k-curve, CLOCK"
python scripts/pcoc_partial_convergence_sweep.py CLOCK >> "$LOG" 2>&1 \
  || { say "FATAL: k-curve failed"; exit 1; }
say "k-curve done"

# ---- 3. reversal sweep -----------------------------------------------------
say "STAGE 3: ER_reversal PCOC sweep"
bash scripts/run_pcoc_all.sh ER_reversal >> "$LOG" 2>&1 \
  || { say "FATAL: reversal sweep failed"; exit 1; }
say "reversal sweep done"

say "QUEUE COMPLETE"
