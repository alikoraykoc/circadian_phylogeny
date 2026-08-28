# Resume here

Last updated 2026-08-28. Everything below is committed.

## What is running

`scripts/finish_queue.sh` runs the remaining analyses strictly one at a time and
stops rather than building on a bad stage. Watch it with:

```
tail -f results/finish_queue.log
```

| stage | what | status |
|---|---|---|
| 0 | spike-in control, `scripts/spikein_run.sh` | PCOC in progress; TDG09 18/18 done |
| 1 | score the control, `spikein_evaluate.py` | queued |
| 2 | re-run the k-curve on corrected scenarios | queued |
| 3 | ER_reversal PCOC sweep | queued |

If the machine is restarted, relaunch with:

```
nohup caffeinate -i bash scripts/spikein_run.sh >> results/spikein/run.log 2>&1 &
nohup caffeinate -i bash scripts/finish_queue.sh > /dev/null 2>&1 &
```

Both skip completed genes, so restarting costs nothing.

**Docker Desktop must be running.** If it is not, PCOC's container exits
instantly, and because its output goes to `/dev/null` every gene looks like it
succeeded in under a second while producing nothing. That happened on 2026-08-28
and silently skipped 12 of 18 genes. Both scripts now preflight the daemon and
verify a results table per gene, so this fails loudly now, but start Docker
first: `open -a Docker`.

**Do not edit a shell script while it is running.** Bash re-reads the file as it
executes, so an edit shifts byte offsets and corrupts the running instance. This
killed the TDG09 leg mid-run on 2026-08-28 (the file itself was valid).

## The one thing that matters

**Until the control reports, the STRENGTH of the null result is not
established.** Its direction is well supported by six methods; how large a
signal we could have detected is not.

Do not read `results/spikein/ANSWER_KEY.csv`. `spikein_evaluate.py` is the only
thing that should open it.

### How to read the outcome

- **PCOC recovers most of the 36 `all` sites** means the plumbing carries signal
  end to end and the null is interpretable.
- **PCOC recovers few or none** means signal is dying somewhere between alignment
  and detector, and every null result is void until it is located.
- **The `few` level (3 of 10 lineages)** is the important one for the parsimony
  test, which is the primary evidence precisely because it should not need
  near-universal convergence.
- **TDG09's false positive count** finally adjudicates it. On real data it flags
  885 of 3,820 testable sites while nothing corroborates any of them.

The first version of this control was faulty (it planted into 41 of 60 species,
making the residue a new consensus rather than convergence) and all conclusions
from it were withdrawn. The current version plants into day-active species only,
verified against the sealed key: mean diurnal gap 1.000 / 0.816 / 0.374 by
level, zero nocturnal carriers. See commit `d1fb5c7`.

## Documents

| file | what it is |
|---|---|
| `RESULTS.md` | the results section. Section 11 is marked RUNNING pending the control |
| `METHODS.md` | full methods, with software versions |
| `PROGRESS_REPORT.md` | section 0A is the read-first technical log |
| `EXPLAINER.html` | plain-language account for a biologist |
| `shared_results/` | the actual result files, since `results/` is gitignored |

## Known-stale text to fix when the queue finishes

1. **`RESULTS.md` section 4, the k-curve.** Computed 2026-08-13 20:47; the
   re-rooting fix landed 21:28 the same evening, so it rests on superseded
   scenario files. A provenance caveat is in place. When stage 2 finishes,
   compare against the old curve and say explicitly whether the shape changed
   rather than silently replacing numbers. `PROGRESS_REPORT.md` 0A.5b needs the
   same treatment.
2. **`RESULTS.md` section 11** is a placeholder until stage 1 reports.

## Deliberately not queued, and why

- **ARD sensitivity sweeps** (~12 h). ARD places a diurnal ancestral placental
  mammal, contradicting the literature the study rests on; it is a robustness
  check on a modelling choice, not a hypothesis test. Calibrate before sweeping
  if it is ever run.
- **RELAX remainder.** The test does not run reliably on this data. 9 of 18
  converged, 5 gave up after 3 attempts each, all failing in ancestral
  reconstruction on deeply divergent lineages (`Lagorchestes_hirsutus` 9x,
  `Phalanger_gymnotis` 3x, `Choloepus_hoffmanni` 3x, `Dasypus_novemcinctus` 1x).
  Adding more non-converging genes does not fix that. Its one significant result,
  NPAS2 K = 0.352, replicated 1 time in 11 attempts and is an artefact.
- **RELAX taxon-dropping experiment.** User deferred. If dropping those four
  lineages does let RELAX converge, **the removal must be declared in the
  methods, not done quietly**, and the result reported as conditional on it.

## Standing context

- **TDG09 disagrees with everything else** and the consensus rule handles it:
  885 sites flagged, 0 corroborated, consensus set empty. Never report its hits
  without that context.
- **Two genes cannot answer the question.** CSNK1D has no position at which two
  or more independent lineages changed; ARNTL has two. Their nulls are
  uninformative rather than negative.
- **Ten defects were found, none of which raised an error.** All were unchecked
  assumptions that two representations of the same object agreed.
  `scripts/lib_checks.py` asserts those boundaries; `scripts/test_lib_checks.py`
  proves each check fires on the real historical bug (15/15 pass). Run it after
  any refactor, and run `scripts/spikein_quick.sh` (about 5 minutes) after
  touching anything in the tree, scenario or alignment path.
