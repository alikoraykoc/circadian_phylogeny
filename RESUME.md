# Resume here

Paused 2026-08-14. Everything is committed and pushed; the working tree is clean.

## The one thing that matters

**The end-to-end validation control is half finished, and until it reports, the
STRENGTH of the null result is not established.** Its direction is well supported
by six methods; how large a signal we could have detected is not.

Resume it with:

```
bash scripts/spikein_run.sh          # skips the 7 genes already done
python scripts/spikein_evaluate.py   # scores against the sealed answer key
```

PCOC is 7 of 18 genes in. It is slow on spiked alignments, roughly 25 to 50
minutes per gene, because planted signal gives the optimiser something to fit
rather than being rejected immediately. Budget several hours, then TDG09 (~20
min) and the model-free methods (~15 min) follow automatically.

Do not read `results/spikein/ANSWER_KEY.csv` before running the evaluation.
`spikein_evaluate.py` is the only thing that should open it.

### What the control will tell you

- **PCOC recovers most of the 36 `all` level sites** &rarr; the plumbing carries
  signal end to end and the null is interpretable.
- **PCOC recovers few or none** &rarr; signal is dying somewhere between alignment
  and detector, and every null result is void until it is located.
- **The `few` level (3 of 10 lineages)** is the important one for the parsimony
  test, which is the primary evidence precisely because it should not need
  near-universal convergence.
- **TDG09's false positive count** finally adjudicates it. On real data it flags
  885 of 3820 sites while nothing corroborates any of them.

Note the first version of this control was faulty (it planted into 41 of 60
species, making the residue a new consensus rather than convergence) and all
conclusions from it were withdrawn. The current version plants into day-active
species only, verified: gap 1.000 / 0.800 / 0.346 by level, zero nocturnal
carriers. See commit `d1fb5c7`.

## Then, in priority order

1. **Re-run the k-curve** (`scripts/pcoc_partial_convergence_sweep.py CLOCK`).
   The published curve in PROGRESS_REPORT.md 0A.5b was computed 2026-08-13 20:07,
   BEFORE the rooting fix in `72b42b2`, and its `blocks.csv` still shows the old
   25-branch event structure. It underpins the claim that PCOC only detects
   near-universal convergence, which frames how the entire PCOC null is read.
   This can invalidate text that is already written, so it goes first.

2. **ER_reversal sweep** (`bash scripts/run_pcoc_all.sh ER_reversal`). Calibration
   is done and gives power 0.996 to 1.000 across all 18 genes, so the sweep is
   interpretable rather than uninformative. Closes a gap open since the start:
   reversions to nocturnality are part of the premise and have never been tested.
   Report CSNK1D separately, it retains only 2 reversal events against 5 elsewhere.

3. **Fast regression control** (plan item 2C, not yet built). Reduce the spike-in
   to 2 genes and 6 sites, about 5 minutes, and run it after any change to the
   tree, scenario or alignment path. This is what would have caught most of the
   ten defects at the moment they were introduced.

4. **Extend the k-curve** to CSNK1E (longest branches) and ARNTL (shortest), and
   raise replicates near k=8 to 10 where it is noisy.

5. **Deferrable:** ARD sensitivity sweeps (calibration first), the RELAX taxon
   dropping experiment (see below), finishing the RELAX remainder at 9 of 18.

## Standing context you will want

- **RELAX is not reliably estimable on this data.** 9 of 18 converged, 5 gave up
  after 3 attempts. It crashes in ancestral reconstruction naming deeply
  divergent lineages (`Lagorchestes_hirsutus` 9x, `Phalanger_gymnotis` 3x,
  `Choloepus_hoffmanni` 3x, `Dasypus_novemcinctus` 1x). The parked experiment is
  to test on one gene whether dropping those four lets it converge. **If it
  works, the taxon removal must be declared in the methods, not done quietly.**
  Its one significant result, NPAS2 K = 0.352 at p < 0.0001, replicated 1 time in
  11 attempts and is an artefact.

- **TDG09 disagrees with everything else** and the consensus rule handles it:
  885 sites flagged, 0 corroborated. Do not report its hits without that context.

- **Two genes cannot answer the question.** CSNK1D has zero positions with two or
  more independent changes; ARNTL has two. Their nulls are uninformative rather
  than negative.

- **Ten defects were found, none of which raised an error.** All were unchecked
  assumptions that two representations of the same object agreed. `lib_checks.py`
  now asserts those boundaries and `test_lib_checks.py` proves each check fires on
  the real historical bug. Run it after any refactor.

## Documents

- `PROGRESS_REPORT.md` section 0A is the read-first technical log.
- `EXPLAINER.html` is the plain-language account for a biologist.
- `shared_results/` holds the actual result files, since `results/` is gitignored.
