# Resume here

Last updated 2026-08-29. **All planned analyses are complete and pushed.**

## State

Nothing is running. The working tree is clean and `main` is level with `origin`.

| analysis | status |
|---|---|
| ER_gain PCOC sweep + calibration | done, null |
| ER_reversal PCOC sweep + calibration | done, 1 hit, artefact |
| ARD_gain PCOC sweep + calibration | done, 2 hits, same artefact |
| Partial-convergence k-curve | re-run on corrected scenarios, curve confirmed |
| End-to-end spike-in control | **passed**, 90 sites, PCOC 35/36 with 0 false positives |
| TDG09, RERconverge, Contrast-FEL | done, null |
| RELAX | not estimable on this data, closed |
| Boundary check suite | 19/19 |

## Read these first

| file | what it is |
|---|---|
| `RESULTS.md` | the results section, complete |
| `METHODS.md` | full methods with software versions |
| `PROGRESS_REPORT.md` | section 0A is the technical log |
| `EXPLAINER.html` | plain-language account for a biologist |
| `shared_results/` | the actual result tables, since `results/` is gitignored |

## The one thing to know before touching anything

**Ten defects were found in this pipeline and not one raised an error.** Every one
produced a confident, plausible, wrong number. They were all the same shape: an
unchecked assumption that two representations of the same object agreed.

`scripts/lib_checks.py` asserts those boundaries, and
`scripts/test_lib_checks.py` proves each check fires on the real historical bug
(19/19). Run it after any refactor. Run `scripts/spikein_quick.sh` (about 5
minutes) after touching anything in the tree, scenario or alignment path.

Two operational rules learned the hard way:

- **Start Docker Desktop before any PCOC run.** A dead daemon made all 18 genes
  "succeed" in under a second while producing nothing, because the container's
  output goes to `/dev/null`. Both runners now preflight the daemon and verify a
  results table per gene, but start it first: `open -a Docker`.
- **Never edit a shell script while it is running.** Bash re-reads the file as it
  executes, so an edit shifts byte offsets and corrupts the running instance.
  This killed a TDG09 leg mid-run.

## How to read a PCOC hit

The study produced three sites above threshold, all in RORB, all artefacts. The
diagnostic is `check_pcoc_hit_credible` in `lib_checks.py`, and the load-bearing
criterion is the **profile-change component**, not the combined posterior.

All three had `PC = 0.500` exactly, meaning the posterior rested entirely on the
one-change term, which asks only whether substitutions occurred on the declared
branches. RORB trimmed column 1 maps to untrimmed column 204, where 20 of 60
species' annotated proteins simply begin, so it carries 14 residues while its
neighbours carry one. Occupancy is not homology.

One of the three (ARD_gain site 2) would have passed the residue-count and
diagnostic-gap tests. Only the PC test caught it.

## What was deliberately not done, and why

- **ARD_reversal.** Scenarios are built (18 files, 10 events / 44 branches on
  CLOCK); add it to the `SETS` line in `scripts/ard_queue.sh` to run it. Under a
  reconstruction that already places a diurnal ancestral mammal, the reversal
  direction tests a scenario the study argues against twice over.
- **RELAX taxon-dropping.** Closed. Dropping the four lineages RELAX crashes on
  would select on the outcome rather than the biology; the nine genes that
  converged are already the subset that happened to converge.

## Known limits to state in any write-up

- **Two genes cannot answer the question.** CSNK1D has no position where two or
  more independent lineages changed; ARNTL has two. Uninformative, not negative.
- **The ARD sensitivity analysis is weaker than ER.** It declares 55 of 60 leaves
  convergent, leaving thin ancestral contrast. Consistent with the ER result
  rather than independent confirmation of it.
- **The 30/30 diurnal-nocturnal sample is artificial.** Real mammals are
  predominantly nocturnal. This is what pulls ARD to a diurnal ancestor.
- **TDG09 flags 23 percent of testable sites** and the spike-in measured its false
  positive rate at 24.7 percent. Never report its hits without that context.
