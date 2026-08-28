# Results

Convergent molecular evolution in mammalian circadian genes and independent
transitions in diel activity.

Methods, parameters and software versions are in `METHODS.md`. Section numbers
below refer to it. Result tables cited by filename live in `shared_results/`.

---

## 1. Headline

**There is no detectable convergent amino acid signal associated with
independent gains of diurnality in any of the 18 circadian genes tested.**

Six independent lines of evidence, resting on different assumptions, agree:

| line of evidence | what it assumes | result |
|---|---|---|
| Direct parsimony substitution count | nothing | same-residue rate 0.310 against a null of 0.337; 0 phenotype-specific sites |
| Model-free diagnostic residue screen | nothing | at the null in all 18 genes |
| PCOC posteriors | profile-shift model, all 10 lineages converging | 0 of 11,727 sites; maximum posterior 0.098 |
| Contrast-FEL | codon model | 0 of 15,349 sites with different dN/dS |
| RERconverge | rate model | no gene significant; all adjusted p > 0.67 |
| TDG09 | site-specific fitness model | 885 of 3,820 sites flagged, 0 corroborated by any other method |

The two model-free lines carry the claim in full generality, because they
declare no convergent set and therefore assume nothing about how many lineages
must converge together. The PCOC line must be read narrowly (Section 4).

## 2. Evolutionary scenario

Diel activity has significant phylogenetic signal (Pagel's lambda = 0.671,
p = 0.010), so the convergence question is well posed.

The primary equal-rates reconstruction recovers a **nocturnal ancestral placental
mammal without that state being imposed**, consistent with the nocturnal
bottleneck hypothesis, and yields **10 independent gains of diurnality** and
**5 reversals to nocturnality** (15 transition branches: 10 on tip branches,
5 internal). Ten independent origins is comfortably above the threshold at which
convergence methods have usable power.

The all-rates-different reconstruction, reported as a sensitivity analysis,
instead places a **diurnal** ancestor at the base of the placental radiation and
reframes the history as 6 gains and 10 reversals. This is the artefact the
balanced 30/30 taxon sample was expected to produce (Methods Section 1), and it
is why the equal-rates model is primary. Any result quoted under ARD inherits a
premise that contradicts the published mammalian literature.

**Hemiplasy is not driving anything.** Of the 15 transition branches, none sits
on a branch that is weak on both concordance axes (gCF < 50% and sCFL < 50%).
The two closest cases (gCF 56.2 / sCFL 48.1; gCF 50.0 / sCFL 42.0) are weak on
one axis only. Gene-tree discordance therefore cannot plausibly be manufacturing
false convergence at any transition.

## 3. Primary evidence: direct count of convergent substitutions

This is the most informative result in the study, because it counts what
actually happened rather than fitting a model to it, and because it declares no
convergent set and so carries none of PCOC's restrictions.

**A confound had to be removed first.** Raw counts initially appeared to show a
convergent excess (PER2 61 observed against 34.2 null, p = 0.015; PER1 21 against
10.9, p = 0.040). This is an artefact: real transition branches carry far more
substitutions of *any* kind than random branches of matched clade size (NR1D1:
58 opportunity sites against a null of 9.6, a 6-fold excess), because the
reconstruction places transitions on real, often long branches, and more
substitutions mechanically produce more coincidences. Conditioning on
opportunity removes the excess completely:

| gene | raw p (confounded) | rate-conditioned p |
|---|---|---|
| PER2 | 0.015 | 0.418 |
| PER1 | 0.040 | 0.628 |
| CLOCK | 0.085 | 0.600 |
| NPAS2 | 0.129 | 0.302 |

**Pooled across all 18 genes: 648 sites changed in two or more independent
diurnal lineages, and 201 of those changed to the same residue. Observed rate
0.310, against a length-matched null of 0.337.** No excess. The lowest
rate-conditioned p across all 18 genes is 0.282 (CSNK1E); no gene approaches
significance.

**Per-site testing, and the filter it forced.** Twelve of the 201 same-residue
sites reach q <= 0.05. **None survives inspection.** Nine of the twelve carry a
residue as common in nocturnal species as in diurnal ones; the clearest is
CLOCK 675, residue V in two diurnal and two nocturnal species, a
diurnal-minus-nocturnal frequency gap of -0.002, with two of the four carriers
being nocturnal marsupials.

The cause is structural rather than a coding error. The permutation test asks
only whether two or more transition branches changed to the same residue more
often than randomly placed branches would. At a tolerant site a residue arises
repeatedly across the tree, some of those changes land on transition branches by
chance, and because the real transitions are a phylogenetically clustered set
while the null scatters branches more widely, the coincidence looks unusual.
That is homoplasy unrelated to diel activity.

Requiring the convergent residue to be phenotype-specific as well as
statistically unusual:

| filter | sites |
|---|---|
| same residue in 2+ independent diurnal lineages | 201 |
| ... and q <= 0.05 | 12 |
| ... and phenotype-specific (diurnal gap >= 0.25) | **0** |

The largest gap among the twelve is 0.232 (PER1 876), below threshold; the
median is 0.067.

## 4. PCOC, and the limit on how its zero can be read

PCOC returned **zero sites above threshold in all 18 genes**. The highest
posterior observed anywhere in 11,727 columns is **0.0979** (RORB); the next
highest are NPAS2 0.0027 and PER3 0.0010, and every other gene's maximum is
0.0000. This is not a marginal miss.

Power calibration on each gene's own tree and scenario gives **power 1.000 and
false positive rate 0.0000 for every gene at every threshold from 0.70 to 0.99**.
Because power and FPR are saturated across the whole range, the calibration does
not identify a threshold; the strictest value (0.99) was adopted because it costs
no power. This is moot for the present result, since the highest observed
posterior is 0.098.

**Two caveats limit what this zero means, and both are load-bearing.**

**(a) The simulated power is conditional.** `pcoc_sim` builds the convergent
shift with a `OneChange` model that conditions on a substitution occurring on
every transition branch, so the simulated substitution is guaranteed regardless
of branch length. Detection uses the same model, so the calibration is
internally fair, but it establishes only that *given* a convergent substitution
occurred in each lineage, PCOC finds it. It says nothing about whether there was
evolutionary time for those substitutions to occur.

**(b) PCOC only detects near-universal convergence on this design.** The sweep
declares all 10 gains convergent whatever the truth. Simulating convergence in
only *k* of the 10 lineages and detecting with the full declared scenario:

| converging lineages of 10 | power at posterior 0.9 |
|---|---|
| 2 | 0.000 |
| 3 | 0.000 |
| 4 | 0.007 |
| 5 | 0.020 |
| 6 | 0.088 |
| 7 | 0.250 |
| 8 | 0.548 |
| 9 | 0.463 |
| 10 | 1.000 |

A direct check confirms the mechanism: at k = 2, detection under the declared
10-event scenario found **0 of 600** planted sites, while the *same data*
detected under the correct 2-event scenario found **600 of 600**. Perfect
detection against total blindness, on identical data.

The controlling variable is the number of **falsely declared** events, not
branch length. Equal-branch comparisons separate them: k = 2 with 25 branches and
8 false events gives power 0.000, while k = 8 with 37 branches and 2 false events
gives 0.545. PCOC's convergent model requires the derived profile on every
declared branch, so each declared lineage that did not converge actively
penalises the fit. The cost is contaminating evidence, not missing evidence.

**Consequence.** The PCOC result means *no site at which all 10 diurnal lineages
converged together*. That is a much narrower and biologically less likely claim
than *no convergence*. Convergence confined to two or three lineages would be
invisible to PCOC here, which is exactly why the parsimony count of Section 3,
which examined every site where two or more lineages changed, carries the general
negative.

## 5. Model-free residue screen

Removing the evolutionary model entirely and scoring each column by how
diagnostic its best residue is for diurnality, against a null that re-places
diurnal clades of matched size at random on the same gene tree:

**Nothing, in all 18 genes.** The observed maximum score sits on the null
throughout, and in several genes below it (RORB 0.223 against a null of 0.305;
CSNK1D 0.250 against 0.326). Only 2 of 11,727 columns exceed a 0.6 diagnostic
score, against a null expectation of about the same. The lowest p is CRY1 at
0.043, which does not survive correction across 18 genes.

## 6. Selection

**Contrast-FEL: 0 of 15,349 codon sites** show a significant difference in dN/dS
between the diurnal-transition branch set and the rest of the tree, at FDR 0.05
either within gene or across all 18 genes. The smallest uncorrected p anywhere is
2.2e-4, which is unremarkable across 15,349 tests.

**RELAX is not reliably estimable on this dataset and its results should not be
quoted as findings.** Nine of 18 genes converged; five gave up after three
attempts, all failing inside ancestral reconstruction on deeply divergent
lineages (*Lagorchestes hirsutus* nine times, *Phalanger gymnotis* and *Choloepus
hoffmanni* three times each, *Dasypus novemcinctus* once). Its one nominally
significant result, NPAS2 with K = 0.352 at p < 1e-6, **replicated once in eleven
attempts** and is treated as an optimiser artefact on a flat likelihood surface
rather than as evidence of relaxed selection.

## 7. Gene-level rate association

RERconverge finds **no gene with relative evolutionary rates significantly
associated with diel activity**. All 18 adjusted p-values exceed 0.67. The
strongest raw signal is BHLHE40 (rho = -0.178, p = 0.051, adjusted p = 0.677),
which is what one expects as the best of 18 null tests.

Notably, the sign of the correlation is negative in 13 of 18 genes, meaning
diurnal lineages tend to show slightly *slower* relative rates, the opposite of
what accelerated adaptive evolution would predict, though nowhere near
significance.

## 8. Multi-method consensus

TDG09 flagged **885 of 3,820 testable sites** at FDR 0.05, or 23 percent. That
figure is an outlier among the six methods by orders of magnitude, and it does
not survive integration: of the 885 sites, **every one is supported by exactly
one method, and none reaches the two-method threshold for high confidence.**

| n methods supporting | sites |
|---|---|
| 1 | 885 |
| 2 or more | **0** |

Per gene, TDG09 flags between 3 (CSNK1D) and 196 (PER2) sites, while PCOC and
Contrast-FEL flag zero everywhere. **The consensus set is empty.** TDG09's hits
should not be reported without this context; the most parsimonious reading is
that its per-site test is poorly calibrated on alignments of this depth and
divergence, not that it detected 885 real fitness shifts that five other methods
missed.

## 9. Genes that cannot answer the question

Two of the 18 genes are opportunity-starved rather than negative, and their
nulls should not be read as evidence of absence:

- **CSNK1D**: 2.0 expected usable sites; in the observed data, **no position at
  all** has two or more independent diurnal lineages changing. There is nothing
  for a convergence test to examine.
- **ARNTL**: 2.1 expected usable sites; two positions observed.

CSNK1D also has the poorest taxon occupancy (40 of 60 species) and the shortest
mean branch lengths after ARNTL. The remaining 16 genes have ample opportunity
and their nulls are meaningful.

## 10. Validation of the pipeline itself

Ten defects were found over the course of this project. **Not one raised an
error.** Every one produced a confident, plausible, wrong number, and several sat
in code that had already been reviewed. They share a single shape: an unchecked
assumption that two representations of the same object agreed. Two of them
inverted the study's premise outright, by having the pipeline assume a diurnal
ancestral mammal.

The most consequential was a **re-rooting defect**: IQ-TREE emits unrooted trees,
the species tree is rooted on the marsupial outgroup, and descendant-tip-set
matching failed silently on three nodes per gene, one of them a transition node,
which were then dropped from the convergent scenario. After the fix, 175 of 175
scenario groups validate node by node, and the maximum PCOC posterior on real
data fell from 0.3011 to 0.0979. The pre-fix number was the largest apparent
signal in the study and it was an artefact of a dropped transition branch.

`scripts/lib_checks.py` now asserts these interface conditions at every handoff,
raising rather than warning. `scripts/test_lib_checks.py` contains 15 cases, each
reconstructing an actual historical defect from this project and asserting both
that the check fires on the broken input and that it passes on the corrected
input. All 15 pass.

**An end-to-end positive control is reported in Section 11.** Component-level
validation does not establish that an assembled pipeline carries signal, and
neither does `pcoc_sim`, which validates PCOC against PCOC's own model using its
own scenario handling and therefore cannot detect a scenario that was built
incorrectly.

## 11. End-to-end positive control

STATUS: RUNNING. This section is the gating result for the strength, though not
the direction, of everything above, and will be completed when the control
finishes.

Ninety convergent sites were planted into the real alignments across all 18
genes, at three difficulty levels, and the entire pipeline was re-run unchanged
over the spiked data using the same trees, the same scenario files and the same
scripts.

| level | converging events | sites planted | mean diurnal gap | mean species changed |
|---|---|---|---|---|
| `all` | all 10 | 36 | 1.000 | 28.6 |
| `most` | 7 of 10 | 18 | 0.816 | 23.4 |
| `few` | 3 of 10 | 36 | 0.374 | 10.7 |

Zero nocturnal species carry a planted residue at any level.

**A first version of this control was faulty and every conclusion drawn from it
was withdrawn.** It planted into all species descending from each event rather
than the diurnal descendants, which is not equivalent: an event's membership
includes descendant nodes whose subtrees contain nocturnal species from nested
reversals. That planted into 41 of 60 species rather than the intended 30, making
the residue the majority state across the tree in 36 of 36 `all` sites. A new
consensus residue is not convergence, so the detectors were being asked to find
something that was not there, and PCOC was correct to miss it.

## 12. What this study does and does not establish

**Establishes.** Across 18 core circadian genes, 60 mammals and 10 independent
origins of diurnality, there is no convergent amino acid signal detectable by six
methods spanning model-free counting, profile-shift modelling, site-specific
fitness modelling, codon-level selection and relative-rate association. The
strongest single statement is the parsimony count: of 648 sites where two or more
independent diurnal lineages changed, the same-residue rate is at chance, and not
one statistically unusual site is also phenotype-specific.

**Does not establish.** That diel activity evolves without molecular change. The
leading alternative, that adaptation to a diurnal niche acts through *regulatory*
change (promoter and enhancer evolution, expression level, phase relationships)
rather than through protein-coding substitution, is untested here and is the
natural next experiment. Nor does it establish anything about genes outside this
18-gene panel, in particular light-input pathways (opsins, retinal signalling)
and clock-output effectors, where the phenotype-relevant variation may sit.

**Cannot answer.** CSNK1D and ARNTL, for lack of evolutionary opportunity
(Section 9).

**Conditions on.** An artificially balanced 30/30 diurnal-nocturnal sample, which
inflates the apparent rate of diurnality relative to real mammals and is what
drives the ARD reconstruction to an implausible diurnal ancestor. The equal-rates
reconstruction is robust to this, but the sampling should be stated in any
report of these results.
