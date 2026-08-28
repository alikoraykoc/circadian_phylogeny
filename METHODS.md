# Methods

Convergent molecular evolution in mammalian circadian genes and independent
transitions in diel activity.

All analysis code, configuration and result tables are in the project
repository. Numbered scripts (`scripts/00` through `scripts/13`) run in the
order given in `README.md`; helper scripts named in each subsection below are
in the same directory.

---

## 1. Taxon sampling and phenotype coding

Sixty mammalian species were sampled, spanning all major placental orders plus a
marsupial outgroup. Diel activity was coded as a binary character from the
published compilations of Maor et al. (2017) and Bennie et al. (2014):
30 species scored diurnal and 30 scored nocturnal (`data/diel_activity.csv`).
No species in the final set carried an intermediate (crepuscular or cathemeral)
score, so no intermediate states were collapsed.

Species were selected to maximise the number of independent transitions between
activity states while retaining phylogenetic spread. **This produces a sample
that is balanced 30/30, which real mammals are not**, mammals being
predominantly nocturnal. The consequence is that a reconstruction fitted to
these tips alone is pulled toward a diurnal ancestor. This is the reason the
ancestral state reconstruction is reported under two models rather than assumed
(Section 5), and it is stated here because it conditions every downstream
result.

## 2. Species tree

The topology was pruned from the mammalian supertree of Upham et al. (2019).
This tree is the **fixed scaffold for every analysis in the study**: the diel
scenario is mapped onto it once, so that transition branches sit at identical
nodes across all 18 genes and per-gene results are directly comparable.

Three distinct tree sets are used, and conflating them is the single most
consequential thing a reader could get wrong:

1. **The Upham species tree**, fixed topology, used as the scaffold.
2. **Per-gene fixed-topology trees**, the Upham topology carrying gene-specific
   branch lengths in substitutions per site (Section 4). These are what PCOC,
   TDG09, RERconverge and HyPhy consume.
3. **Unconstrained per-gene ML trees and a supermatrix tree**, used as support
   tools only: orthology quality control, gene concordance factors for the
   hemiplasy screen, topology confirmation, and a uniform-model fallback branch
   length set. These are never used as an analysis backbone.

Forcing the species topology removes topology-error artefacts; the concordance
step (Section 6) then guards against the remaining hemiplasy artefact.

## 3. Sequences, alignment and trimming

Per-gene protein sequences were collected for 18 core circadian genes covering
all functional modules of the transcription-translation feedback loop:

| Module | Genes |
|---|---|
| Positive arm | CLOCK, NPAS2, ARNTL |
| Negative arm | PER1, PER2, PER3, CRY1, CRY2 |
| Auxiliary loop | NR1D1, NR1D2, RORA, RORB, RORC |
| Post-translational regulators | CSNK1D, CSNK1E, FBXL3 |
| Output and light response | BHLHE40, BHLHE41 |

Taxon occupancy is incomplete for some genes; per-gene alignments contain
between 40 (CSNK1D) and 60 species. Alignment tip labels were verified to match
species tree tips exactly before any tree was built (`scripts/00_check_inputs.sh`);
one label mismatch (`Gorilla_gorilla_gorilla` against `Gorilla_gorilla`) was
corrected across all 18 files. Label mismatch is the commonest cause of silent
failure in fixed-topology analyses and is checked rather than assumed.

Alignments were trimmed with **trimAl v1.5.rev1** using the `-automated1`
heuristic (`scripts/01_trim_alignments.sh`). For each gene, `-colnumbering`
recorded, for every retained column, its index in the untrimmed alignment. That
mapping, combined with the *Homo sapiens* alignment row, is what converts a
trimmed column index back to a human residue number when sites are reported.

| gene | taxa | untrimmed cols | trimmed cols |
|---|---|---|---|
| CLOCK | 59 | 991 | 841 |
| NPAS2 | 58 | 1145 | 814 |
| ARNTL | 60 | 778 | 620 |
| PER1 | 60 | 1643 | 1233 |
| PER2 | 60 | 1726 | 1208 |
| PER3 | 56 | 2122 | 897 |
| CRY1 | 59 | 670 | 582 |
| CRY2 | 60 | 665 | 592 |
| NR1D1 | 60 | 760 | 610 |
| NR1D2 | 60 | 660 | 562 |
| RORA | 58 | 602 | 466 |
| RORB | 60 | 665 | 459 |
| RORC | 60 | 790 | 518 |
| CSNK1D | 40 | 792 | 401 |
| CSNK1E | 52 | 564 | 509 |
| FBXL3 | 60 | 442 | 426 |
| BHLHE40 | 60 | 666 | 407 |
| BHLHE41 | 58 | 768 | 582 |

Total analysed protein alignment: 11,727 trimmed columns.

Coding sequences for the selection track were retrieved separately and matched
to each protein by **transcript accession, not by gene name**, since isoform
mismatch is the dominant failure mode in codon back-translation. Sequences came
from RefSeq where a curated transcript existed and from `miniprot` alignment to
the genome otherwise; every CDS was verified to translate to its protein
alignment row, and sequences failing that check (length not a multiple of three,
internal stop codons, translation mismatch) were dropped and recorded
(`codon_export/cds_qc_report_ALL.csv`). Codon alignments were produced by
back-translating the protein alignment with **pal2nal**.

## 4. Branch lengths on the fixed topology

Branch lengths were estimated with **IQ-TREE 3.1.1** under the fixed Upham
topology (`-te`), which optimises lengths without searching topology
(`scripts/05_branchlengths_fixed.sh`). Where a gene lacked species present in
the species tree, the tree was pruned to that gene's taxa first.

Two length sets were produced, for two different reasons:

- **Per-gene best-fit models (`-m MFP`)** for PCOC, TDG09 and the selection
  track. These analyse each gene independently with no cross-gene comparison, so
  each gene should get its own best-fit branch lengths.
- **A single uniform model across all genes** for RERconverge. The model was
  selected by running ModelFinder once on the concatenation of all 18 trimmed
  alignments (Q.MAMMAL+F+R6). RERconverge normalises overall rate scale but not
  model-driven differences in branch length *shape*: if one gene were fitted with
  a four-category rate model and another with six, the resulting long-to-short
  branch ratios would differ for reasons unrelated to biology, and that
  difference would enter the cross-gene rate comparison as a confound.

Mean branch length spans a 37-fold range across genes, from 0.0035 subs/site in
ARNTL to 0.130 in CSNK1E. This range is why power is calibrated per gene on each
gene's own tree (Section 8) rather than once for the dataset.

**Re-rooting.** IQ-TREE emits unrooted trees. The species tree is rooted on the
marsupial outgroup, and the two rootings induce different clade decompositions
at nodes near the root, so descendant-tip-set matching between a gene tree and
the species tree fails silently on those nodes. Every per-gene tree is therefore
re-rooted on the marsupial outgroup before scenario construction
(`scripts/reroot_gene_trees.py`), with the unrooted original retained. This is
not cosmetic: before it was introduced, three nodes per gene, including one
transition node, failed to match and were dropped from the convergent scenario
without any error being raised.

## 5. Ancestral state reconstruction and definition of convergent events

Ancestral diel states were reconstructed with **corHMM 2.8** on the fixed
species tree, using marginal reconstruction with one rate category and a free
root (`scripts/07_scenario.R`). States are coded 1 = nocturnal (ancestral),
2 = diurnal (derived). Phylogenetic signal in the trait was measured with
Pagel's lambda in **phytools 2.5.2** (lambda = 0.671, p = 0.010), confirming
that activity state is phylogenetically structured rather than randomly
distributed, which is a prerequisite for the convergence question to be
well posed.

Two models were fitted:

- **ER (equal rates), the primary analysis.** Best AIC. Recovers a **nocturnal**
  ancestral placental mammal without that state being imposed, consistent with
  the nocturnal bottleneck hypothesis. Yields **10 independent gains of
  diurnality** and **5 reversals to nocturnality**.
- **ARD (all rates different), reported as a sensitivity analysis.** Places a
  **diurnal** ancestor at the base of the placental radiation and reframes the
  history as 6 gains and 10 reversals. This is the artefact the balanced 30/30
  sample was expected to produce (Section 1), and it is the reason ER is primary
  rather than the reverse.

**Direction is never merged.** PCOC fits a single convergent amino acid profile
per run. Placing branches that move toward opposite phenotypes into one
convergent class is incoherent and self-cancelling, so gains of diurnality and
reversals to nocturnality are analysed as separate convergent classes throughout.

**Convergent event membership.** For each transition onto a child node, the
convergent event comprises that node plus every descendant branch that remains
in the derived state. Descent stops at any node that reverts, so a nested
reversal is excluded from the convergent group rather than swept into it. The
transition node is always the first member, which is what PCOC's scenario syntax
requires.

**Node identity is keyed by descendant tip set, not by node number.** `ape`,
`ete3`/IQ-TREE and PCOC each number internal nodes differently, and a node
number is meaningless once a tree is pruned to a gene's taxa. `07_scenario.R`
therefore publishes `node_tipsets_ER.csv`, giving every node's reconstructed
state alongside its sorted descendant tip set, and every downstream consumer
(PCOC scenario construction, TDG09 node labelling, HyPhy foreground labelling,
the hemiplasy screen) joins on that key rather than re-deriving states.

## 6. Hemiplasy control

Incomplete lineage sorting can make ancestral polymorphism look like convergence
when a gene tree differs from the species tree. Gene and site concordance
factors were computed in IQ-TREE (`scripts/06_concordance.sh`) from the
unconstrained per-gene trees against the species tree, and each diel transition
branch was matched to its concordance values by descendant tip set
(`scripts/08_flag_discordant_branches.py`).

A branch is flagged only when **both** gCF < 50% and sCFL < 50%, since either
alone is usually estimation noise; missing values never flag. Of the 15 ER
transitions (10 on tip branches, 5 internal), **none was flagged**. The two
closest cases (gCF 56.2 / sCFL 48.1, and gCF 50.0 / sCFL 42.0) are weak on one
axis only.

## 7. Convergence detection: PCOC

PCOC (Rey et al. 2018) was run with `pcoc_det.py` from the
`carinerey/pcoc` Docker image (digest `sha256:11ea18fb9b96...`), taking each
gene's re-rooted fixed-topology tree, its trimmed protein alignment, and the
ER gain-of-diurnality scenario (`scripts/run_pcoc_all.sh`). Gamma-distributed
rate heterogeneity was enabled (`--gamma`) and the posterior reporting threshold
was set to 0 (`-f 0.0`) so that the full posterior distribution over sites was
retained rather than only sites passing a cutoff.

Scenario strings were built by `scripts/prep_pcoc_scenarios.py` from the
convergent event table, with the transition node placed first in each
`/`-separated group as PCOC requires, and validated node by node against the
gene tree (175 of 175 groups correct after the re-rooting fix of Section 4).

Four scenario sets were built, one per combination of reconstruction model and
transition direction, and each is calibrated and swept independently because
event counts differ enough between them to change power (counts shown for
CLOCK):

| set | events | convergent branches | role |
|---|---|---|---|
| `ER_gain` | 10 | 55 | primary analysis |
| `ER_reversal` | 5 | 19 | reversals to nocturnality, the second convergent class |
| `ARD_gain` | 6 | 64 | ASR model sensitivity |
| `ARD_reversal` | 10 | 44 | built but not run (see below) |

Gains and reversals are never merged into one set, for the reason given in
Section 5. The two ARD sets are sensitivity analyses on the reconstruction model
rather than independent hypothesis tests: ARD places a diurnal ancestral
placental mammal, so a null under ARD adds robustness, while a positive under
ARD would have to be reported as conditional on a reconstruction this study
otherwise argues against.

Only `ARD_gain` was run. The gain direction is where the reconstruction model
actually changes the question, since ER and ARD disagree about whether diurnality
is derived at all; testing the reversal direction under a model that already
places a diurnal ancestor examines a scenario the study argues against twice
over, at a cost of nine hours of computation. `ARD_reversal` scenarios were built
and retained so the analysis can be run if a reviewer asks.

## 8. PCOC power calibration

The project's design rules forbid a fixed posterior threshold. Each gene was
therefore calibrated on **its own tree and its own scenario** with `pcoc_sim.py`
(`scripts/run_pcoc_sim.sh`, summarised by `scripts/summarize_pcoc_sim.py`):
100 simulated convergent sites and 100 null sites per profile couple, 10 profile
couples per gene, all 18 genes.

**One implementation detail materially affects the result.** `pcoc_sim` treats
`-m` only as a *pool* of candidate events and simulates a random subset whenever
`-c` is left at its default of 0 (`pcoc_sim.py:496`; deletion happens in
`events_placing.py:placeNTransitionsInTree`). A pilot run silently calibrated on
7 of the 10 real events. The production script passes `-c` equal to each gene's
real event count, verified by diffing the simulated annotated tree against the
scenario node by node.

Power was 1.000 and false positive rate 0.0000 for every gene at every threshold
tested (0.70 to 0.99). Because power and FPR are saturated across the whole
range, **the calibration does not identify a threshold**; the strictest value
(0.99) was adopted on the grounds that it costs no power. This is reported
explicitly rather than presented as an optimised choice.

**A limitation that must be read with the power figure.** `pcoc_sim` builds the
convergent shift with the Bio++ `OneChange` model on every transition branch
(`bpp_lib.py`: `MODEL_OC='OneChange(model=$(MODEL_C))'`), which *conditions on*
at least one substitution occurring there. The simulated substitution is
guaranteed regardless of branch length. Since `pcoc_det` fits the same model,
simulation and detection are matched and the calibration is internally fair, but
the reported power is **conditional**: given that a convergent substitution
occurred in each lineage, PCOC finds it. It says nothing about whether there was
evolutionary time for those substitutions to occur. Sections 10 and 11 address
that separately.

## 9. PCOC sensitivity to partial convergence

The step-9 sweep declares all 10 gains of diurnality convergent, whatever the
truth. To measure the cost of that assumption, `scripts/pcoc_partial_convergence_sweep.py`
simulates convergence in only *k* of the 10 lineages and then detects with the
full declared 10-event scenario, exactly as the real analysis does, for
k = 2..10 with independent draws of *which* lineages at each k.

Power at posterior 0.9 rises from 0.000 at k = 2 or 3, through 0.020 at k = 5
and 0.250 at k = 7, to 1.000 only at k = 10. A direct check confirms the
mechanism: with k = 2, detection under the declared 10-event scenario found 0 of
600 planted sites, while the *same data* detected under the correct 2-event
scenario found 600 of 600.

The controlling variable is the number of **falsely declared** events, not
branch length; equal-branch comparisons separate the two (k = 2 with 25 branches
and 8 false events gives power 0.000, while k = 8 with 37 branches and 2 false
events gives 0.545). PCOC's convergent model requires the derived profile on
every declared branch, so each declared lineage that did not converge actively
penalises the fit. The cost is contaminating evidence, not missing evidence.

**Consequence.** A PCOC zero on this design means *no site at which all 10
diurnal lineages converged together*. It does not mean *no convergence*. This is
why the model-free tests of Sections 10 and 11, which declare no convergent set
at all, carry the general claim.

## 10. Model-free diagnostic residue screen

`scripts/model_free_convergence_screen.py` removes the evolutionary model
entirely. For each alignment column the statistic is

```
score = max over residues a of ( freq(a | diurnal) - freq(a | nocturnal) )
```

which is 1.0 for a perfectly phenotype-diagnostic column and near 0 for a column
whose residues ignore the phenotype. Columns with fewer than 5 ungapped residues
in either group are skipped.

**The null is not a label shuffle.** Diurnality arose in a small number of
clustered clades, so shuffling tip labels destroys the phylogenetic clumping and
would make almost any clade-restricted residue appear significant. Instead the
null re-places clades *of the observed sizes* at random positions on the same
gene tree (1000 permutations per gene), so the null carries the same
autocorrelation structure and only the phenotype assignment is randomised.

## 11. Direct parsimony count of convergent substitutions

`scripts/observed_convergent_substitutions.py` is the study's primary evidence,
because it counts what actually happened rather than fitting a model to it, and
because it declares no convergent set and so carries none of the restriction
described in Section 9.

Fitch parsimony assigns an amino acid to every internal node per column, so a
substitution on a branch is simply parent state not equal to child state. Five
features of the implementation are load-bearing:

1. **It conditions on opportunity.** Raw counts are confounded. Real transition
   branches carry far more substitutions of *any* kind than random branches of
   matched clade size (NR1D1: 58 opportunity sites against a null of 9.6),
   because the reconstruction places transitions on real, often long branches,
   and more substitutions mechanically produce more same-residue coincidences.
   The test statistic is therefore the *fraction* of sites changing in two or
   more independent lineages that landed on the same residue, not the count.

2. **The null is matched on branch length as well as clade size.** Null branches
   must fall within a factor of 2.0 of the real branch length. Conditioning on
   opportunity corrects the symptom; length-matched sampling removes the cause.
   Both are reported and they agree.

3. **Ancestral-state uncertainty is propagated.** Fitch leaves ties, and
   resolving them deterministically silently commits to one history. Each site
   is resolved 10 times with random tie-breaking.

4. **Individual sites are tested, not only the aggregate.** The pooled rate is
   well powered against a pervasive signal but blind to a handful of real sites,
   which is the more likely biology. Every site receives a permutation p-value
   (2000 permutations) and a Benjamini-Hochberg q-value across all genes at
   FDR 0.05.

5. **The convergent residue must be phenotype-specific.** Points 1 to 4 are not
   sufficient. Each candidate site also carries `diurnal_gap`, the frequency of
   the convergent residue among diurnal species minus its frequency among
   nocturnal species, and a site is reported as convergent only if it is both
   statistically unusual and phenotype-specific (gap >= 0.25). The necessity of
   this filter is empirical: without it, 12 sites reached q <= 0.05, and 9 of
   them carried a residue as common in nocturnal species as in diurnal ones
   (CLOCK 675, residue V in two diurnal and two nocturnal species, gap -0.002).
   The permutation test alone asks only whether two or more transition branches
   changed to the same residue more often than randomly placed branches would;
   at a tolerant site a residue arises repeatedly across the tree, some of those
   changes land on transition branches by chance, and because the real
   transitions are a phylogenetically clustered set while the null scatters
   branches more widely, the coincidence looks unusual. That is homoplasy
   unrelated to diel activity, not convergence.

For efficiency, the parsimony reconstruction does not depend on which branches
are called transitions, so it is computed once per site and each permutation
becomes a lookup; this is what makes 2000 permutations affordable.

## 12. Site-specific fitness shifts: TDG09

TDG09 v1.1.2 (Java 24) tests, per site, whether amino acid fitness differs
between two groups of lineages (`scripts/10_tdg09.sh`). Inputs were built by
`scripts/prep_tdg09_inputs.py`, which labels every node of the per-gene tree
with its reconstructed state **read from `node_tipsets_ER.csv`** rather than
derived by propagating from an assumed root, and which raises rather than warns
if any node requires a majority-rule fallback label.

Two details determine whether the analysis answers the intended question:

- **Group argument order sets the root group.** `-groups No Di` declares the
  nocturnal group ancestral. The reverse order makes TDG09 assume a diurnal
  ancestral mammal, inverting the study's premise. This is not documented
  prominently in the tool and was initially wrong.
- **Results must be parsed from the `FullResults:` block**, not from the
  `# Result{...}` progress log, which reports `lrt=0.0` at every site regardless
  of the real value. A parser reading the progress log returns a confident,
  entirely wrong table.

Per-site p-values were corrected within each gene by Benjamini-Hochberg.

## 13. Gene-level rate association: RERconverge

**RERconverge 0.3.0** (`scripts/11_rerconverge.R`) tests which genes have
relative evolutionary rates associated with diel activity. It consumes the
uniform-model fixed-topology trees of Section 4. Relative rates were computed
with `getAllResiduals(transform = "sqrt")`. The foreground was defined as
diurnal tips with `clade = "terminal"`, the conservative choice, which uses tip
branches only rather than also including ancestral branches leading into diurnal
clades. Association was tested with `correlateWithBinaryPhenotype`
(`min.sp = 10`, `min.pos = 2`), and p-values corrected across the 18 genes.

## 14. Selection: Contrast-FEL and RELAX

Codon-level tests were run in **HyPhy 2.5.93** (`scripts/12_selection_hyphy.sh`)
on the codon alignments of Section 3.

Foreground labelling and tree preparation are handled by
`scripts/prep_selection_inputs.py`, which does three things that matter:

- **Branch lengths come from the per-gene substitution trees, not from the Upham
  time tree.** The Upham tree is scaled in millions of years, and feeding branch
  lengths near 46 into HyPhy's GTR fitting crashes it.
- **All-gap codon columns are stripped**, which is required because HyPhy raises
  an assertion failure on them; 13 of 18 genes contained such columns (PER3 482,
  CSNK1D 239).
- **A per-gene coordinate map (`<gene>.codonmap.csv`) is written**, so that
  post-strip codon indices can be mapped back to untrimmed protein columns and
  then to human residues. Without it, reported site numbers would be silently
  shifted.

Contrast-FEL was run with `--branch-set Foreground` to test, per site, whether
dN/dS differs between the diurnal-transition branch set and the rest of the
tree. p-values were corrected both within gene and across all genes.

**RELAX is not reliably estimable on this dataset and its results are reported
as such, and no attempt was made to force convergence by removing taxa.** Nine of 18 genes converged; five gave up after three attempts, all
failing inside ancestral reconstruction on deeply divergent lineages
(*Lagorchestes hirsutus* nine times, *Phalanger gymnotis* and *Choloepus
hoffmanni* three times each, *Dasypus novemcinctus* once). Its one nominally
significant result (NPAS2, K = 0.352, p < 0.0001) replicated once in eleven
attempts and is treated as an optimiser artefact rather than a finding.

## 15. Multi-method consensus

`scripts/13_consensus.py` integrates the site-level methods. It parses PCOC
posteriors, TDG09 `FullResults` LRTs and Contrast-FEL JSON, maps every site
through the trimmed-column to untrimmed-column to human-residue chain using the
trimAl `colnumbering` files and the *Homo sapiens* alignment row, and requires
support from at least two methods for a site to be called high confidence.

Two checks are enforced at parse time: every mapped coordinate must resolve
within range, and each parser must return as many rows as its input claims. The
row check validates rows *seen*, not rows *retained*, which is the distinction
that catches a parser silently reading the wrong block of a results file.

## 16. End-to-end positive control

A null result is not a finding until the pipeline is shown capable of detecting
a signal. Component-level validation does not establish this, and neither does
`pcoc_sim`, which validates PCOC against PCOC's own model using its own scenario
handling and therefore cannot detect a scenario that was built incorrectly.

`scripts/spikein_generate.py` plants known convergent sites into the **real**
alignments at recorded positions; the entire pipeline then runs unchanged over
the spiked data, using the same trees, the same scenario files and the same
scripts (`scripts/spikein_run.sh`). Only the alignment differs. Recovery is
scored against a sealed answer key by `scripts/spikein_evaluate.py`, which is
the only script permitted to read it.

Sites are planted by writing a residue that is **absent at that column** into
the **diurnal species descending from a chosen set of gain events**. Candidate
columns must be variable but not saturated (2 to 6 residues present, at most 10%
gaps). Three difficulty levels are planted, chosen to span the sensitivity range
identified in Section 9:

| level | converging events | purpose |
|---|---|---|
| `all` | all 10 | PCOC should recover these |
| `most` | 7 of 10 | borderline for PCOC |
| `few` | 3 of 10 | PCOC should miss; model-free methods should not |

Every planted site is verified before it is kept: the residue must be
essentially absent from nocturnal species (frequency <= 0.02), the diurnal minus
nocturnal frequency gap must exceed 0.10, and at least two independent events
must contribute. Sites failing any criterion are reverted and another column is
tried.

**A first version of this control was faulty and its conclusions were
withdrawn.** It wrote the residue into every species descending from each event,
which is not equivalent to writing it into the diurnal descendants: an event's
membership includes descendant nodes whose subtrees contain nocturnal species
from nested reversals. That planted into 41 of 60 species rather than the
intended 30, making the planted residue the majority state across the tree in
36 of 36 sites at the `all` level. A new consensus residue is not convergence,
so the control was asking the detectors to find something that was not there.
The corrected version plants into diurnal descendants only, verified: mean
diurnal-nocturnal gap 1.000 (`all`), 0.816 (`most`) and 0.374 (`few`), with
zero nocturnal carriers at every level. The design spread turns the control
into a measurement of which method recovers which difficulty level, and it
measures TDG09's false positive rate directly.

## 17. Pipeline validation

Ten defects were found over the course of this project. **Not one raised an
error.** Every one produced a confident, plausible, wrong number, and several
sat in code that had already been reviewed. They share a single shape: an
unchecked assumption that two representations of the same object agreed, namely
`ape` node numbering against `ete3` numbering, species-tree clades against
gene-tree clades under a different rooting, node labels in memory against node
labels on disk, a command-line flag's argument order against a tree's root
state, alignment columns before pruning against after, and a tool's progress log
against its results table. Two of them inverted the study's premise outright, by
having the pipeline assume a diurnal ancestral mammal.

`scripts/lib_checks.py` therefore asserts these interface conditions at every
handoff, raising rather than warning, since a warning is what allowed the
re-rooting defect to persist:

- `check_same_taxa` alignment headers exactly equal tree tips;
- `check_same_rooting` two trees agree on the outgroup clade;
- `check_scenario_against_tree` every scenario node exists, each group's first
  entry is its transition node, and total branch count matches the event table;
- `check_coordinate_map` a coordinate map is total over its domain and lands in
  range;
- `check_parsed_rows` a parser returned as many rows as its input claims;
- `check_no_fallbacks` no node required a majority-rule fallback label.

`scripts/test_lib_checks.py` contains 15 cases, each reconstructing an actual
historical defect from this project (using the retained pre-fix artefacts where
available, for example `ARNTL.treefile.unrooted`) and asserting both that the
corresponding check fires on the broken input and that it passes on the
corrected input. A check that has never failed is not known to work.

## 18. Software and reproducibility

| tool | version | use |
|---|---|---|
| IQ-TREE | 3.1.1 | branch lengths, concordance factors, model selection |
| trimAl | v1.5.rev1 | alignment trimming, column numbering |
| corHMM | 2.8 | ancestral state reconstruction |
| phytools | 2.5.2 | Pagel's lambda |
| ape | 5.8.1 | tree manipulation in R |
| RERconverge | 0.3.0 | gene-level rate association |
| PCOC | `carinerey/pcoc` Docker, `sha256:11ea18fb9b96...` | convergence detection and simulation |
| TDG09 | 1.1.2 (Java 24.0.2) | site-specific fitness shifts |
| HyPhy | 2.5.93 | Contrast-FEL, RELAX |
| pal2nal | - | codon back-translation |
| miniprot | - | CDS retrieval where RefSeq lacked a transcript |
| R | 4.4.1 | corHMM, RERconverge |
| Python | 3.11.3 (ete3 3.1.3, numpy 1.25.2, plotly 5.18.0) | scenario construction, model-free tests, consensus |

Figures were produced with Plotly and exported as vector PDF. All scripts are
resumable and skip completed genes, so an interrupted run can be restarted
without recomputation. Random seeds are fixed where stochastic (`SEED = 20260814`
for spike-in site selection; `set.seed(1)` for corHMM fitting).

## 19. Statistical reporting

Multiple testing was controlled with Benjamini-Hochberg FDR at 0.05, applied
within gene for TDG09 and Contrast-FEL and across all 18 genes for the parsimony
site test and RERconverge. Permutation nulls preserve phylogenetic clumping by
re-placing clades of the observed sizes rather than shuffling tip labels, in both
the model-free screen (Section 10) and the parsimony test (Section 11).
Where a test is uninformative rather than negative, this is stated: ARNTL and
CSNK1D have too little evolutionary opportunity to answer the question
(2.1 and 2.0 expected usable sites respectively; CSNK1D has no position at which
two or more independent lineages changed at all), and their nulls should not be
read as evidence of absence.
