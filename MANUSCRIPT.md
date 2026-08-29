# Convergent molecular evolution in circadian genes is not associated with independent transitions to diurnality in mammals

Manuscript Methods and Results. Internal technical documentation, including
pipeline validation history, is in `METHODS.md`, `RESULTS.md` and
`PROGRESS_REPORT.md`; this file contains only the manuscript sections.

---

# Methods

## Taxon sampling and phenotype data

Sixty mammalian species were sampled across the placental orders together with a
marsupial outgroup. Diel activity was coded as a binary character from published
compilations (Maor et al. 2017; Bennie et al. 2014), giving 30 diurnal and 30
nocturnal species. No sampled species carried an intermediate crepuscular or
cathemeral score, so no intermediate states were collapsed. Taxa were chosen to
maximise the number of independent transitions between activity states while
retaining phylogenetic spread. This yields a diel sample that is deliberately
balanced and therefore does not reflect the nocturnal predominance of mammals as
a whole; the consequences are addressed in the ancestral state reconstruction
below.

## Sequence data, alignment and trimming

Protein sequences were collected for 18 core circadian genes spanning the
functional modules of the transcription-translation feedback loop: the positive
arm (*CLOCK*, *NPAS2*, *ARNTL*), the negative arm (*PER1*, *PER2*, *PER3*,
*CRY1*, *CRY2*), the auxiliary loop (*NR1D1*, *NR1D2*, *RORA*, *RORB*, *RORC*),
post-translational regulators (*CSNK1D*, *CSNK1E*, *FBXL3*) and the output
repressors (*BHLHE40*, *BHLHE41*). Taxon occupancy is incomplete for some genes;
per-gene alignments contain between 40 and 60 species (Table 1).

Alignment headers were verified against species tree tip labels before any
downstream analysis. Alignments were trimmed with trimAl v1.5.rev1 under the
`-automated1` heuristic, retaining for each gene a map from every trimmed column
to its index in the untrimmed alignment. That map, combined with the *Homo
sapiens* alignment row, converts trimmed column indices to human residue
numbering for reporting. The trimmed dataset comprises 11,727 columns.

Coding sequences for the selection analyses were retrieved by matching transcript
accession rather than gene name, since isoform mismatch is the principal failure
mode in codon back-translation. Sequences came from RefSeq where a curated
transcript existed and from `miniprot` genome alignment otherwise. Every coding
sequence was verified to translate to its corresponding protein alignment row,
and sequences failing that check were excluded. Codon alignments were generated
by back-translating protein alignments with pal2nal.

## Phylogenetic framework

The topology was pruned from the mammalian supertree of Upham et al. (2019) and
held fixed for every analysis, so that transition branches occupy identical nodes
across all 18 genes. Branch lengths were estimated in IQ-TREE v3.1.1 under this
fixed topology (`-te`), which optimises lengths without topology search. Where a
gene lacked taxa present in the species tree, the tree was pruned to that gene's
taxa before estimation.

Two branch length sets were produced. Analyses treating genes independently
(PCOC, TDG09, the selection tests) used per-gene best-fit models selected by
ModelFinder. RERconverge, which compares rates across genes, used a single
uniform model (Q.MAMMAL+F+R6) selected by running ModelFinder once on the
concatenation of all 18 trimmed alignments; differing rate-heterogeneity models
across genes would otherwise introduce systematic differences in branch length
shape unrelated to biology. Mean branch length spans a 37-fold range across genes
(0.0035 substitutions per site in *ARNTL* to 0.130 in *CSNK1E*), which is why
detection power was calibrated per gene rather than once for the dataset.

Because IQ-TREE returns unrooted trees while the species tree is rooted on the
marsupial outgroup, and the two rootings induce different clade decompositions
near the root, every per-gene tree was re-rooted on the outgroup before scenario
construction.

Unconstrained per-gene maximum likelihood trees and a partitioned supermatrix
tree were estimated separately and used only for orthology quality control,
concordance factor estimation and topology confirmation. They were not used as an
analysis backbone.

## Ancestral state reconstruction and definition of convergent events

Ancestral diel states were reconstructed with corHMM v2.8 using marginal
reconstruction, one rate category and a free root. Phylogenetic signal in the
trait was assessed with Pagel's lambda in phytools v2.5.2.

Two models were fitted. The equal-rates (ER) model was the primary analysis, and
the all-rates-different (ARD) model was retained as a sensitivity analysis on the
reconstruction. Because gains of diurnality and reversals to nocturnality move
toward opposite phenotypes, they were analysed as separate convergent classes
throughout; PCOC fits a single convergent amino acid profile per run, so merging
them would be self-cancelling.

For each transition, the convergent event comprises the transition branch plus
every descendant branch remaining in the derived state. Descent terminates at any
node that reverts, so nested reversals are excluded from the convergent group.
Node identity was keyed by descendant tip set rather than node index, since node
numbering differs between `ape`, `ete3`/IQ-TREE and PCOC and is not preserved
under pruning. A node state table keyed by tip set was published once and read by
every downstream consumer.

## Hemiplasy control

Gene and site concordance factors were computed in IQ-TREE from the unconstrained
per-gene trees against the species tree, and each diel transition branch was
matched to its concordance values by descendant tip set. A branch was flagged
when both gCF and sCFL fell below 50 percent; either alone is commonly estimation
noise, and missing values never flagged.

## Convergence detection

PCOC (Rey et al. 2018) was run via `pcoc_det.py` from the `carinerey/pcoc`
container, taking each gene's re-rooted fixed-topology tree, its trimmed protein
alignment and the relevant convergent scenario. Gamma-distributed rate
heterogeneity was enabled and the posterior reporting threshold was set to zero
so that the full posterior distribution over sites was retained.

The posterior threshold was calibrated per gene on each gene's own tree and
scenario using `pcoc_sim.py`, with 100 simulated convergent sites and 100 null
sites per profile couple and 10 profile couples per gene. The number of simulated
events was set explicitly to each gene's real event count, since the simulator
otherwise samples a random subset of the supplied events. Calibration was
performed independently for each of the three scenario sets analysed (Table 2).

Two properties of the calibration bound its interpretation. First, `pcoc_sim`
constructs the convergent shift with a `OneChange` model that conditions on a
substitution occurring on every transition branch, so simulated power is
conditional on a substitution having occurred rather than reflecting whether
there was evolutionary time for one. Second, the detection step declares all
events convergent regardless of the truth. To quantify the cost of that
assumption, convergence was simulated in only k of the 10 gain lineages and
detected with the full declared scenario, for k = 2 to 10 with three independent
draws of which lineages at each k.

## Site-specific fitness shifts

TDG09 v1.1.2 was used to test, per site, whether amino acid fitness differs
between nocturnal and diurnal lineages. Tree nodes were labelled from the
published node state table rather than by propagating states from an assumed
root, and the group argument order was set so that the nocturnal group is
ancestral. Results were parsed from the `FullResults` block; per-site p-values
were corrected within gene by the Benjamini-Hochberg procedure.

## Model-free residue screen

For each alignment column the statistic was the maximum over residues of
freq(residue | diurnal) minus freq(residue | nocturnal), which equals 1.0 for a
perfectly phenotype-diagnostic column. Columns with fewer than five ungapped
residues in either group were skipped. Because diurnality is phylogenetically
clumped, the null was generated not by shuffling tip labels but by re-placing
clades of the observed sizes at random positions on the same gene tree, 1000
permutations per gene, so that the null retains the same autocorrelation
structure and only the phenotype assignment is randomised.

## Direct counting of convergent substitutions

Fitch parsimony was used to assign amino acids to internal nodes, so that a
substitution on a branch is simply a difference between parent and child states.
Because transition branches carry more substitutions of any kind than randomly
chosen branches of matched clade size, raw counts are confounded; the test
statistic was therefore the proportion of sites changing in two or more
independent lineages that changed to the same residue. Null branches were
additionally matched on branch length to within a factor of two. Parsimony ties
were resolved ten times per site with random tie-breaking to propagate ancestral
state uncertainty. Each site received a permutation p-value from 2000
permutations and a Benjamini-Hochberg q-value across all genes.

Sites were additionally required to be phenotype-specific: the frequency of the
convergent residue among diurnal species minus its frequency among nocturnal
species was required to exceed 0.25. Without this criterion the test identifies
homoplasy unrelated to diel activity, because a residue arising repeatedly across
the tree will by chance fall on some transition branches, and the real
transitions form a phylogenetically clustered set while the null scatters
branches more widely.

## Relative evolutionary rates

RERconverge v0.3.0 was used to test gene-level association between relative
evolutionary rate and diel activity, consuming the uniform-model fixed-topology
trees. Relative rates were computed with a square-root transform. The foreground
was defined as diurnal tip branches only, the conservative choice. Association
was tested with a binary phenotype correlation requiring at least 10 species and
2 foreground branches per gene, and p-values were corrected across the 18 genes.

## Selection analyses

Codon-level analyses were run in HyPhy v2.5.93 on the codon alignments. Branch
lengths were taken from the per-gene substitution trees rather than the
time-calibrated species tree, and all-gap codon columns were removed, with a
per-gene coordinate map retained so that post-filtering codon indices map back to
protein and human residue coordinates. Contrast-FEL tested, per site, whether the
ratio of nonsynonymous to synonymous substitution rates differs between the
diurnal-transition branch set and the remainder of the tree, with correction
applied both within gene and across genes. RELAX was run to test for relaxed or
intensified selection on the same branch set.

## Multi-method consensus

Per-site results from PCOC, TDG09 and Contrast-FEL were mapped through the
trimmed to untrimmed to codon to human residue coordinate chain and integrated,
with support from at least two methods required for a site to be considered high
confidence.

## Positive control

Because a null result is uninformative unless the analysis is shown capable of
detecting a signal, and because simulation under a method's own model cannot
detect an incorrectly constructed scenario, known convergent sites were planted
into the real alignments and the entire analysis was repeated unchanged over the
modified data, using identical trees, scenario files and scripts.

Sites were created by writing a residue absent at that column into the diurnal
species descending from a chosen set of gain events. Candidate columns were
required to be variable but not saturated (two to six residues present, at most
10 percent gaps). Ninety sites were planted across all 18 genes at three
difficulty levels: all 10 gain events converging, 7 of 10, and 3 of 10. Each site
was verified before retention: the planted residue had to be essentially absent
from nocturnal species, the diurnal minus nocturnal frequency gap had to exceed
0.10, and at least two independent events had to contribute. Recovery was scored
against a withheld key.

## Software and data availability

IQ-TREE v3.1.1; trimAl v1.5.rev1; corHMM v2.8; phytools v2.5.2; ape v5.8.1;
RERconverge v0.3.0; PCOC (`carinerey/pcoc` container); TDG09 v1.1.2; HyPhy
v2.5.93; pal2nal; miniprot; R v4.4.1; Python v3.11.3 with ete3 v3.1.3. Figures
were produced with Plotly and exported as vector PDF. Analysis code, scenario
files and result tables are available at the project repository.

---

# Results

## Evolutionary scenario

Diel activity shows significant phylogenetic signal (Pagel's lambda = 0.671,
p = 0.010). The equal-rates reconstruction recovers a nocturnal ancestral
placental mammal without that state being imposed, consistent with the nocturnal
bottleneck hypothesis, and infers **10 independent gains of diurnality and 5
reversals to nocturnality** across 15 transition branches, 10 terminal and 5
internal (Figure 1, Table 2).

The all-rates-different reconstruction instead places a diurnal ancestor at the
base of the placental radiation, rooting diurnality at a 54-tip clade, and
reframes the history as 6 gains and 10 reversals. This is the expected
consequence of the deliberately balanced diel sample and is the reason the
equal-rates model is treated as primary.

No transition branch falls on a branch that is weak on both concordance axes
(gCF below 50 percent and sCFL below 50 percent). The two closest cases are weak
on one axis only. Gene tree discordance therefore cannot plausibly generate false
convergence at any transition in this dataset.

## The analysis detects planted convergence

Ninety convergent sites planted into the real alignments were recovered as
follows (Table 3, Figure 2).

PCOC flagged exactly 50 of 11,727 sites on the modified data, and all 50 were
planted: 35 of 36 sites at which all 10 lineages converged, and 15 of 18 at 7 of
10, with **no false positives**. Signal therefore survives scenario construction,
tree handling and detection intact.

The phenotype-specific parsimony test recovered all 36 sites at 10 lineages, all
18 at 7 of 10, and 23 of 36 at only 3 of 10, with 26 false positives. Requiring
phenotype specificity reduced false positives from 195 to 26, an 87 percent
reduction, at no cost to recovery at the two higher levels.

PCOC recovered none of the 36 sites at 3 of 10 lineages, as its design predicts.
Sensitivity to convergence confined to few lineages therefore rests on the
parsimony test, which retains approximately 64 percent recovery at that level.

TDG09 recovered 36 of 36, 17 of 18 and 18 of 36 respectively, but produced **879
false positives** among 3,852 testable sites, a rate of 24.7 percent.

## No convergent signal accompanies gains of diurnality

Across all 18 genes, **648 sites changed in two or more independent diurnal
lineages, and 201 of those changed to the same residue**, an observed rate of
0.310 against a branch-length-matched null of 0.337 (Figure 3). No gene shows a
rate-conditioned p-value below 0.282.

Twelve of the 201 same-residue sites reached q <= 0.05, but none is
phenotype-specific. Nine of the twelve carry a residue as common in nocturnal as
in diurnal species; the clearest example carries the residue in two diurnal and
two nocturnal species, a frequency gap of -0.002, with two of the four carriers
being nocturnal marsupials. The largest diurnal minus nocturnal gap among the
twelve is 0.232 and the median is 0.067, so **no site is both statistically
unusual and phenotype-specific**.

The model-free residue screen agrees. Observed maximum diagnostic scores sit on
the permutation null in all 18 genes and below it in several (*RORB* 0.223
against 0.305; *CSNK1D* 0.250 against 0.326). Only 2 of 11,727 columns exceed a
diagnostic score of 0.6, against a comparable null expectation, and the lowest
p-value (*CRY1*, 0.043) does not survive correction across 18 genes.

PCOC returned **no site above threshold in any gene**. The highest posterior
observed across 11,727 columns is 0.098, and every other gene maximum is 0.003 or
below. Calibration on these scenarios gives power 1.000 and a false positive rate
of 0.0000 at every threshold tested, so this is not a marginal result.

## PCOC detects only near-universal convergence

Simulating convergence in only k of the 10 declared lineages and detecting with
the full scenario shows power rising from 0.000 at k = 2 and 3, through 0.012 at
k = 5 and 0.258 at k = 7, reaching 1.000 only at k = 10 (Figure 4). A direct
check confirms the mechanism: at k = 2, detection under the declared 10-event
scenario recovered 0 of 600 planted sites, while the same data detected under the
correct 2-event scenario recovered 600 of 600.

The controlling variable is the number of falsely declared events rather than the
quantity of sequence change available. Equal-branch comparisons separate them:
k = 2 with 25 branches and 8 falsely declared events gives power 0.000, while
k = 8 with 37 branches and 2 falsely declared events gives 0.545. PCOC's model
requires the derived profile on every declared branch, so each declared lineage
that did not converge penalises the fit.

The PCOC result should therefore be read as excluding convergence shared by
nearly all diurnal lineages, a narrower claim than excluding convergence
generally. The general claim rests on the parsimony and residue-screen results,
which declare no convergent set.

Notably, PCOC recovered 83 percent of planted sites at 7 of 10 lineages while
this curve gives 0.258 at k = 7. The two measurements plant different strengths
of signal: the simulation applies a probabilistic profile shift, whereas the
planted sites carry an unambiguous novel residue in every diurnal descendant.
PCOC's sensitivity to partial convergence therefore depends strongly on how clean
the convergent substitution is.

## Reversals to nocturnality

The reversal direction, tested here for the first time, returned **1 site of
11,727 above threshold** (*RORB*, trimmed column 1, posterior 0.99999963); every
other gene maximum is 0.119 or below and 11 of 18 are exactly zero.

That site is an alignment artefact. The trimmed alignment's first column maps to
untrimmed column 204, and 20 of the 60 sequences have their annotated protein
begin at exactly that position, making the column an initiator residue for a
third of the dataset and an internal residue for the rest. It carries 14 distinct
residues, whereas the five columns immediately following each carry a single
residue in 59 of 60 species. The 13 leaves descending from the five reversal
events carry eight different residues, giving a best diagnostic score of 0.137
against a null expectation of 0.305 for that gene.

PCOC's posterior decomposition identifies the mechanism. At this site the
profile-change component is 0.500, exactly chance, while the one-change component
is 0.9999994. The entire posterior derives from the term asking whether
substitutions occurred on the declared branches, which is trivially satisfied at
a hypervariable column, while the term testing for a shared derived amino acid
preference contributes nothing.

## Sensitivity to the ancestral-state model

Repeating the analysis under the all-rates-different reconstruction, for which
calibration gives power 0.990 to 1.000 and a false positive rate of 0.0000,
returned 2 sites of 11,727 above threshold, both in *RORB* and both in the same
ragged N-terminal region (columns 1 and 2, posteriors 0.998 and 0.9996). Both
carry a profile-change component of exactly 0.500 and are rejected on the same
grounds. All other gene maxima are 0.284 or below, and 12 of 18 are exactly zero.

This sensitivity analysis is weaker than the primary one: the ARD reconstruction
declares 55 of 60 leaves and 54 percent of branches convergent, against 41 leaves
and 47 percent under ER, leaving comparatively little ancestral contrast. The ARD
result is therefore consistent with the primary analysis rather than an
independent confirmation of it.

## Selection and rate association

Contrast-FEL identified **no site of 15,349** with a significant difference in
the nonsynonymous to synonymous rate ratio between the diurnal-transition branch
set and the remainder of the tree, at a 5 percent false discovery rate either
within gene or across genes; the smallest uncorrected p-value is 2.2e-4.

RELAX did not run reliably on this dataset. Nine of 18 genes converged and five
failed after three attempts each, all failing during ancestral reconstruction on
deeply divergent lineages. Among the nine, only *NPAS2* survives correction
(K = 0.352, q < 1e-5), but this result replicated once in eleven attempts and is
attributed to a flat likelihood surface rather than to relaxed selection. Because
the nine converged genes are the subset that happened to converge rather than a
random sample, no inference is drawn from this analysis.

RERconverge identified **no gene** whose relative evolutionary rate is
significantly associated with diel activity; all adjusted p-values exceed 0.67.
The strongest raw signal is *BHLHE40* (rho = -0.178, p = 0.051, adjusted
p = 0.677). The correlation sign is negative in 13 of 18 genes, indicating
marginally slower relative rates in diurnal lineages, the opposite of the
expectation under accelerated adaptive evolution, though far from significance.

## Cross-method integration

TDG09 flagged 885 of 3,820 testable sites at a 5 percent false discovery rate,
23.2 percent, a figure exceeding every other method by orders of magnitude. **No
flagged site is supported by any second method**, and the high-confidence
consensus set is empty (Figure 5). The positive control provides the direct
comparison: on modified alignments the same procedure produced false positives at
24.7 percent, statistically indistinguishable from its real-data flagging rate.
Its recovery of planted sites is high, so the limitation is specificity rather
than sensitivity.

## Genes with insufficient opportunity

Two genes cannot address the question. *CSNK1D* contains no position at which two
or more independent diurnal lineages changed, and *ARNTL* contains two, against
expectations of 2.0 and 2.1 usable sites respectively from branch lengths alone.
*CSNK1D* also has the lowest taxon occupancy (40 of 60). Their null results are
uninformative rather than negative and are reported separately throughout.

---

## Tables and figures

**Table 1. Gene set, taxon occupancy and alignment dimensions.**

Trimmed columns are those retained by trimAl `-automated1`. TDG09 testable sites are the variable columns for which a likelihood ratio could be computed. Codons analysed are those remaining after removal of all-gap columns.

| Gene | Module | Taxa | Untrimmed columns | Trimmed columns | TDG09 testable sites | Codons analysed |
|---|---|---|---|---|---|---|
| CLOCK | Positive arm | 59 | 991 | 841 | 177 | 991 |
| NPAS2 | Positive arm | 58 | 1145 | 814 | 351 | 1120 |
| ARNTL | Positive arm | 60 | 778 | 620 | 50 | 750 |
| PER1 | Negative arm | 60 | 1643 | 1233 | 409 | 1584 |
| PER2 | Negative arm | 60 | 1726 | 1208 | 801 | 1704 |
| PER3 | Negative arm | 56 | 2122 | 897 | 637 | 1640 |
| CRY1 | Negative arm | 59 | 670 | 582 | 80 | 666 |
| CRY2 | Negative arm | 60 | 665 | 592 | 211 | 661 |
| NR1D1 | Auxiliary loop | 60 | 760 | 610 | 144 | 731 |
| NR1D2 | Auxiliary loop | 60 | 660 | 562 | 156 | 660 |
| RORA | Auxiliary loop | 58 | 602 | 466 | 25 | 602 |
| RORB | Auxiliary loop | 60 | 665 | 459 | 34 | 609 |
| RORC | Auxiliary loop | 60 | 790 | 518 | 210 | 753 |
| CSNK1D | Post-translational | 40 | 792 | 401 | 12 | 553 |
| CSNK1E | Post-translational | 52 | 564 | 509 | 174 | 564 |
| FBXL3 | Post-translational | 60 | 442 | 426 | 34 | 442 |
| BHLHE40 | Output repressors | 60 | 666 | 407 | 132 | 599 |
| BHLHE41 | Output repressors | 58 | 768 | 582 | 183 | 720 |

**Table 2. Convergent scenarios, calibrated detection power and outcome.**

Branch and leaf counts are on the 60-taxon species tree (118 branches). Power and false positive rate are the worst value across the 18 genes at the calibrated posterior threshold of 0.99. All sites above threshold were subsequently rejected as alignment artefacts.

| Reconstruction | Direction | Events | Convergent branches | Percent of branches | Convergent leaves | Minimum gene power | Worst gene FPR | Sites above threshold |
|---|---|---|---|---|---|---|---|---|
| ER | Gains of diurnality | 10 | 55 | 46.6 | 41 | 1.000 | 0.0000 | 0 |
| ER | Reversals to nocturnality | 5 | 19 | 16.1 | 13 | 0.996 | 0.0000 | 1 |
| ARD | Gains of diurnality | 6 | 64 | 54.2 | 55 | 0.990 | 0.0000 | 2 |
| ARD | Reversals to nocturnality | 10 | 44 | 37.3 | 29 |  |  | not run |

**Table 3.** Positive control recovery by method and difficulty level, with false
positive counts. Data: `shared_results/spikein/spikein_scorecard.csv`.

**Figure 1.** Species tree of the 60 sampled mammals with diel activity states at
the tips and reconstructed transitions marked: gains of diurnality in blue,
reversals to nocturnality in orange. `results/figures/fig1_tree.pdf`.

**Figure 2.** Positive control recovery by method and difficulty level, with each
method's false positive count given in the legend.
`results/figures/fig2_positive_control.pdf`.

**Figure 3.** Observed against null same-residue rates per gene from the parsimony
analysis, genes ordered by the observed minus null difference.
`results/figures/fig3_parsimony_rates.pdf`.

**Figure 4.** PCOC power as a function of the number of converging lineages.
Existing: `results/figures/partial_convergence_CLOCK.pdf`.

**Figure 5.** Per-gene method tallies and consensus support.
Existing: `results/figures/consensus_by_gene.pdf`.

**Supplementary Figure S1.** PCOC power calibration per gene, three scenario
sets. Existing: `results/figures/pcoc_sim_power_{ER_gain,ER_reversal,ARD_gain}.pdf`.

**Supplementary Figure S2.** Selection analysis summary.
Existing: `results/figures/selection_summary.pdf`.

---

# Supplementary Methods S1. Distinguishing convergence from substitution rate in PCOC output

## The problem

PCOC reports a combined posterior probability that a site has undergone
convergent evolution on a declared set of branches. That posterior is a
combination of two components which the software also reports separately:

- **PC**, the profile-change component, which asks whether the site's amino acid
  preference profile shifted to a shared derived profile on the declared
  branches. This is the component that tests convergence in the sense usually
  intended.
- **OC**, the one-change component, which asks whether at least one substitution
  occurred on the declared branches.

The combined posterior can approach 1.0 when OC alone is near-certain and PC sits
at chance. Such a site has experienced substitutions on the convergent branches
but shows no shared derived preference, which is not convergence. Because the
combined posterior is what is conventionally thresholded, these sites pass
filtering.

This is not a rare edge case in practice. Every site that crossed a calibrated
posterior threshold of 0.99 in the present study, across three scenario sets and
35,181 site tests, had this signature.

## Worked example

*RORB* trimmed column 1 reached a posterior of 0.99999963 in the reversal
analysis, with PC = 0.500 and OC = 0.9999994. Inspection showed the column to be
an alignment artefact rather than a convergent site:

- The trimmed column maps to untrimmed column 204, and 20 of the 60 sequences
  have their annotated protein begin at exactly that position. The column is the
  initiator residue for a third of the dataset and an internal residue for the
  remainder.
- It carries 14 distinct amino acids, whereas the five columns immediately
  following each carry a single residue in 59 of 60 sequences.
- The 13 leaves descending from the declared convergent events carry eight
  different residues. The best achievable diagnostic score, defined as
  freq(residue | convergent) minus freq(residue | background), is 0.137, against
  a null expectation of 0.305 for that gene.

The column had been retained by trimAl because it contains only one gap. High
occupancy is not evidence of homology: those sequences are not aligned at that
position, they begin there.

## Recommended screen

We suggest three criteria, applied to any site passing a posterior threshold.
They are independent, and in our data no single criterion was sufficient.

1. **PC must not be at chance.** We used PC >= 0.8. A site whose posterior derives
   from OC alone should not be reported as convergent regardless of the combined
   value. This was the only criterion that rejected all three artefactual sites we
   encountered.
2. **The column must not be hypervariable.** We used at most 8 distinct residues.
   Convergence entails independent lineages arriving at the same residue, so a
   column with many states is more consistent with a ragged terminus or a
   misaligned region.
3. **The declared convergent leaves must share a residue.** We required a best
   diagnostic score of at least 0.20. Without this, a site can pass on the
   strength of substitution counts alone.

The value of applying all three is illustrated by the second artefactual site we
found, *RORB* column 2 in the ARD analysis, which had 8 distinct residues and a
diagnostic score of 0.278. Both fell inside the tolerances of criteria 2 and 3,
and only criterion 1 rejected it.

## Implementation

The screen is implemented as `check_pcoc_hit_credible` in `scripts/lib_checks.py`
and is exercised by four regression cases built from the real *RORB* column in
`scripts/test_lib_checks.py`.

## Scope

These thresholds were chosen against a single dataset of 18 genes and 60 taxa and
should be treated as a starting point rather than as calibrated values. The
general recommendation, that PC and OC be inspected separately rather than only
their combination, does not depend on the specific cutoffs.
