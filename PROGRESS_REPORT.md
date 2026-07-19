# Circadian Gene Convergent Evolution Analysis: Progress Report

**Date**: 2026-07-14 (updated; original draft 2026-07-07)
**Status**: Steps 00-08, 10, 11 complete. Step 09 (PCOC) primary run (ER / gains
of diurnality) complete; three sensitivity sweeps scaffolded and queued but not
run. See the update log immediately below before reading the older sections.

---

## 0. Update log, 2026-07-14 session (READ FIRST)

This session found and fixed a chain of defects in the ancestral reconstruction
and the PCOC scenario, then ran the corrected primary PCOC analysis. Several
claims in the older sections (5, 6, 8 in particular) were written before these
fixes and are now superseded. Where that is the case, the old text is kept for
traceability with a pointer to this section, not deleted.

### 0.1 What was wrong

1. **Ancestral state reconstruction inverted the project's premise.** Step 07 used
   corHMM with `model = "ARD"` and a free root. Refitting and inspecting the node
   states (which the old script never saved) showed it reconstructed the placental
   common ancestor as **diurnal**. Under that reconstruction, most of the 16
   "transitions" were reversals *to* nocturnality, not gains of diurnality. This
   contradicts the nocturnal-bottleneck premise in CLAUDE.md.

   Root cause: the taxon sample is exactly **30 diurnal / 30 nocturnal**. Real
   mammals are overwhelmingly nocturnal, so a balanced sample under an asymmetric
   model has no pull toward a nocturnal root. Model comparison confirms the issue:

   | ASR config | Placental ancestor | Gains | Reversals | AIC |
   |---|---|---|---|---|
   | ARD, free root (old) | diurnal | 6 | 10 | 79.27 |
   | **ER, free root (now primary)** | **nocturnal** | **10** | **5** | **77.82** |
   | ARD, root = nocturnal | nocturnal | 10 | 5 | 79.2 |
   | ER, root = nocturnal | nocturnal | 10 | 4 | 77.6 |

   The equal-rates (ER) model fits better (lower AIC, simpler) AND recovers the
   nocturnal ancestor without a root prior being imposed, so the bottleneck is a
   result rather than an assumption. ER is now the primary model; ARD is retained
   as a documented sensitivity analysis.

2. **PCOC scenario merged opposite phenotypes into one convergent class.** The old
   builder put every transition, gains and reversals together, into a single PCOC
   class. PCOC fits one convergent amino-acid profile per run, so branches moving
   toward opposite phenotypes cancel. Gains and reversals are now separate runs.

3. **PCOC scenario omitted the descendant branches of each event.** PCOC's `-m`
   format is `transition_node,other_convergent_nodes / next_event / ...`: the first
   node in a group is the transition branch and the rest are the branches that stay
   in the derived state. PCOC does NOT propagate the state downward. The old builder
   emitted one lone node per event, so for every internal transition the entire
   diurnal clade below the stem, including the tips carrying the convergent residues,
   was modelled as ancestral. That deletes exactly the signal PCOC looks for.

4. **PCOC scenario silently dropped events on pruned gene trees.** Node matching used
   exact tip-set equality; gene trees are pruned (CSNK1D has 40 of 60 taxa), so
   clades whose tip set no longer matched were dropped without warning (16 events
   became 14 groups).

   Also fixed earlier in the session (previous commit 3eb5dce): step 08 hemiplasy
   flagging was a silent no-op (wrong tree file, wrong column name, non
   rooting-invariant bipartitions), and step 06 gained likelihood-based site
   concordance (sCFL) alongside gCF.

### 0.2 Corrected state

- **Step 07 (`07_scenario.R`, rewritten):** fits ER (primary) and ARD (sensitivity),
  records each transition's DIRECTION, saves full node states, and for each event
  computes the contiguous clade of branches that remain in the derived state
  (descent stops at nested reversals). Primary ER reconstruction: nocturnal root
  (P = 0.546), **10 independent gains of diurnality (55 convergent branches)** and
  5 reversals (19 branches). Outputs per model: `node_states_*.csv`,
  `transition_branches_*.csv`, `convergent_events_*.csv`; `transition_branches.csv`
  is a copy of the ER version for backward compatibility.

- **Step 08 (hemiplasy flagging), re-run on the ER transitions:** 15 transitions
  (10 tip, 5 internal), **0 flagged** (criterion gCF < 50 AND sCFL < 50). Edges 4
  and 16 sit near the corner but are not low on both axes.

- **PCOC scenarios (`prep_pcoc_scenarios.py`, rewritten):** emit proper PCOC event
  groups, match nodes across trees by descendant-tip set, intersect each event with
  the taxa present in each gene, and report survivors per gene. Written under
  `results/pcoc/scenarios/<model>_<direction>/` for the four combinations
  ER_gain, ER_reversal, ARD_gain, ARD_reversal. Verified independently: the scenario
  marks all 30 diurnal tips convergent, 0 nocturnal tips, 0 diurnal tips missed;
  PCOC node numbering is reproducible against `pcoc_num_tree.py`.

- **PCOC runner (`run_pcoc_all.sh` / `run_pcoc_sets.sh`):** one scenario set per
  invocation, resumable (skips genes with an existing `results.tsv`), cooperative
  STOP flag at `results/pcoc/STOP`, and `-f 0.0` so the posterior for every site is
  retained and the significance threshold stays a post-hoc, `pcoc_sim`-calibrated
  choice rather than being baked into the run.

### 0.3 Primary PCOC result (ER, gains of diurnality)

**All 18 genes complete. Zero convergent sites.** No site in any gene reaches even
posterior 0.5. The single highest posterior anywhere is RORB site 1 at 0.30 (and
site 1 is the first alignment column, i.e. a likely edge artifact); every other gene
is at ~1e-3 or below. Blank rows in the results files are sites whose convergent-model
posterior underflowed to 0, i.e. strongly rejected, not untested; all 11,727 sites
were evaluated.

**Per-gene PCOC results (ER, gains of diurnality).** `scored` = sites with a
representable positive posterior; blanks are sites whose convergent-model posterior
underflowed to 0 (strongly rejected), so all `sites` were evaluated. `max_post` is
the single highest PCOC posterior in the gene.

| gene | sites | scored | max_post | >=0.5 | >=0.8 |
|------|------:|-------:|---------:|------:|------:|
| CLOCK | 841 | 598 | 0.0000 | 0 | 0 |
| NPAS2 | 814 | 756 | 0.0001 | 0 | 0 |
| ARNTL | 620 | 620 | 0.0000 | 0 | 0 |
| PER1 | 1233 | 760 | 0.0000 | 0 | 0 |
| PER2 | 1208 | 1107 | 0.0000 | 0 | 0 |
| PER3 | 897 | 631 | 0.0016 | 0 | 0 |
| CRY1 | 582 | 52 | 0.0000 | 0 | 0 |
| CRY2 | 592 | 529 | 0.0000 | 0 | 0 |
| NR1D1 | 610 | 582 | 0.0000 | 0 | 0 |
| NR1D2 | 562 | 561 | 0.0000 | 0 | 0 |
| RORA | 466 | 17 | 0.0000 | 0 | 0 |
| RORB | 459 | 66 | **0.3011** | 0 | 0 |
| RORC | 518 | 478 | 0.0001 | 0 | 0 |
| CSNK1D | 401 | 3 | 0.0000 | 0 | 0 |
| CSNK1E | 509 | 331 | 0.0000 | 0 | 0 |
| FBXL3 | 426 | 176 | 0.0000 | 0 | 0 |
| BHLHE40 | 407 | 407 | 0.0000 | 0 | 0 |
| BHLHE41 | 582 | 309 | 0.0000 | 0 | 0 |
| **TOTAL** | **11727** | **7983** | - | **0** | **0** |

(Low `scored` counts, e.g. CSNK1D 3/401, RORA 17/466, track how conserved the gene
is: an invariant column cannot show a convergent shift, so its posterior underflows.
These are the most constrained genes in the set, not failed runs.)

**The ER transition set feeding this run (15 branches, from `07_scenario.R`).**
Only the 10 `gain` events are in the ER_gain convergent class; the 5 reversals are a
separate (not-yet-run) class.

| edge | direction | from -> to | child | desc. tips |
|------|-----------|------------|-------|-----------:|
| 3 | gain | nocturnal -> diurnal | (internal) | 24 |
| 4 | reversal | diurnal -> nocturnal | (internal) | 2 |
| 16 | reversal | diurnal -> nocturnal | (internal) | 8 |
| 20 | gain | nocturnal -> diurnal | Suricata_suricatta | 1 |
| 27 | gain | nocturnal -> diurnal | Acinonyx_jubatus | 1 |
| 34 | reversal | diurnal -> nocturnal | Tapirus_terrestris | 1 |
| 46 | reversal | diurnal -> nocturnal | Tragelaphus_eurycerus | 1 |
| 60 | gain | nocturnal -> diurnal | (internal) | 8 |
| 67 | reversal | diurnal -> nocturnal | Aotus_nancymaae | 1 |
| 79 | gain | nocturnal -> diurnal | Ochotona_princeps | 1 |
| 82 | gain | nocturnal -> diurnal | (internal) | 4 |
| 89 | gain | nocturnal -> diurnal | Octodon_degus | 1 |
| 98 | gain | nocturnal -> diurnal | Arvicanthis_niloticus | 1 |
| 104 | gain | nocturnal -> diurnal | Rhynchocyon_petersi | 1 |
| 112 | gain | nocturnal -> diurnal | Myrmecobius_fasciatus | 1 |

This is a clean negative for the inputs as given. It is NOT yet interpretable as
"no convergence exists," for one specific reason:

**Power is uncalibrated, and there is a concrete reason to doubt it.** The scenario
marks 30 of 60 tips (half the tree) as convergent. PCOC is designed for convergent
taxa being a minority against an ancestral background; with half the tree convergent,
the ancestral and convergent profiles are each fit to half the data and the contrast
PCOC relies on may be too weak to detect anything. CLAUDE.md already requires the
posterior threshold to come from `pcoc_sim` calibration on the actual tree and
scenario; that step also answers the power question. Until it is run, the correct
statement is "no convergent sites detected, sensitivity unknown."

### 0.4 What is unaffected

- **RERconverge (step 11)** defines its foreground from tip states only
  (`clade = "terminal"`), never touching the ASR. Its results (section 7) stand.

### 0.5 What still needs doing (carried into section 10)

- **`pcoc_sim` power/FPR calibration FIRST**, before anything else. If power is poor,
  the three sensitivity sweeps would just reproduce an uninformative null three more
  times (~15 h of compute).
- **TDG09 (step 10) should be re-run.** `prep_tdg09_inputs.py` reconstructs internal
  node states by flipping at each transition branch starting from a nocturnal root;
  it consumed the pre-fix (ARD, diurnal-ancestor) transition set, so its internal
  branch group assignments are inconsistent with the corrected ER reconstruction.
  Tip (terminal) group labels come straight from `diel_activity.csv` and are fine;
  the internal-branch labels are the concern. The section 6 numbers below predate
  the fix and should be regenerated before use.
- Resume the sensitivity sweeps once power is confirmed: `rm results/pcoc/STOP` then
  relaunch `run_pcoc_sets.sh` (it auto-skips the finished ER_gain set).

### 0.6 Deprecated artifacts

The first (pre-fix) PCOC run is archived at
`results/pcoc/_deprecated_broken_scenario/` with a README explaining why it is
invalid (it reported 1 site above 0.8 across all genes, an artifact of the merged
class and the omitted descendants). Do not use it.

---

## 1. Project overview

This project tests whether convergent molecular evolution in 18 mammalian circadian
clock genes is associated with independent evolutionary transitions in diel activity
pattern. Mammals are ancestrally nocturnal; diurnality has evolved independently
multiple times across the mammalian phylogeny. If this behavioral convergence is
accompanied by convergent changes at the molecular level, specific amino acid
positions in circadian genes should show correlated substitution patterns on the
branches where diel transitions occurred.

### Species sampling

60 mammalian species were sampled (30 diurnal, 30 nocturnal), spanning all major
mammalian orders. Species were selected to maximize the number of independent
transitions between activity states while retaining phylogenetic balance. The
species tree is pruned from Upham et al. (2019), which provides the fixed
topological backbone for all analyses.

### Gene set

18 core circadian clock genes spanning all major functional modules:

| Module | Genes | Function |
|--------|-------|----------|
| Positive arm (TTFL activators) | CLOCK, NPAS2, ARNTL | CLOCK-BMAL1 heterodimer drives transcription |
| Negative arm (period/cryptochrome) | PER1, PER2, PER3, CRY1, CRY2 | Repress CLOCK-BMAL1 activity |
| Auxiliary loop (REV-ERB/ROR) | NR1D1, NR1D2, RORA, RORB, RORC | Stabilize oscillation via RORE elements |
| Post-translational regulators | CSNK1D, CSNK1E, FBXL3 | Phosphorylation and degradation timing |
| Output/light response | BHLHE40, BHLHE41 | DEC1/DEC2: light-inducible repressors |

---

## 2. Pipeline architecture

The analysis follows a three-track design that converges at the consensus step:

```
Track 1: CONVERGENCE (site-level)
  PCOC (profile change) + TDG09 (fitness shift)
  -> per-site posterior probabilities / LRT statistics

Track 2: RATE ASSOCIATION (gene-level)
  RERconverge (relative evolutionary rates)
  -> per-gene correlation with diel phenotype

Track 3: SELECTION (site-level, deferred)
  Contrast-FEL + RELAX (requires CDS sequences)
  -> per-site evidence for differential selection

                   |
                   v
          CONSENSUS (step 13)
  Sites supported by >= 2 methods are high-confidence
```

### Key design decisions

1. **Three trees, three roles**: The Upham species tree is the fixed scaffold for
   every analysis. Per-gene trees use the Upham topology with gene-specific branch
   lengths (not independent ML searches). Unconstrained gene trees are used only
   for QC (orthology checks, concordance factors).

2. **Two-track model strategy**: PCOC and TDG09 get per-gene best-fit models (MFP)
   because they analyze genes independently and need accurate branch lengths.
   RERconverge gets a single uniform model (Q.MAMMAL+F+R6, selected by supermatrix
   MFP) because it compares rates across genes, and different rate-heterogeneity
   models would introduce shape confounds in cross-gene comparisons.

3. **Hemiplasy control**: Diel transitions on low gene-concordance-factor (gCF)
   branches are flagged, because gene-tree/species-tree discordance at those
   branches can generate false convergence signals via hemiplasy.

---

## 3. Data preparation results (steps 00-01)

### Input validation (step 00)
- All 18 protein alignments confirmed present
- Species tree has 60 tips
- One label mismatch detected and fixed: `Gorilla_gorilla_gorilla` in alignments
  vs `Gorilla_gorilla` in the tree. Corrected across all 18 alignment files.

### Alignment trimming (step 01)
- trimAl `-automated1` heuristic applied to all 18 alignments
- Column numbering maps saved for downstream human-residue mapping
- Alignment lengths after trimming range from ~400 (CSNK1D) to ~1,200 (PER2) columns

---

## 4. Tree construction results (steps 03-05)

### Gene tree QC (step 03)
- Unconstrained ML gene trees built for all 18 genes
- Used for orthology verification and concordance factor computation
- No paralogy or contamination issues detected

### Supermatrix analysis (step 04)

Three products:

**(a) Unconstrained partitioned ML tree**: Strongest in-house topology estimate.
Robinson-Foulds distance to Upham: 17.5% (10 conflicting bipartitions out of 57).
All conflicts occur at known difficult nodes in mammalian phylogenetics:

- Laurasiatheria internal relationships
- Afrotheria-Xenarthra placement
- Some rodent internal arrangements

No unexpected conflicts were found, confirming the Upham topology is appropriate
for this dataset.

**(b) Concatenated MFP model selection**: Single-model MFP on the full concatenation
selected **Q.MAMMAL+F+R6** as the best-fit model under BIC. This is the uniform
model used for all RERconverge gene trees.

- Q.MAMMAL: an empirical amino acid replacement matrix derived from mammalian
  protein evolution (better fit than WAG or LG for this dataset)
- +F: empirical state frequencies from the data
- +R6: free-rate model with 6 categories (more flexible than discrete gamma)

**(c) Constrained fallback tree**: Upham topology with concatenated branch lengths
under Q.MAMMAL+F+R6, available as a fallback if any per-gene tree estimation fails.

### Fixed-topology branch lengths (step 05)

For each of the 18 genes, two tree sets were estimated on the fixed Upham topology:

- **PCOC/TDG09 trees** (per-gene MFP): Each gene gets its own best-fit substitution
  model and branch lengths. These are used by convergence detection methods that
  analyze genes independently.

- **RERconverge trees** (uniform Q.MAMMAL+F+R6): All genes estimated under the same
  model. This ensures that cross-gene rate comparisons in RERconverge are not
  confounded by model-driven differences in branch length ratios.

All trees confirmed to have small-decimal branch lengths (substitutions/site),
no zero-length branches, and correct rooting matching the Upham tree.

---

## 5. Support and hemiplasy control (steps 06-08)

### Gene concordance factors (step 06)

Gene concordance factor (gCF) measures what fraction of individual gene trees
support each branch in the species tree. Site concordance factor does the same at
the site level. Step 06 now computes the likelihood-based site concordance factor
(**sCFL**, a second IQ-TREE run merged by branch ID) rather than the deprecated
parsimony sCF; `--scfl` cannot be combined with `--gcf`, hence two runs. Step 08
uses both.

**Results**:
- All 56 internal branches have gCF < 50% (maximum 16%, median 9%)
- This is expected with only 18 protein loci; gene concordance factors require
  hundreds of loci to consistently exceed 50%
- The values are still informative for **relative** hemiplasy risk: branches with
  gCF of 3% are much more suspect than those with 16%

### Ancestral state reconstruction (step 07)

> SUPERSEDED, see section 0.1. This described the ARD free-root reconstruction,
> which placed a DIURNAL placental ancestor and mixed 6 gains with 10 reversals.
> The primary reconstruction is now ER (nocturnal ancestor, 10 gains, 5 reversals).
> The lambda result below is unchanged; the transition framing is not.

corHMM (all-rates-different model) was used to reconstruct the history of diel
activity transitions on the species tree.

**Key results**:

- **Pagel's lambda = 0.671 (p = 0.010)**: Diel activity has significant
  phylogenetic signal. It is not randomly distributed across the tree; closely
  related species tend to share activity patterns. This is a prerequisite for
  the convergent evolution analyses to be meaningful.

- **16 transition branches identified**: These are the branches where the
  ancestral state reconstruction infers a shift between nocturnal and diurnal
  (or vice versa). 16 independent events is a strong convergence scenario;
  PCOC typically needs at least 4-5 transitions to have adequate statistical
  power.

- The transitions include both gains of diurnality (the primary convergent
  phenotype) and reversals back to nocturnality, as expected from the complex
  evolutionary history of mammalian activity patterns.

### Hemiplasy flagging (step 08)

Hemiplasy occurs when incomplete lineage sorting causes a gene tree to differ
from the species tree, making ancestral polymorphism look like convergent
evolution. This is most likely on short internal branches where gene-tree
discordance is high.

**Method**: Each transition branch was matched to its concordance factor by
aligning descendant-tip bipartitions between the ape (R) node numbering system
and IQ-TREE's internal numbering (canonicalized to be rooting-invariant).

**Results (corrected, see section 0.2):**
- The original version of this step was a silent no-op (three bugs fixed earlier
  in the session; see section 0.1). The numbers previously reported here ("1 of 16
  flagged, Atlantogenata gCF 5.6%") came from the broken run and are withdrawn.
- Re-run on the corrected ER transition set: **15 transitions (10 tip, 5 internal),
  0 flagged** under the criterion gCF < 50% AND sCFL < 50%. A branch is flagged only
  when weak on BOTH axes, since either alone is usually estimation noise. Missing
  (NA) values never flag.
- Edges 4 (gCF 56.2 / sCFL 48.1) and 16 (gCF 50.0 / sCFL 42.0) sit near the corner
  but are not low on both axes, so no diel transition is on a branch where gene-tree
  discordance could plausibly manufacture false convergence.

---

## 6. TDG09 results (step 10): site-specific amino acid fitness shifts

> CAVEAT (see section 0.5): these results predate the ASR fix. TDG09's internal
> branch group assignments were derived from the superseded ARD/diurnal-ancestor
> transition set and should be regenerated on the ER reconstruction before use.
> Tip (terminal) group labels come from `diel_activity.csv` and are unaffected;
> internal-branch labels are the concern. The numbers below are retained for
> continuity, not as final results.

### Method

TDG09 (Tamuri, Dos Reis, and Goldstein 2009) tests whether the amino acid fitness
landscape at each alignment site differs between two groups of lineages (here:
diurnal vs nocturnal). It fits two models per site:

- **WAG+ssF** (site-specific frequencies): A single set of amino acid frequencies
  for the entire tree.
- **WAG+lssF** (lineage-specific site frequencies): Separate amino acid frequencies
  for diurnal and nocturnal branches.

The likelihood ratio test (LRT = 2 * [lnL_lssF - lnL_ssF]) measures whether
allowing different frequencies significantly improves the fit.

### Results summary

| Gene | Total sites | Variable sites | LRT > 3.84 (p<0.05) | LRT > 6.63 (p<0.01) |
|------|-------------|----------------|----------------------|----------------------|
| CLOCK | 841 | 177 | 57 | 50 |
| NPAS2 | 814 | 351 | 136 | 116 |
| ARNTL | 620 | 50 | 28 | 26 |
| PER1 | 1,233 | 409 | 184 | 154 |
| PER2 | 1,208 | 801 | 403 | 303 |
| PER3 | 897 | 637 | 287 | 219 |
| CRY1 | 582 | 80 | 45 | 31 |
| CRY2 | 592 | 211 | 38 | 32 |
| NR1D1 | 610 | 144 | 73 | 62 |
| NR1D2 | 562 | 156 | 61 | 51 |
| RORA | 466 | 25 | 15 | 14 |
| RORB | 459 | 34 | 12 | 9 |
| RORC | 518 | 210 | 92 | 70 |
| CSNK1D | 401 | 12 | 4 | 4 |
| CSNK1E | 509 | 174 | 99 | 74 |
| FBXL3 | 426 | 34 | 14 | 12 |
| BHLHE40 | 407 | 132 | 51 | 48 |
| BHLHE41 | 582 | 183 | 89 | 82 |
| **Total** | **11,727** | **3,820** | **1,688** | **1,357** |

### Top 10 sites across all genes

| Rank | Gene | Site | LRT | Residue pattern |
|------|------|------|-----|-----------------|
| 1 | PER2 | 860 | 211.96 | A:15, L:8, P:2, S:6, T:21 |
| 2 | PER2 | 855 | 209.30 | A:18, N:5, P:6, T:23, V:2 |
| 3 | CSNK1E | 387 | 191.75 | A:9, R:4, G:2, S:34, V:2 |
| 4 | PER3 | 707 | 164.67 | N:10, D:28, E:9, G:2, S:4 |
| 5 | PER3 | 678 | 143.09 | A:3, I:7, L:15, T:2, V:24 |
| 6 | CRY2 | 17 | 121.68 | R:6, C:2, G:30, P:5, S:11 |
| 7 | PER1 | 685 | 111.40 | A:15, G:2, P:8, S:2, T:32 |
| 8 | PER3 | 892 | 108.52 | A:5, I:16, L:5, V:29 |
| 9 | PER2 | 155 | 103.88 | N:18, G:8, K:4, S:27 |
| 10 | BHLHE40 | 191 | 96.58 | A:7, G:5, L:4, P:2, S:42 |

### Interpretation

1. **PER2 dominates the signal**: 403 nominally significant sites, including the two
   strongest hits genome-wide (LRT > 200). PER2 is the longest gene in the dataset
   (1,208 columns) with 66% variable sites, but the signal density is still
   disproportionately high. PER2 is the primary period gene and a core component
   of the negative feedback loop; it is plausible that activity-pattern shifts
   place distinct selective pressures on PER2 residues involved in protein-protein
   interactions (CRY binding, CK1 phosphorylation sites).

2. **PER3 and CSNK1E are strongly enriched**: PER3 has three sites in the overall
   top 10. CSNK1E site 387 (LRT = 191.75) is the third-strongest hit; CK1epsilon
   phosphorylates PER proteins to control their nuclear entry timing, and variation
   in this kinase could directly modulate circadian period length. The PER3 + CSNK1E
   co-enrichment is biologically coherent: they are in the same functional pathway.

3. **Highly conserved genes show minimal signal**: CSNK1D (4 significant sites,
   only 12 variable columns) and RORA (15 significant, 25 variable) are under
   extremely strong purifying selection. CSNK1D and CSNK1E are paralogs, but CSNK1D
   is far more constrained, consistent with its broader role beyond circadian timing
   (Wnt signaling, DNA repair).

4. **Auxiliary loop genes (RORC) show unexpected enrichment**: RORC has 92 significant
   sites despite being a 518-column gene. ROR-gamma is primarily known for immune
   function (Th17 cell differentiation), but its circadian role in peripheral tissues
   may be under-appreciated. This could reflect tissue-specific adaptation in
   diurnal species.

5. **Important caveat**: These are raw, uncorrected p-values. The 44% hit rate among
   variable sites (1,688/3,820) will shrink substantially after FDR correction and
   especially after requiring confirmation by a second method (PCOC) at the
   consensus step. The current numbers represent the upper bound of the signal.

---

## 7. RERconverge results (step 11): gene-level rate association

### Method

RERconverge tests whether the relative evolutionary rate of each gene is correlated
with a binary phenotype (diurnal vs nocturnal) across the phylogeny. It computes
relative evolutionary rates (RERs) by regressing each gene's branch lengths against
the genome-wide average, then correlates these residuals with the phenotype using a
Kendall rank correlation that accounts for phylogenetic non-independence.

A significant negative Rho means diurnal lineages evolve the gene more slowly
(stronger purifying selection). A positive Rho means diurnal lineages evolve it
faster (relaxed constraint or positive selection).

### Results

| Gene | Rho | N (branches) | Raw P | Adjusted P |
|------|-----|--------------|-------|------------|
| BHLHE40 | -0.178 | 82 | 0.051 | 0.677 |
| NPAS2 | -0.124 | 103 | 0.127 | 0.677 |
| ARNTL | -0.149 | 68 | 0.137 | 0.677 |
| CSNK1D | +0.204 | 34 | 0.156 | 0.677 |
| BHLHE41 | -0.105 | 80 | 0.257 | 0.677 |
| PER2 | -0.085 | 114 | 0.270 | 0.677 |
| NR1D2 | -0.092 | 82 | 0.311 | 0.677 |
| CLOCK | -0.080 | 85 | 0.370 | 0.677 |
| CRY2 | -0.077 | 92 | 0.374 | 0.677 |
| NR1D1 | -0.070 | 89 | 0.426 | 0.677 |
| FBXL3 | +0.076 | 67 | 0.453 | 0.677 |
| CRY1 | +0.061 | 77 | 0.514 | 0.677 |
| RORB | -0.061 | 71 | 0.534 | 0.677 |
| PER1 | -0.048 | 101 | 0.555 | 0.677 |
| RORA | +0.068 | 49 | 0.571 | 0.677 |
| RORC | -0.042 | 105 | 0.602 | 0.677 |
| CSNK1E | +0.029 | 84 | 0.744 | 0.788 |
| PER3 | -0.005 | 103 | 0.952 | 0.952 |

### Interpretation

1. **No gene reaches significance after multiple-testing correction** (all adjusted
   P > 0.67). This is the central result: there is no evidence for wholesale
   acceleration or deceleration of any circadian gene in diurnal vs nocturnal
   lineages.

2. **BHLHE40 is the strongest candidate** (Rho = -0.178, raw P = 0.051). The
   negative correlation means diurnal lineages evolve BHLHE40 (DEC1) more slowly,
   suggesting stronger purifying selection. BHLHE40 is a light-inducible
   transcriptional repressor that participates in resetting the clock after light
   exposure. Stronger constraint in diurnal species is biologically plausible:
   organisms active during the day rely more heavily on light-mediated clock
   entrainment, so the light-responsive pathway may be under tighter functional
   constraint. BHLHE40 also ranks in the TDG09 top 10 (site 191, LRT = 96.58),
   suggesting that while the gene overall is not accelerated, specific sites within
   it may show convergent shifts.

3. **The negative trend is consistent across most genes**: 13 of 18 genes show
   negative Rho values (slower evolution in diurnal lineages). This is not
   individually significant, but the pattern is consistent with a mild global
   increase in purifying selection on the circadian system in diurnal species,
   who may depend more tightly on precise clock function due to their reliance
   on daylight cues.

4. **CSNK1D is the only gene with a notable positive correlation** (Rho = +0.204,
   P = 0.156), meaning diurnal lineages evolve it faster. However, N = 34
   (many branches dropped due to missing species) makes this estimate unreliable.

5. **Why no significance?** With only 18 genes and ~60 species, RERconverge has
   limited power to detect gene-level rate shifts. The method is designed for
   genome-wide scans with hundreds to thousands of genes. Our result does not
   rule out convergent evolution in these genes; it specifically rules out
   large-scale rate changes. The convergent signal, if present, is at the
   site level (as TDG09 suggests), not the gene level.

---

## 8. Pending analyses

### PCOC (step 09): profile-based convergence detection

PCOC (Rey et al. 2018) tests whether amino acid preference profiles shift on
convergent branches. Unlike TDG09 (which compares frequency distributions),
PCOC uses a mixture-model framework (C10/C60 profiles from Le et al.) to
estimate the posterior probability that each site underwent a convergent
profile change.

**Status (updated, see section 0.3): primary run complete, result is a null.** The
ER / gains-of-diurnality set ran all 18 genes on the M1 under Rosetta (~4 min to
~74 min per gene; runtime does not track alignment length). Zero convergent sites at
any threshold; highest posterior anywhere is RORB site 1 at 0.30. The three
sensitivity sweeps (ER_reversal, ARD_gain, ARD_reversal) are scaffolded with scenarios
built, but were intentionally held (STOP flag) pending the power question below. See
section 0.1 for the four scenario bugs that invalidated the first attempt, and
`results/pcoc/_deprecated_broken_scenario/` for that archived run.

**Before trusting the null: run `pcoc_sim`.** The scenario marks half the tree as
convergent, which may leave PCOC underpowered (section 0.3). `pcoc_sim` on this tree
and scenario gives both the calibrated posterior threshold and the power estimate
that says whether the null is meaningful.

**Why PCOC matters**: TDG09 and PCOC test overlapping but distinct hypotheses. TDG09
asks "do the amino acid frequencies differ between groups?" (a frequentist LRT).
PCOC asks "did the amino acid preference profile change on the convergent branches?"
(a Bayesian posterior). Sites significant in both methods are the strongest
candidates for genuine convergent evolution, as they pass two independent
statistical frameworks.

### Selection track (step 12): differential selection pressure

Contrast-FEL and RELAX from HyPhy test whether specific sites or genes are under
different selection regimes in diurnal vs nocturnal lineages. This requires
codon (nucleotide) alignments, which are not yet available.

**Plan**: Retrieve CDS sequences using the same transcript accessions as the
protein alignments, verify translation matches, build codon alignments with
pal2nal. Run only on candidate genes identified by the convergence track.

### Consensus (step 13): multi-method integration

The final step merges per-site evidence across all methods:
- PCOC posteriors (above calibrated threshold)
- TDG09 LRTs (above significance after FDR correction)
- Contrast-FEL results (if available)

Sites supported by >= 2 independent methods are classified as high-confidence
convergent sites. These are then mapped to human residue numbering using the
trimAl column maps and the Homo_sapiens alignment row, enabling structural
interpretation.

---

## 9. Preliminary biological picture

Even before PCOC and the consensus step, several patterns are emerging:

1. **The convergence signal is site-level, not gene-level**. RERconverge finds no
   significant gene-wide rate shifts, but TDG09 finds hundreds of sites with
   amino acid fitness differences between diurnal and nocturnal lineages. This is
   consistent with stabilizing selection on overall protein function combined with
   adaptive or compensatory changes at specific residues.

2. **The PER genes are the primary targets**. PER2 and PER3 together account for
   690 of 1,688 nominally significant TDG09 sites (41%). The Period proteins are
   the rate-limiting step in the circadian negative feedback loop; changes in their
   interaction surfaces with CRY, CK1, and nuclear transport machinery could shift
   the phase or period of the clock in ways relevant to activity timing.

3. **The PER-CK1 axis shows coordinated signals**. CSNK1E (CK1-epsilon) has
   the third-strongest individual site hit (LRT = 191.75), and PER2/PER3 are the
   two most enriched genes. Since CK1-epsilon phosphorylates PER proteins to
   regulate their stability and nuclear entry, convergent changes in both the
   kinase and its substrates could represent co-evolution of the phospho-timer
   that sets circadian period.

4. **BHLHE40 bridges the gene-level and site-level results**. It is the only gene
   approaching significance in RERconverge (marginally slower evolution in diurnal
   lineages) and also has a top-10 TDG09 site. As the primary light-responsive
   repressor, it is positioned at the interface between environmental light input
   and the core oscillator, exactly where diel activity pattern would be expected
   to exert selective pressure.

5. **One transition branch (Atlantogenata) is suspect**. Any convergent sites
   driven primarily by this branch should be flagged in the consensus step.
   The remaining 15 transitions are on well-supported branches.

---

## 10. Next steps (priority order)

Revised for the corrected state (see section 0.5). The first item gates the rest:
if PCOC is underpowered on this scenario, the convergence track needs a design
change (e.g. taxon resampling toward the natural nocturnal majority), not more runs.

1. **Run `pcoc_sim` power/FPR calibration** on the actual ER tree and gains
   scenario. This decides whether the ER_gain null is meaningful and sets the
   posterior threshold. Do this before item 2.
2. If power is adequate: resume the three sensitivity sweeps
   (`rm results/pcoc/STOP`; relaunch `run_pcoc_sets.sh`, which skips ER_gain).
   If power is poor: reconsider taxon sampling (the 30/30 balance is the root
   issue) rather than running more sweeps.
3. **Re-run TDG09 (step 10)** on the corrected ER reconstruction, then parse with
   FDR correction (its current numbers in section 6 predate the fix).
4. Run consensus (step 13) to identify sites supported by >= 2 methods. Note: with
   a PCOC null, the current consensus would rest on TDG09 alone, which is why items
   1 and 3 come first.
5. Retrieve CDS for top candidate genes and run the selection track (steps 02, 12).
6. Map high-confidence sites to human residue numbering and onto protein structures
   (PDB/AlphaFold) for interpretation against domains, interaction surfaces, and PTM
   sites.

### Handoff notes for the next person

- Environment: conda env `phylo`; PCOC via Docker (`carinerey/pcoc`, runs under
  Rosetta on Apple Silicon). `results/` is gitignored and regenerated by the
  pipeline; pull the scripts and re-run, do not expect results in the repo. The
  key result numbers are embedded in section 0.3 so they survive the handoff.
- Run order and per-step detail are in README.md; open decisions in TODO.md.
- The `07_scenario.R` -> `prep_pcoc_scenarios.py` -> `run_pcoc_sets.sh` chain is the
  corrected convergence-scenario path. `08_flag_discordant_branches.py` is the
  hemiplasy control. These four plus `06_concordance.sh` are what changed this
  session.
- Watch the 30/30 taxon balance: it is the upstream cause of both the ASR
  model-sensitivity and the half-the-tree-convergent power concern.
