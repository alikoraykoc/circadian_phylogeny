# To-do list

Ordered by phase. `[x]` = done, `[ ]` = pending, `[?]` = decision needed.

## Phase 0 — inputs and setup
- [x] Species tree pruned from Upham (`data/species_tree.nwk`)
- [x] Per-gene protein alignments (18 genes, `data/alignments/`)
- [ ] Diel activity table `data/diel_activity.csv` (columns: species, activity)
      Source: Maor et al. 2017 / Bennie et al. 2014 for coding.
- [ ] Edit `config.sh`: exact `REF_SPECIES` tip label, gene basenames, extensions
- [ ] Build conda env; pull PCOC Docker; install RERconverge; place `tdg09.jar`
- [?] Do you have CDS (nucleotide) per gene? If yes -> selection track is live.
      If no -> defer track C; convergence + rate tracks run without it.

## Phase 1 — alignment prep
- [ ] `01_trim_alignments.sh` (trimAl `-automated1`, keep `-colnumbering` maps)
- [ ] Save trimmed-column -> human-residue mapping for final reporting
- [ ] `02_codon_alignments.sh` (pal2nal) — only if CDS available

## Phase 2 — trees
- [ ] `03_genetrees_qc.sh` — unconstrained ML per gene (QC + gCF input)
- [ ] MANUAL: inspect each QC tree; record keep / re-align / drop per gene
- [ ] `04_branchlengths_fixed.sh` — fixed-topology branch lengths (the backbone)
- [ ] SANITY: branch lengths in subs/site, no zero-length branches
- [ ] Re-root each per-gene tree to the Upham rooting
- [ ] `05_supermatrix.sh` — partitioned ML tree + constrained fallback lengths
- [ ] Compare supermatrix topology vs Upham, especially on transition branches
- [ ] `06_concordance.sh` — gCF / sCF on the Upham tree

## Phase 3 — phenotype scenario
- [?] Decide binary coding rule for cathemeral / crepuscular (edit `07_scenario.R`)
- [ ] `07_scenario.R` — ASR (corHMM), transition branches, PCOC scenario,
      HyPhy-labelled tree, Pagel's lambda
- [ ] Check transition count; if low (< ~4-5 independent events), revisit sampling
- [ ] `08_flag_discordant_branches.py` — implement bipartition matching, flag
      transition branches on low-gCF branches (hemiplasy risk)
- [?] Decide: exclude flagged branches from the scenario, or keep-and-caveat

## Phase 4 — analyses (parallel)
- [ ] `09_pcoc.sh` — run `pcoc_sim` FIRST for power/threshold, then `pcoc_det`
- [ ] Record the calibrated posterior threshold; use it (not a default 0.8)
- [ ] `10_tdg09.sh` — build the two-group file from the scenario, run per gene
- [ ] `11_rerconverge.R` — set foreground tips, run association, rank genes
- [ ] `12_selection_hyphy.sh` — Contrast-FEL + RELAX (only if codon alignments)

## Phase 5 — interpretation
- [ ] Ancestral sequence reconstruction (IQ-TREE `--ancestral` or FastML)
- [ ] Map candidate sites onto AlphaFold models; annotate functional domains
      (PAS-A/B, bHLH, CRY photolyase/FAD pocket, CSNK1 kinase sites, ROR/REV-ERB LBD)
- [ ] PHACT / PHACTboost: score functional impact of transition-branch substitutions
- [ ] PHACE: test coevolution among candidate sites and interacting partners

## Phase 6 — synthesis
- [ ] `13_consensus.py` — merge PCOC + TDG09 + Contrast-FEL, map to human residues,
      flag sites supported by >= 2 methods
- [ ] Multiple-testing correction across the site x gene space
- [ ] Figures as vector PDF (Plotly): consensus-by-gene, per-gene site tracks,
      annotated tree with transition branches
- [ ] Write up: which circadian module (core loop / stabilizing loop / kinases)
      carries the convergent signal, and the hemiplasy caveat where relevant

## Open decisions log
- [?] Convergence direction: derived diurnality, reversals to nocturnality, or both
- [?] Foreground definition for RERconverge (terminal vs ancestral branches)
- [?] gCF threshold for the hemiplasy flag (default 50%)
- [?] Whether to run aBSREL / BUSTED / MEME as a positive-selection add-on
