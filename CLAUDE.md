# CLAUDE.md

Project context for Claude Code. Read this before editing anything.

## Project

Testing whether convergent molecular evolution in 18 mammalian circadian genes is
associated with independent transitions in diel activity (derived diurnality, plus
reversals to nocturnality). The output is a set of high-confidence convergent sites
and phenotype-associated genes, interpreted against protein structure and function.

Genes: CLOCK NPAS2 ARNTL PER1 PER2 PER3 CRY1 CRY2 NR1D1 NR1D2 RORA RORB RORC
CSNK1D CSNK1E FBXL3 BHLHE40 BHLHE41.

## Environment Constraints

When installing software that requires `sudo` (e.g., Docker, brew cask installs),
provide the exact commands for the user to run in their external terminal rather
than attempting to run them directly. Claude's terminal cannot handle interactive
sudo prompts.

## Pipeline Execution

For long-running bioinformatics tools (Docker containers, RERconverge, PCOC), check
if the process is still running before assuming failure. Use `docker ps`,
`ps aux | grep`, or log file tailing to monitor progress.

## Status snapshot

Ready:
- Species tree pruned from Upham et al. 2019 (`data/species_tree.nwk`).
- Per-gene protein alignments for all 18 genes (`data/alignments/`).
- Full pipeline scaffold: 14 numbered scripts, config, conda env.

Pending:
- `data/diel_activity.csv` not yet built (source: Maor et al. 2017 / Bennie et al. 2014).
- No CDS (nucleotide) sequences yet. Selection track is deferred (see below).
- Two scripts are deliberate skeletons awaiting real implementation: `08` and `13`.

## Key design decisions (do NOT change without asking the user)

1. Three trees, three distinct roles:
   - Upham species tree is the FIXED scaffold for every analysis. The diel scenario
     is mapped onto it once so transition branches sit at identical nodes across genes.
   - Per-gene trees are the Upham topology with gene-specific branch lengths
     (subs/site), produced with IQ-TREE `-te`. NOT independent ML searches.
   - Unconstrained gene trees + supermatrix tree are SUPPORT tools only (orthology
     QC, gene-concordance factors, topology confirmation, fallback branch lengths).
     Never used as the analysis backbone.
2. Convergence direction: mammals are ancestrally nocturnal, so the convergent state
   is derived diurnality (and reversals to nocturnality). Transition branches are read
   off the ASR, never hand-coded from tip states.
3. One shared substitution model (LG+G4) across all fixed-topology gene trees, so
   RERconverge's cross-gene rate comparison is not confounded by model differences.
4. Three analysis tracks run in parallel, then reconverge: convergence (PCOC + TDG09),
   rate signal (RERconverge), selection (Contrast-FEL + RELAX). A site is high
   confidence only if backed by >= 2 methods.
5. Hemiplasy control: diel transitions that fall on low-gCF branches are flagged,
   because gene-tree/species-tree discordance there can manufacture false convergence.

## Repo layout

```
config.sh                 shared vars (paths, GENES, REF_SPECIES, model, threads)
data/                     species_tree.nwk, diel_activity.csv, alignments/, cds/
scripts/00..13            numbered pipeline steps (see README run order)
results/                  all outputs (created by config.sh)
env/environment.yml       conda env
```

Run order and per-step descriptions are in README.md. TODO.md is the working
checklist with the open decisions log.

## Immediate next tasks (priority order)

1. `scripts/08_flag_discordant_branches.py`: implement bipartition matching. Match
   transition branches (from corHMM, step 07) to IQ-TREE concordance branches
   (step 06) by descendant-tip sets, since both reference the same Upham topology.
   Join gCF/sCF, flag transition branches with gCF below threshold. Use ete3.
   This is the load-bearing hemiplasy control.
2. CDS retrieval script (new, for the selection track): fetch CDS by the SAME
   transcript/accession ID used for each protein (not by gene name), verify each CDS
   translates to its alignment row, drop or report mismatches. Isoform matching is the
   only real pitfall; handle it up front. Recommended sequencing: run PCOC + RERconverge
   first, then pull CDS ONLY for candidate genes and run selection on those.
3. `scripts/13_consensus.py`: implement the per-tool parsers (PCOC posteriors, TDG09
   LRTs, Contrast-FEL JSON) and the trimmed-column to human-residue mapping using the
   `*.colnumbering.txt` files plus the REF_SPECIES alignment row.

## Conventions and constraints (must follow)

- Visualization: Plotly only, never matplotlib.
- All figures: vector PDF, never PNG. (Plotly `write_image(..., ".pdf")` via kaleido.)
- No em dashes anywhere, including code comments and docs. Use commas, semicolons,
  or restructured sentences.
- Prefer automated command-line pipelines. Python for stats and visualization; R only
  where a method requires it (corHMM for ASR, RERconverge).
- Never assume tip-label consistency; verify alignment headers match tree tips exactly.
  Mismatched labels are the top cause of IQ-TREE `-te` failures.

## Environment and external tools

```
conda env create -f env/environment.yml && conda activate circadian-conv
docker pull carinerey/pcoc                          # PCOC (not on conda)
R: remotes::install_github("nclark-lab/RERconverge") # RERconverge (from GitHub)
place tdg09.jar in tools/                            # TDG09 (Java, from its release)
```

## Invariants and gotchas

- Per-gene branch lengths must be small decimals (subs/site), NOT the million-year
  scale of the Upham tree. No zero-length branches. Re-root each per-gene tree to the
  Upham rooting before scenario building.
- PCOC posterior threshold comes from `pcoc_sim` power/FPR calibration on the actual
  tree and scenario, not a fixed default 0.8.
- Selection track (steps 02, 12) needs codon alignments from pal2nal, which need CDS.
  Convergence and rate tracks do not; they run on the protein data as-is.
- If the transition count from step 07 is low (< ~4-5 independent events), PCOC power
  is limited; that is a signal to revisit taxon sampling, not to push on.
- CLI flags for PCOC / TDG09 / HyPhy vary by version; check `-h` against the installed
  build rather than trusting the scaffolded flags verbatim.

## Grounding references

PCOC method: Rey et al. 2018 (MBE). Mammalian activity data and ancestral nocturnality:
Maor et al. 2017, Bennie et al. 2014. Hemiplasy and false convergence: Mendes and Hahn.
RERconverge: Clark lab (nclark-lab). Interpretation tools: PHACT/PHACTboost and PHACE
(Kuru and Adebali).
