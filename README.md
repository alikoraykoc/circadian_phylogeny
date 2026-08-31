# Circadian-gene convergence and mammalian diel activity

Testing whether convergent molecular evolution in circadian genes is associated
with independent transitions in diel activity (derived diurnality, plus reversals
to nocturnality) across mammals.

## The three trees, and why each exists

This is the single most important design point.

1. **Upham species tree (the scaffold).** Fixed topology for every downstream
   analysis. The diel scenario is mapped onto it once, so transition branches sit
   at the same nodes across all genes. Estimated from far more data than any one
   circadian gene, so it is the better genealogy on nearly every branch.
2. **Per-gene fixed-topology trees (the working input).** Same Upham topology,
   with gene-specific branch lengths in substitutions/site. This is what PCOC,
   TDG09, and RERconverge consume. Not an independent ML search per gene.
3. **Unconstrained per-gene trees + supermatrix tree (support tools only).**
   Used to police the analysis, never to scaffold it:
   - gene trees -> orthology QC and gene-concordance factors (hemiplasy screen);
   - supermatrix -> topology confirmation against Upham, and a uniform-model
     branch-length fallback.

Forcing the species topology removes topology-error artifacts; the concordance
step (06) then guards against the remaining hemiplasy artifact by flagging diel
transitions that fall on low-concordance branches.

## Analysis tracks (run in parallel, then reconverge)

- **Convergence:** PCOC + TDG09 (+ substitution counting) -> consensus sites.
- **Rate signal:** RERconverge -> phenotype-associated genes.
- **Selection (codon):** Contrast-FEL + RELAX. Needs codon alignments (pal2nal),
  so this track is gated on having CDS. Convergence and rate tracks do not.

Interpretation (ASR, AlphaFold structural mapping, PHACT/PHACE) then turns the
cross-method consensus into candidate functional sites.

## Run order

```
scripts/00_check_inputs.sh
scripts/01_trim_alignments.sh
scripts/02_codon_alignments.sh        # optional, selection track only
scripts/03_genetrees_qc.sh
scripts/04_supermatrix.sh              # picks RER uniform model
scripts/05_branchlengths_fixed.sh      # builds PCOC (MFP) + RER (uniform) trees
scripts/06_concordance.sh
Rscript scripts/07_scenario.R
python  scripts/08_flag_discordant_branches.py
scripts/09_pcoc.sh
scripts/10_tdg09.sh
Rscript scripts/11_rerconverge.R
scripts/12_selection_hyphy.sh          # optional, selection track only
python  scripts/13_consensus.py
```

## Setup

```
conda env create -f env/environment.yml && conda activate circadian-conv
docker pull carinerey/pcoc
# in R: remotes::install_github("nclark-lab/RERconverge")
# download tdg09.jar into tools/
```

Fill in the data/ directory (species_tree.nwk, diel_activity.csv, alignments/,
and optionally cds/), then edit config.sh (tip labels, gene list, reference).
Every script has `EDIT` / `TODO` markers where project-specific choices go.
See TODO.md for the working checklist.

## Checks

Run these after touching anything in the tree, scenario or alignment path, and
before sharing results.

```
python scripts/test_lib_checks.py          # unit tests for the interface guards
python scripts/check_reported_numbers.py   # write-ups still match the data
bash   scripts/spikein_quick.sh            # ~5 min end-to-end positive control
```

`check_reported_numbers.py` exists because three separate documents were found
quoting numbers that had been superseded by a later run: the parsimony counts
(213 against the current 201), the partial-convergence k-curve, and a stale set
of figures. It recomputes each headline number from the data, fails if a
write-up no longer states it, and flags any known-superseded number that appears
without a marker saying so. Superseded results live under
`shared_results/superseded/`, which has a README giving the provenance of each.
