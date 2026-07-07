#!/usr/bin/env bash
# 03_genetrees_qc.sh
# Unconstrained ML tree per gene. TWO purposes only:
#   (a) QC: a taxon placed wildly wrong flags paralog contamination / bad alignment.
#   (b) input loci trees for the gene-concordance-factor step (06).
# These trees are NOT the analysis backbone and are NOT fed to PCOC/TDG09.
set -euo pipefail
source "$(dirname "$0")/../config.sh"

: > "$RES/genetrees_qc/all_gene_trees.treefile"   # start empty collection

for g in $GENES; do
  echo "== unconstrained ML: $g =="
  iqtree -s "$RES/trim/$g.trim.$AAEXT" \
          -m MFP -B 1000 \
          -pre "$RES/genetrees_qc/$g" -T "$THREADS" -redo
  cat "$RES/genetrees_qc/$g.treefile" >> "$RES/genetrees_qc/all_gene_trees.treefile"
done

echo
echo "QC trees + all_gene_trees.treefile ready."
echo "MANUAL CHECK: open each *.treefile and compare against the species tree."
echo "Decide per gene: keep / re-align / drop (record decisions in TODO.md)."
