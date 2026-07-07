#!/usr/bin/env bash
# 06_concordance.sh
# Annotate every branch of the Upham tree with gene concordance factor (gCF) and
# site concordance factor (sCF). Low-gCF branches are where gene-tree/species-tree
# discordance concentrates, i.e. the branches most exposed to hemiplasy-driven
# false convergence. These values gate the scenario in step 08.
set -euo pipefail
source "$(dirname "$0")/../config.sh"

iqtree -t "$SPTREE" \
        --gcf "$RES/genetrees_qc/all_gene_trees.treefile" \
        -s "$RES/supermatrix/parts" --scf 100 \
        --prefix "$RES/concordance/concord" -T 4 -redo

echo
echo "Outputs:"
echo "  concord.cf.tree      Upham tree with gCF/sCF annotations"
echo "  concord.cf.stat      per-branch gCF/sCF table (feeds step 08)"
