#!/usr/bin/env bash
# 05_supermatrix.sh
# Concatenated, partitioned analysis. Two products, both SUPPORT tools:
#   (a) supermat: unconstrained partitioned ML tree = strongest in-house topology
#       estimate; cross-check against Upham, especially on transition branches.
#   (b) supermat_fixed: Upham topology with branch lengths from the full
#       concatenation = the uniform-model FALLBACK length set.
set -euo pipefail
source "$(dirname "$0")/../config.sh"

# Stage trimmed alignments into one directory so IQ-TREE can auto-partition them.
stage="$RES/supermatrix/parts"
mkdir -p "$stage"
for g in $GENES; do cp "$RES/trim/$g.trim.$AAEXT" "$stage/$g.$AAEXT"; done

echo "== (a) unconstrained partitioned ML supermatrix =="
iqtree2 -p "$stage" -m MFP+MERGE -B 1000 \
        -pre "$RES/supermatrix/supermat" -T "$THREADS" -redo

echo "== (b) constrained: Upham topology + concatenated branch lengths (fallback) =="
iqtree2 -p "$stage" -te "$SPTREE" -m "$FIXED_MODEL" \
        -pre "$RES/supermatrix/supermat_fixed" -T "$THREADS" -redo

echo
echo "COMPARE supermat.treefile with the Upham tree. Any conflict on a branch that"
echo "carries a diel transition is a priority to resolve before trusting calls there."
