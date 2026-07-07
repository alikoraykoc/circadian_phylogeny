#!/usr/bin/env bash
# 04_supermatrix.sh
# Concatenated analysis. Three products:
#   (a) supermat: unconstrained partitioned ML tree = strongest in-house topology
#       estimate; cross-check against Upham, especially on transition branches.
#   (b) supermat_concat: single-model MFP on the concatenation to pick the uniform
#       RERconverge model (written to best_model.txt for step 04).
#   (c) supermat_fixed: Upham topology + concatenated branch lengths under that
#       model = the uniform-model FALLBACK length set.
set -euo pipefail
source "$(dirname "$0")/../config.sh"

# Stage trimmed alignments into one directory so IQ-TREE can auto-partition them.
stage="$RES/supermatrix/parts"
mkdir -p "$stage"
for g in $GENES; do cp "$RES/trim/$g.trim.$AAEXT" "$stage/$g.$AAEXT"; done

echo "== (a) unconstrained partitioned ML supermatrix =="
iqtree -p "$stage" -m MFP+MERGE -B 1000 \
        -pre "$RES/supermatrix/supermat" -T "$THREADS" -redo

echo "== (b) concatenated MFP for RERconverge uniform model =="
# Concatenate all trimmed alignments into a single FASTA, then run MFP once.
# This picks the single best amino-acid model for the full dataset.
concat="$RES/supermatrix/concat.faa"
python3 -c "
import os
from collections import defaultdict

aln_dir = '$stage'
ext = '.$AAEXT'

# Read each alignment, gap-pad missing species
all_species = set()
gene_seqs = {}  # {fname: {sp: seq}}
gene_lens = {}  # {fname: alignment_length}

for fname in sorted(os.listdir(aln_dir)):
    if not fname.endswith(ext): continue
    seqs = {}
    cur = None
    for line in open(os.path.join(aln_dir, fname)):
        line = line.strip()
        if line.startswith('>'):
            cur = line[1:]
            seqs[cur] = ''
        elif cur:
            seqs[cur] += line
    all_species.update(seqs.keys())
    alen = len(next(iter(seqs.values())))
    gene_seqs[fname] = seqs
    gene_lens[fname] = alen

# Concatenate with gap-padding for missing taxa
concat_seqs = defaultdict(str)
for fname in sorted(gene_seqs):
    seqs = gene_seqs[fname]
    alen = gene_lens[fname]
    for sp in all_species:
        concat_seqs[sp] += seqs.get(sp, '-' * alen)

with open('$concat', 'w') as out:
    for sp in sorted(concat_seqs):
        out.write(f'>{sp}\n{concat_seqs[sp]}\n')
"
iqtree -s "$concat" -m MFP -n 0 \
        -pre "$RES/supermatrix/supermat_concat" -T "$THREADS" -redo

# Extract best model and write it for step 04
best_model=$(grep "^Best-fit model" "$RES/supermatrix/supermat_concat.iqtree" \
             | sed 's/.*BIC: //')
echo "$best_model" > "$RER_MODEL_FILE"
echo "RERconverge uniform model: $best_model (written to $RER_MODEL_FILE)"

echo "== (c) constrained: Upham topology + concatenated branch lengths (fallback) =="
iqtree -s "$concat" -te "$SPTREE" -m "$best_model" \
        -pre "$RES/supermatrix/supermat_fixed" -T "$THREADS" -redo

echo
echo "COMPARE supermat.treefile with the Upham tree. Any conflict on a branch that"
echo "carries a diel transition is a priority to resolve before trusting calls there."
