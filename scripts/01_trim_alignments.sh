#!/usr/bin/env bash
# 01_trim_alignments.sh
# Trim each protein alignment with trimAl and retain a column-number mapping so
# that trimmed-column positions can later be translated back to the reference
# (human) residue numbering used when reporting candidate convergent sites.
set -euo pipefail
source "$(dirname "$0")/../config.sh"

for g in $GENES; do
  in="$ALN/${g}${AASUFFIX}.$AAEXT"
  out="$RES/trim/$g.trim.$AAEXT"
  colmap="$RES/trim/$g.colnumbering.txt"

  echo "== trimAl: $g =="
  # -automated1 chooses a heuristic suited to ML phylogenetics.
  # -colnumbering prints, for each RETAINED column, its index in the ORIGINAL
  #   alignment. Combined with the reference row this gives trimmed->human mapping.
  trimal -in "$in" -out "$out" -automated1 -colnumbering > "$colmap"
done

echo
echo "Trimmed alignments in $RES/trim/"
echo "NEXT: build_ref_coordinate_map.py turns each *.colnumbering.txt plus the"
echo "$REF_SPECIES row into a trimmed-column -> human-residue table (see 13 step)."
