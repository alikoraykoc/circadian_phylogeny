#!/usr/bin/env bash
# 02_codon_alignments.sh   (OPTIONAL: only needed for the selection track, step 12)
# Back-translate each trimmed protein alignment into a codon alignment using the
# matching CDS. HyPhy (Contrast-FEL, RELAX) needs codon-level input.
#
# Requirement: for every gene, $CDS/<gene>.fna must contain the SAME species as
# the protein alignment, with headers matching the alignment headers. The CDS
# sequences must translate to the (untrimmed) protein sequences.
set -euo pipefail
source "$(dirname "$0")/../config.sh"

have_cds=1
for g in $GENES; do
  [[ -f "$CDS/$g.$CDSEXT" ]] || { echo "No CDS for $g"; have_cds=0; }
done
if [[ $have_cds -eq 0 ]]; then
  echo "CDS incomplete. The selection track (Contrast-FEL / RELAX) is DEFERRED."
  echo "Convergence and rate-signal tracks do not need this step."
  exit 0
fi

for g in $GENES; do
  # NOTE: pal2nal maps a PROTEIN ALIGNMENT to unaligned CDS. Use the UNTRIMMED
  # protein alignment here so codon coordinates stay consistent, then trim the
  # codon alignment afterwards if desired, OR feed the trimmed protein alignment
  # if your pal2nal build supports gapped input cleanly. Check with -h.
  pal2nal.pl "$ALN/$g.$AAEXT" "$CDS/$g.$CDSEXT" -output fasta > "$RES/codon/$g.codon.fasta"
  echo "codon alignment: $g"
done

echo "Codon alignments in $RES/codon/"
