# config.sh
# Shared variables for the circadian-gene / diel-activity convergence pipeline.
# Source this at the top of every shell script: source ../config.sh

# ---- paths (edit to your absolute paths if not running from scripts/) ----
export PROJ="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export DATA="$PROJ/data"
export ALN="$DATA/alignments"        # per-gene protein alignments, one file per gene
export CDS="$DATA/cds"               # per-gene CDS (nucleotide), only for the selection track
export RES="$PROJ/results"
export SPTREE="$DATA/species_tree.nwk"   # Upham tree pruned to your taxa
export DIEL="$DATA/diel_activity.csv"    # columns: species,activity

# ---- reference for coordinate mapping ----
export REF_SPECIES="Homo_sapiens"    # EDIT: exact tip label used in your alignments/tree

# ---- gene list (18 genes) ----
# EDIT: file basenames must match, e.g. $ALN/CLOCK.faa
export GENES="CLOCK NPAS2 ARNTL PER1 PER2 PER3 CRY1 CRY2 NR1D1 NR1D2 RORA RORB RORC CSNK1D CSNK1E FBXL3 BHLHE40 BHLHE41"

# ---- alignment file extensions ----
export AAEXT="faa"          # protein alignment extension
export CDSEXT="fna"         # CDS extension

# ---- compute ----
export THREADS="AUTO"       # or an integer for IQ-TREE / HyPhy

# ---- models ----
# One shared model across all fixed-topology gene trees keeps RERconverge's
# cross-gene rate comparison clean. LG+G4 is the defensible default here.
export FIXED_MODEL="LG+G4"

# ---- helper ----
mkdir -p "$RES"/{trim,codon,genetrees_qc,branchlengths,supermatrix,concordance,scenario,pcoc,tdg09,rerconverge,selection,consensus,figures}
