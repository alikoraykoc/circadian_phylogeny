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
# EDIT: file basenames must match, e.g. $ALN/CLOCK_aligned.fa
export GENES="CLOCK NPAS2 ARNTL PER1 PER2 PER3 CRY1 CRY2 NR1D1 NR1D2 RORA RORB RORC CSNK1D CSNK1E FBXL3 BHLHE40 BHLHE41"

# ---- alignment file extensions ----
export AASUFFIX="_aligned"  # suffix between gene name and extension
export AAEXT="fa"           # protein alignment extension
export CDSEXT="fna"         # CDS extension

# ---- compute ----
export THREADS="AUTO"       # or an integer for IQ-TREE / HyPhy

# ---- models (two tracks) ----
# PCOC / TDG09 trees: per-gene MFP. These are independent per-gene analyses
# with no cross-gene comparison, so each gene gets its best-fit branch lengths.
export PCOC_MODEL="MFP"

# RERconverge trees: one uniform model across all genes, data-driven from the
# supermatrix (step 05). RER normalizes scale but NOT model-driven shape
# differences (long-vs-short branch ratios from different rate-heterogeneity
# choices). A shared model prevents that confound.
# Step 05 writes the winner to $RES/supermatrix/best_model.txt; step 04 reads it.
export RER_MODEL_FILE="$RES/supermatrix/best_model.txt"
export RER_MODEL_FALLBACK="Q.MAMMAL+R4"

# ---- helper ----
mkdir -p "$RES"/{trim,codon,genetrees_qc,branchlengths,supermatrix,concordance,scenario,pcoc,tdg09,rerconverge,selection,consensus,figures}
