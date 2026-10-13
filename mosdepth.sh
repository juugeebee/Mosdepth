#!/bin/bash
# =============================================================================
# run_mosdepth_coverage_v2.sh
# Calcul de couverture sur regions BED avec mosdepth, matrice + normalisation
# =============================================================================
set -uo pipefail

# =============================================================================
# PARAMETRES : a adapter
# =============================================================================
BAM_DIR='/media/jbogoin/Data3/Il-str-mRNA-D_NOVASEQ-Run6/BAM'
BED="/home/jbogoin/Downloads/Liste_genes_RNAseq_v2_final_06102026_exons.bed"
OUTDIR="${BAM_DIR}/mosdepth_results_exons"
FINAL="${BAM_DIR}/coverage_matrix_exons.tsv"
NORMALIZE_PY="$HOME/SCRIPTS/Mosdepth/normalize_coverage.py"
JOBS=18
LOG="${BAM_DIR}/mosdepth_pipeline.log"

# =============================================================================
# INIT
# =============================================================================
mkdir -p "$OUTDIR"
: > "$LOG"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

# =============================================================================
# ENVIRONNEMENT CONDA
# =============================================================================
# set +u : les scripts d'activation conda r�f�rencent des variables non d�finies
set +u
eval "$(conda shell.bash hook)" || { log "ERREUR : conda introuvable."; exit 1; }
conda activate mosdepth || { log "ERREUR : env conda 'mosdepth' introuvable."; exit 1; }
set -u
log "Env conda actif : ${CONDA_DEFAULT_ENV:-aucun}"

# =============================================================================
# VERIFICATION OUTILS
# =============================================================================
log "=== Verification de l'environnement ==="
for tool in mosdepth samtools parallel python3; do
    if ! command -v "$tool" &>/dev/null; then
        log "ERREUR : '$tool' introuvable."; exit 1
    fi
    log "  OK : $tool"
done

# =============================================================================
# VERIFICATION BED
# =============================================================================
log "=== Verification BED ==="
if [ ! -f "$BED" ]; then
    log "ERREUR : BED introuvable : $BED"; exit 1
fi
NB_REGIONS=$(grep -Evc '^(#|track|browser|$)' "$BED" || true)
log "  BED : $BED ($NB_REGIONS regions)"

# =============================================================================
# LISTE DES BAM
# =============================================================================
log "=== Inventaire des BAM ==="
BAM_LIST="${BAM_DIR}/all_bam.txt"
find "$BAM_DIR" -maxdepth 1 -name '*.bam' | sort > "$BAM_LIST"
NB_BAM=$(wc -l < "$BAM_LIST")
log "  BAM trouves : $NB_BAM"
if [ "$NB_BAM" -eq 0 ]; then
    log "ERREUR : aucun BAM dans $BAM_DIR."; exit 1
fi

# =============================================================================
# ETAPE 1 : INTEGRITE
# =============================================================================
log "=== Etape 1 : Integrite des BAM ==="

CORRUPT_LIST="${BAM_DIR}/corrupt_bam.txt"
VALID_LIST="${BAM_DIR}/valid_bam.txt"
: > "$CORRUPT_LIST"
: > "$VALID_LIST"

while IFS= read -r bam; do
    if samtools quickcheck "$bam" 2>/dev/null; then
        echo "$bam" >> "$VALID_LIST"
    else
        log "  CORRUPT : $bam"
        echo "$bam" >> "$CORRUPT_LIST"
    fi
done < "$BAM_LIST"

NB_CORRUPT=$(wc -l < "$CORRUPT_LIST")
NB_VALID=$(wc -l < "$VALID_LIST")
log "  Valides : $NB_VALID / $NB_BAM  |  Corrompus exclus : $NB_CORRUPT"

if [ "$NB_VALID" -eq 0 ]; then
    log "ERREUR : aucun BAM valide."; exit 1
fi

# =============================================================================
# ETAPE 2 : INDEXATION
# =============================================================================
log "=== Etape 2 : Indexation ==="
TO_INDEX="${BAM_DIR}/to_index.txt"
: > "$TO_INDEX"
while IFS= read -r bam; do
    [ -f "${bam}.bai" ] || [ -f "${bam%.bam}.bai" ] || echo "$bam" >> "$TO_INDEX"
done < "$VALID_LIST"

NB_TO_INDEX=$(wc -l < "$TO_INDEX")
if [ "$NB_TO_INDEX" -gt 0 ]; then
    log "  Indexation de $NB_TO_INDEX BAM..."
    parallel -j "$JOBS" 'samtools index {}' :::: "$TO_INDEX"
    log "  Indexation terminee."
else
    log "  Tous les BAM sont deja indexes."
fi

# =============================================================================
# ETAPE 3 : MOSDEPTH
# =============================================================================
log "=== Etape 3 : Calcul de couverture ==="

TODO_LIST="${BAM_DIR}/todo_bam.txt"
: > "$TODO_LIST"

while IFS= read -r bam; do
    sample=$(basename "$bam" .bam)
    res="${OUTDIR}/${sample}.regions.bed.gz"
    if [ -f "$res" ] && gzip -t "$res" 2>/dev/null; then
        log "  Deja fait (ignore) : $sample"
    else
        [ -f "$res" ] && log "  Resultat corrompu, recalcul : $sample"
        rm -f "${OUTDIR}/${sample}".*
        echo "$bam" >> "$TODO_LIST"
    fi
done < "$VALID_LIST"

NB_TODO=$(wc -l < "$TODO_LIST")
log "  BAM a traiter : $NB_TODO"

run_mosdepth() {
    local bam="$1"
    local sample
    sample=$(basename "$bam" .bam)
    if mosdepth --by "$BED" "${OUTDIR}/${sample}" "$bam"; then
        echo "DONE: $sample"
    else
        echo "FAILED: $sample"
    fi
}
export -f run_mosdepth
export BED OUTDIR

if [ "$NB_TODO" -gt 0 ]; then
    log "  Lancement de mosdepth ($JOBS jobs paralleles)..."
    parallel -j "$JOBS" run_mosdepth :::: "$TODO_LIST" | tee -a "$LOG"
    log "  mosdepth termine."
fi

# V�rification
FAILED=0
while IFS= read -r bam; do
    sample=$(basename "$bam" .bam)
    if [ ! -f "${OUTDIR}/${sample}.regions.bed.gz" ]; then
        log "  MANQUANT : $sample"
        FAILED=$((FAILED + 1))
    fi
done < "$VALID_LIST"

if [ "$FAILED" -gt 0 ]; then
    log "  ATTENTION : $FAILED resultats manquants."
else
    log "  Tous les resultats sont presents."
fi

# =============================================================================
# ETAPE 4 : MERGE
# =============================================================================
log "=== Etape 4 : Construction de la matrice ==="

MERGE_PY="$HOME/SCRIPTS/Mosdepth/merge_coverage.py"
if [ ! -f "$MERGE_PY" ]; then
    log "ERREUR : $MERGE_PY introuvable."; exit 1
fi

# merge_coverage.py : couverture = derniere colonne, nom de gene = colonne 4 du BED,
# noms d'�chantillons raccourcis. Sortie : chr start end gene <samples>
if ! python3 "$MERGE_PY" "$OUTDIR" "$BED" "$FINAL" 2>&1 | tee -a "$LOG"; then
    log "ERREUR : Echec de merge_coverage.py."; exit 1
fi

rm -f "$TODO_LIST" "$TO_INDEX" "$BAM_LIST"

# =============================================================================
# ETAPE 5 : NORMALISATION
# =============================================================================
log "=== Etape 5 : Normalisation ==="
if [ -f "$NORMALIZE_PY" ]; then
    if python3 "$NORMALIZE_PY" "$FINAL"; then
        log "  Normalisation terminee."
    else
        log "  ATTENTION : Echec de la normalisation."
    fi
else
    log "  ATTENTION : $NORMALIZE_PY introuvable, normalisation ignoree."
fi

# =============================================================================
# RESUME
# =============================================================================
NB_ROWS=$(( $(wc -l < "$FINAL") - 1 ))
NB_COLS=$(( $(head -1 "$FINAL" | tr '\t' '\n' | wc -l) - 4 ))

log "=== Pipeline termine ==="
log "  Matrice : $FINAL"
log "  Regions : $NB_ROWS  |  Samples : $NB_COLS"
log "  BAM corrompus exclus : $NB_CORRUPT  (liste : $CORRUPT_LIST)"
log "  Log : $LOG"