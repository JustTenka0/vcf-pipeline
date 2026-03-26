#!/bin/bash
# =============================================================================
# TRIO GENOMIC ANALYSIS PIPELINE
# Usage: bash pipeline.sh
# The script will prompt for all required paths before starting.
# =============================================================================

set -euo pipefail


# =============================================================================
# INTERACTIVE CONFIGURATION
# =============================================================================

echo "============================================="
echo " TRIO VARIANT CALLING PIPELINE — SETUP"
echo "============================================="
echo ""

# --- Case ID ---
echo "Enter case ID (or press ENTER to use current folder name: $(basename "$PWD")):"
read -r INPUT_CASE
CASE="${INPUT_CASE:-$(basename "$PWD")}"
echo "  → Case ID: ${CASE}"
echo ""

# --- Reads directory ---
echo "Enter the directory containing the FASTQ files"
echo "(or press ENTER to use current directory: $PWD):"
read -r INPUT_READS_DIR
READS_DIR="${INPUT_READS_DIR:-$PWD}"
echo "  → Reads directory: ${READS_DIR}"
echo ""

# Verify that expected FASTQ files actually exist before proceeding
FQ_CHILD="${READS_DIR}/${CASE}_child.fq.gz"
FQ_FATHER="${READS_DIR}/${CASE}_father.fq.gz"
FQ_MOTHER="${READS_DIR}/${CASE}_mother.fq.gz"

echo "Checking FASTQ files..."
for FQ in "${FQ_CHILD}" "${FQ_FATHER}" "${FQ_MOTHER}"; do
    if [[ ! -f "${FQ}" ]]; then
        echo "  ERROR: File not found: ${FQ}"
        echo "  Check the case ID and reads directory, then re-run."
        exit 1
    fi
    echo "  ✓ Found: $(basename "${FQ}")"
done
echo ""

# --- Reference genome ---
echo "Enter the directory containing the reference genome"
echo "(or press ENTER to use current directory: $PWD):"
read -r INPUT_REF_DIR
REF_DIR="${INPUT_REF_DIR:-$PWD}"

echo "Enter the bowtie2 index prefix (press ENTER for: uni):"
read -r INPUT_INDEX
BOWTIE2_INDEX="${REF_DIR}/${INPUT_INDEX:-uni}"

echo "Enter the reference FASTA filename (press ENTER for: universe.fasta):"
read -r INPUT_FASTA
FASTA="${REF_DIR}/${INPUT_FASTA:-universe.fasta}"

# --- Target regions ---
echo "Enter the full path to the BED targets file"
echo "(press ENTER for: ${REF_DIR}/targetsPad100.bed):"
read -r INPUT_TARGETS
TARGETS="${INPUT_TARGETS:-${REF_DIR}/targetsPad100.bed}"

# --- VEP cache ---
echo "Enter the VEP cache directory (press ENTER for: /data/vep_cache):"
read -r INPUT_VEP
VEP_CACHE_DIR="${INPUT_VEP:-/data/vep_cache}"


# --- Inheritance model ---
echo "Select inheritance model:"
echo "  1) Autosomal Recessive  — child AA, father RA, mother RA  [default]"
echo "  2) De novo              — child RA, father RR, mother RR"
echo "  3) Autosomal Dominant   — child RA, at least one parent RA"
echo "  4) X-linked Recessive   — child AA, father RR, mother RA"
echo "  5) Custom               — enter your own bcftools expression"
read -r INH_CHOICE

case "${INH_CHOICE}" in
    2) INHERITANCE_FILTER='GT[0]="RA" && GT[1]="RR" && GT[2]="RR"' ;;
    3) INHERITANCE_FILTER='GT[0]="RA" && (GT[1]="RA" || GT[2]="RA")' ;;
    4) INHERITANCE_FILTER='GT[0]="AA" && GT[1]="RR" && GT[2]="RA"' ;;
    5)
        echo "Enter custom bcftools genotype filter expression:"
        echo "(GT index: 0=child, 1=father, 2=mother)"
        read -r INHERITANCE_FILTER ;;
    *)  INHERITANCE_FILTER='GT[0]="AA" && GT[1]="RA" && GT[2]="RA"' ;;
esac
echo "  → Inheritance filter: ${INHERITANCE_FILTER}"
echo ""



echo "Checking reference files..."
for REF_FILE in "${FASTA}" "${FASTA}.fai" "${TARGETS}"; do
    if [[ ! -f "${REF_FILE}" ]]; then
        echo "  ERROR: File not found: ${REF_FILE}"
        exit 1
    fi
    echo "  ✓ Found: $(basename "${REF_FILE}")"
done

# Controlla almeno un file dell'indice bowtie2
if ! ls "${BOWTIE2_INDEX}".*.bt2 &>/dev/null; then
    echo "  ERROR: Bowtie2 index not found: ${BOWTIE2_INDEX}.*.bt2"
    exit 1
fi
echo "  ✓ Found bowtie2 index: $(basename "${BOWTIE2_INDEX}")"
echo ""

# --- Summary before running ---
echo "============================================="
echo " CONFIGURATION SUMMARY"
echo "============================================="
echo "  Case ID       : ${CASE}"
echo "  Reads dir     : ${READS_DIR}"
echo "  Reference     : ${FASTA}"
echo "  Bowtie2 index : ${BOWTIE2_INDEX}"
echo "  Targets BED   : ${TARGETS}"
echo "  VEP cache     : ${VEP_CACHE_DIR}"
echo "  Inheritance   : ${INHERITANCE_FILTER}"
echo "  Sample order  : child [0] → father [1] → mother [2]"
echo "============================================="
echo ""
echo "Press ENTER to start the pipeline, or Ctrl+C to abort."
read -r


# =============================================================================
# Fixed parameters (edit here if needed, not above)
# =============================================================================

VEP_ASSEMBLY="GRCh37"
MIN_MAP_QUAL=20
MIN_ALT_COUNT=5
MIN_BASE_QUAL=10
MIN_COVERAGE=10
MIN_QUAL=20
MAX_AF_THRESHOLD=0.0001


# =============================================================================
# PART I — DATA SETUP
# =============================================================================

echo "[$(date)] Setting up output directory: ${CASE}_analysis"

mkdir -p "${CASE}_analysis"
cd "${CASE}_analysis"

ln -sf "${BOWTIE2_INDEX}".* .
ln -sf "${FASTA}" .
ln -sf "${FASTA}.fai" .
ln -sf "${TARGETS}" .

# samples.txt: fixed order child → father → mother
# This must match the BAM order passed to freebayes (GT[0], GT[1], GT[2])
printf "child\nfather\nmother\n" > samples.txt


# =============================================================================
# PART II — ALIGNMENT
# =============================================================================

echo "[$(date)] Starting alignment..."

for MEMBER in child father mother; do
    FQ_VAR="FQ_${MEMBER^^}"
    FQ="${!FQ_VAR}"

    echo "  [$(date)] Aligning ${MEMBER}: ${FQ}"

    bowtie2 \
        -U "${FQ}" \
        -x "$(basename "${BOWTIE2_INDEX}")" \
        --rg-id "${MEMBER}" \
        --rg "SM:${MEMBER}" \
    | samtools view -Sb \
    | samtools sort -o "${MEMBER}.bam"

    samtools index "${MEMBER}.bam"
    echo "  [$(date)] Done: ${MEMBER}.bam"
done


# =============================================================================
# PART III — QC
# =============================================================================

echo "[$(date)] Running QC..."

fastqc *.bam

for MEMBER in child father mother; do
    qualimap bamqc \
        -bam "${MEMBER}.bam" \
        --feature-file "$(basename "${TARGETS}")" \
        -outdir "${MEMBER}_qualimap"
done

multiqc .


# =============================================================================
# PART IV — VARIANT CALLING
# Order: child.bam father.bam mother.bam → GT[0]=child GT[1]=father GT[2]=mother
# =============================================================================

echo "[$(date)] Calling variants..."

freebayes \
    -f "$(basename "${FASTA}")" \
    -m "${MIN_MAP_QUAL}" \
    -C "${MIN_ALT_COUNT}" \
    -Q "${MIN_BASE_QUAL}" \
    --min-coverage "${MIN_COVERAGE}" \
    child.bam father.bam mother.bam \
    > "${CASE}.vcf"

bgzip "${CASE}.vcf"
bcftools index "${CASE}.vcf.gz"


# =============================================================================
# PART V — FILTERING AND ANNOTATION
# =============================================================================

echo "[$(date)] Filtering by inheritance pattern and quality..."

bcftools view -R "$(basename "${TARGETS}")" "${CASE}.vcf.gz" \
    | bcftools view -S samples.txt \
    | bcftools view -i "(${INHERITANCE_FILTER})" \
    | bcftools filter -i "QUAL>${MIN_QUAL}" \
    -Ov -o "${CASE}.cand.vcf"

echo "[$(date)] Annotating with VEP..."

vep \
    -i "${CASE}.cand.vcf" \
    -o "${CASE}.vep_annotated.vcf" \
    --vcf --cache --offline \
    --assembly "${VEP_ASSEMBLY}" \
    --dir_cache "${VEP_CACHE_DIR}" \
    --use_given_ref --mane --pick_allele \
    --af --af_1kg --af_gnomade --max_af \
    --sift b --polyphen b

echo "[$(date)] Applying final rarity + impact filter..."

filter_vep \
    -i "${CASE}.vep_annotated.vcf" \
    -o "${CASE}.vep_filtered.vcf" \
    --filter "IMPACT is HIGH and (not MAX_AF or MAX_AF < ${MAX_AF_THRESHOLD})"


# =============================================================================
# PART VI — COVERAGE TRACKS
# =============================================================================

echo "[$(date)] Generating coverage tracks..."

for MEMBER in child father mother; do
    bedtools genomecov \
        -ibam "${MEMBER}.bam" \
        -bg -trackline \
        -trackopts "name=\"${MEMBER}\"" \
        -max 100 \
        > "${MEMBER}Cov.bg"
done


echo ""
echo "============================================="
echo " PIPELINE COMPLETE"
echo "  Final output : ${CASE}.vep_filtered.vcf"
echo "  Coverage     : childCov.bg  fatherCov.bg  motherCov.bg"
echo "============================================="
```


