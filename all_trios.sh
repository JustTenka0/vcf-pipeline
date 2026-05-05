#!/bin/bash

set -euo pipefail


# =============================================================================
# CONFIGURATION — edit this section to match your environment
# =============================================================================

ROOT_DIR="$PWD"                     # run from project root, or set absolute path
INHERITANCE_FILE="${ROOT_DIR}/inheritance.txt"
REF_FASTA="${ROOT_DIR}/chr20.fa"
REF_INDEX="${ROOT_DIR}/chr20"       # bowtie2 index prefix
BED_FILE="${ROOT_DIR}/$(ls "${ROOT_DIR}"/*.padded.bed 2>/dev/null | head -1 | xargs basename)"
VEP_CACHE="/data/vep_cache"
VEP_ASSEMBLY="GRCh38"

# Variant calling thresholds
MIN_MAP_QUAL=20
MIN_ALT_COUNT=5
MIN_BASE_QUAL=10
MIN_COVERAGE=10
QUAL_FILTER=10

# Final variant filter: MAX_AF threshold
MAX_AF=0.0001

# VEP impact levels to retain (HIGH only for strict, add MODERATE for broader search)
VEP_IMPACT_FILTER="IMPACT is HIGH or IMPACT is MODERATE"


# =============================================================================
# PRE-FLIGHT CHECKS
# =============================================================================

echo "============================================="
echo " TRIO EXOME PIPELINE"
echo " Root: ${ROOT_DIR}"
echo "============================================="
echo ""

# Check required shared files
echo "Checking reference files..."
for F in "${REF_FASTA}" "${REF_FASTA}.fai" "${BED_FILE}" "${INHERITANCE_FILE}"; do
    [[ -f "$F" ]] \
        && echo "  ✓ $(basename "$F")" \
        || { echo "  ERROR: file not found: $F"; exit 1; }
done
ls "${REF_INDEX}".*.bt2 &>/dev/null \
    && echo "  ✓ bowtie2 index: $(basename "${REF_INDEX}")" \
    || { echo "  ERROR: bowtie2 index not found: ${REF_INDEX}*.bt2"; exit 1; }

# Detect trio folders automatically (any directory named trio_*)
TRIO_DIRS=( "${ROOT_DIR}"/trio_*/ )
if [[ ${#TRIO_DIRS[@]} -eq 0 ]]; then
    echo "  ERROR: no trio_* directories found in ${ROOT_DIR}"
    exit 1
fi

echo ""
echo "Trios detected: $(printf '%s ' "${TRIO_DIRS[@]}" | xargs -n1 basename | tr '\n' ' ')"
echo ""
echo "Press ENTER to start, Ctrl+C to abort."
read -r


# =============================================================================
# HELPER: read inheritance mode from inheritance.txt
# Returns the bcftools GT filter expression for the given trio folder name
# GT[0]=child  GT[1]=father  GT[2]=mother
# =============================================================================

get_inheritance_filter() {
    local TRIO_NAME="$1"
    local MODE AFFECTED FILT

    # Skip comment lines, match first column to trio name
    MODE=$(grep -v "^#" "${INHERITANCE_FILE}" \
           | awk -v t="${TRIO_NAME}" '$1==t {print $2}' | head -1)
    AFFECTED=$(grep -v "^#" "${INHERITANCE_FILE}" \
               | awk -v t="${TRIO_NAME}" '$1==t {print $3}' | head -1)

    if [[ -z "${MODE}" ]]; then
        echo "ERROR: '${TRIO_NAME}' not found in ${INHERITANCE_FILE}" >&2
        exit 1
    fi

    case "${MODE}" in
        AR)
            # Child homozygous alt (AA), both parents heterozygous (RA)
            FILT='GT[0]="AA" && GT[1]="RA" && GT[2]="RA"'
            ;;
        AD_inherited)
            # Child het (RA), affected parent het (RA), other parent ref (RR)
            if [[ "${AFFECTED}" == "father" ]]; then
                FILT='GT[0]="RA" && GT[1]="RA" && GT[2]="RR"'
            elif [[ "${AFFECTED}" == "mother" ]]; then
                FILT='GT[0]="RA" && GT[1]="RR" && GT[2]="RA"'
            else
                # Affected parent unknown: accept either
                FILT='GT[0]="RA" && (GT[1]="RA" || GT[2]="RA")'
            fi
            ;;
        AD_denovo)
            # Child het (RA), both parents homozygous ref (RR)
            FILT='GT[0]="RA" && GT[1]="RR" && GT[2]="RR"'
            ;;
        *)
            echo "ERROR: unknown mode '${MODE}' for ${TRIO_NAME}" >&2
            exit 1
            ;;
    esac

    echo "${FILT}"
}


# =============================================================================
# HELPER: identify child/father/mother IDs from FASTQ filenames
# Expects files named: SAMPLEID.targets_R1.fq.gz
# Trios must contain exactly 3 samples (6 files total)
# =============================================================================

get_sample_ids() {
    local TRIO_DIR="$1"

    # Collect all unique sample IDs from R1 filenames
    local IDS
    IDS=$(ls "${TRIO_DIR}"*.targets_R1.fq.gz 2>/dev/null \
          | xargs -n1 basename \
          | sed 's/\.targets_R1\.fq\.gz//' \
          | sort)

    local COUNT
    COUNT=$(echo "${IDS}" | wc -l)

    if [[ "${COUNT}" -ne 3 ]]; then
        echo "ERROR: expected 3 samples in ${TRIO_DIR}, found ${COUNT}" >&2
        exit 1
    fi

    echo "${IDS}"
}


# =============================================================================
# MAIN LOOP — one iteration per trio_* directory
# =============================================================================

for TRIO_DIR in "${TRIO_DIRS[@]}"; do

    TRIO_NAME=$(basename "${TRIO_DIR}")

    echo ""
    echo "============================================="
    echo " Processing: ${TRIO_NAME}"
    echo "============================================="

    # --- Get inheritance filter ---
    FILT=$(get_inheritance_filter "${TRIO_NAME}")
    echo "  Inheritance filter : ${FILT}"

    # --- Get sample IDs from FASTQ filenames ---
    # IDs are sorted alphabetically; trios.txt column order (child/father/mother)
    # must be reflected in how you name or sort your files.
    # If your project uses trios.txt, replace this block with a grep on that file.
    SAMPLE_IDS=$(get_sample_ids "${TRIO_DIR}")
    CHILD_ID=$(echo  "${SAMPLE_IDS}" | sed -n '1p')
    FATHER_ID=$(echo "${SAMPLE_IDS}" | sed -n '2p')
    MOTHER_ID=$(echo "${SAMPLE_IDS}" | sed -n '3p')

    echo "  Child  : ${CHILD_ID}"
    echo "  Father : ${FATHER_ID}"
    echo "  Mother : ${MOTHER_ID}"

    # --- Check all FASTQ files ---
    for ID in "${CHILD_ID}" "${FATHER_ID}" "${MOTHER_ID}"; do
        for R in R1 R2; do
            FQ="${TRIO_DIR}${ID}.targets_${R}.fq.gz"
            [[ -f "${FQ}" ]] \
                && echo "  ✓ $(basename "${FQ}")" \
                || { echo "  ERROR: ${FQ} not found"; exit 1; }
        done
    done

    cd "${TRIO_DIR}"

    # samples.txt — order must match BAM order: child[0] father[1] mother[2]
    printf "child\nfather\nmother\n" > samples.txt


    # =========================================================================
    # PART II — ALIGNMENT (paired-end)
    # Skip if BAM already exists (allows re-running from any step)
    # =========================================================================

    for ROLE in child father mother; do
        case "${ROLE}" in
            child)  ID="${CHILD_ID}"  ;;
            father) ID="${FATHER_ID}" ;;
            mother) ID="${MOTHER_ID}" ;;
        esac

        if [[ -f "${ROLE}.bam" ]]; then
            echo "  [SKIP] ${ROLE}.bam already exists"
        else
            echo "  [$(date +%T)] Aligning ${ROLE} (${ID})..."
            bowtie2 \
                -x "${REF_INDEX}" \
                -1 "${TRIO_DIR}${ID}.targets_R1.fq.gz" \
                -2 "${TRIO_DIR}${ID}.targets_R2.fq.gz" \
                --rg-id "${ROLE}" --rg "SM:${ROLE}" \
            | samtools view -Sb \
            | samtools sort -o "${ROLE}.bam"
            samtools index "${ROLE}.bam"
            echo "  [$(date +%T)] Done: ${ROLE}.bam"
        fi
    done


    # =========================================================================
    # PART III — QC
    # Skip if multiqc report already exists
    # =========================================================================

    if [[ -f "${TRIO_NAME}_multiqc_report.html" ]]; then
        echo "  [SKIP] QC report already exists"
    else
        echo "  [$(date +%T)] Running QC..."
        fastqc child.bam father.bam mother.bam

        qualimap bamqc -bam child.bam  --feature-file "${BED_FILE}" \
            -outdir qc_child  --java-mem-size=4G > /dev/null 2>&1
        qualimap bamqc -bam father.bam --feature-file "${BED_FILE}" \
            -outdir qc_father --java-mem-size=4G > /dev/null 2>&1
        qualimap bamqc -bam mother.bam --feature-file "${BED_FILE}" \
            -outdir qc_mother --java-mem-size=4G > /dev/null 2>&1

        multiqc . --force -n "${TRIO_NAME}_multiqc_report"
    fi


    # =========================================================================
    # PART IV — VARIANT CALLING
    # BAM order fixed: child[0]  father[1]  mother[2]
    # =========================================================================

    if [[ -f "${TRIO_NAME}_sorted.vcf.gz" ]]; then
        echo "  [SKIP] VCF already exists"
    else
        echo "  [$(date +%T)] Calling variants..."
        freebayes \
            -f "${REF_FASTA}" \
            -t "${BED_FILE}" \
            -m "${MIN_MAP_QUAL}" \
            -C "${MIN_ALT_COUNT}" \
            -Q "${MIN_BASE_QUAL}" \
            --min-coverage "${MIN_COVERAGE}" \
            child.bam father.bam mother.bam \
            > "${TRIO_NAME}_raw.vcf"

        bcftools sort "${TRIO_NAME}_raw.vcf" -Oz -o "${TRIO_NAME}_sorted.vcf.gz"
        bcftools index "${TRIO_NAME}_sorted.vcf.gz"
    fi


    # =========================================================================
    # PART V — FILTERING AND ANNOTATION
    # =========================================================================

    echo "  [$(date +%T)] Filtering by inheritance pattern..."

    bcftools view -R "${BED_FILE}" "${TRIO_NAME}_sorted.vcf.gz" \
        | bcftools view -S samples.txt \
        | bcftools view -i "(${FILT})" \
        | bcftools filter -i "QUAL>${QUAL_FILTER}" \
        -Ov -o "${TRIO_NAME}_candidates.vcf"

    NCAND=$(grep -v "^#" "${TRIO_NAME}_candidates.vcf" | wc -l)
    echo "  Candidate variants after inheritance filter: ${NCAND}"

    echo "  [$(date +%T)] Annotating with VEP..."

    vep \
        -i "${TRIO_NAME}_candidates.vcf" \
        -o "${TRIO_NAME}_annotated.vcf" \
        --vcf --cache --offline \
        --assembly "${VEP_ASSEMBLY}" \
        --dir_cache "${VEP_CACHE}" \
        --no_fasta \
        --use_given_ref --mane --pick_allele \
        --af --af_1kg --af_gnomade --max_af \
        --sift b --polyphen b \
        --force_overwrite

    echo "  [$(date +%T)] Applying rarity + impact filter..."

    filter_vep \
        -i "${TRIO_NAME}_annotated.vcf" \
        -o "${TRIO_NAME}_filtered_final.vcf" \
        --filter "(${VEP_IMPACT_FILTER}) and (not MAX_AF or MAX_AF < ${MAX_AF})" \
        --force_overwrite

    NFINAL=$(grep -v "^#" "${TRIO_NAME}_filtered_final.vcf" | wc -l)
    echo "  Final candidate variants: ${NFINAL}"


    # =========================================================================
    # PART VI — COVERAGE TRACKS FOR UCSC/IGV
    # =========================================================================

    echo "  [$(date +%T)] Generating coverage tracks..."

    bedtools genomecov -ibam child.bam  -bg -trackline \
        -trackopts "name=\"${TRIO_NAME}_child\""  -max 100 > "${TRIO_NAME}_child.bg"
    bedtools genomecov -ibam father.bam -bg -trackline \
        -trackopts "name=\"${TRIO_NAME}_father\"" -max 100 > "${TRIO_NAME}_father.bg"
    bedtools genomecov -ibam mother.bam -bg -trackline \
        -trackopts "name=\"${TRIO_NAME}_mother\"" -max 100 > "${TRIO_NAME}_mother.bg"

    echo ""
    echo "  ✓ ${TRIO_NAME} COMPLETE"
    echo "    Candidates : ${TRIO_NAME}_filtered_final.vcf (${NFINAL} variants)"
    echo "    Coverage   : ${TRIO_NAME}_child/father/mother.bg"
    echo "    QC         : ${TRIO_NAME}_multiqc_report.html"

    cd "${ROOT_DIR}"

done

echo ""
echo "============================================="
echo " ALL TRIOS COMPLETE"
echo "============================================="
