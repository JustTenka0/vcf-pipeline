#!/bin/bash
set -euo pipefail

# ---------------------------------------------------------
# 1. GLOBAL PATHS SETUP
# ---------------------------------------------------------

# Modify them with your path and files

ROOT_DIR="$PWD"   # main work directory
REF_FASTA="${ROOT_DIR}/chr20.fa"
REF_INDEX="${ROOT_DIR}/chr20" 
BED_FILE="${ROOT_DIR}/chr20_ILMN_Exome_2.0_Plus_Panel.hg38_padded.bed"
VEP_CACHE="/data/vep_cache" 

# ID are the same for every trios
CHILD_ID="HG00427"
FATHER_ID="HG00428"
MOTHER_ID="HG00429"

cd "$ROOT_DIR" || { echo "Errore: directory non trovata"; exit 1; }


for CASE_DIR in trio_1 trio_2 trio_3 trio_4 trio_5; do
    if [ ! -d "$CASE_DIR" ]; then continue; fi
    
    CASE=$(basename "$CASE_DIR")
    echo "#####################################################"
    echo "# ANALYSIS FOR: $CASE"
    echo "#####################################################"
    
    cd "$CASE_DIR"

    # ---------------------------------------------------------
    # PART II & III: MAPPING E QC
    # ---------------------------------------------------------
    for ROLE in child father mother; do
        ID_VAR="${ROLE^^}_ID" # Trasform in CHILD_ID, FATHER_ID, etc.
        SAMPLE_ID="${!ID_VAR}"
        
        if [ ! -f "${ROLE}.bam" ]; then
            echo "[MAPPING] Generation of ${ROLE}.bam..."
            bowtie2 -x "$REF_INDEX" -1 "${SAMPLE_ID}.targets_R1.fq.gz" -2 "${SAMPLE_ID}.targets_R2.fq.gz" \
                --rg-id "$ROLE" --rg "SM:$ROLE" | samtools view -Sb | samtools sort -o "${ROLE}.bam"
            samtools index "${ROLE}.bam"
        else
            echo "[SKIP] ${ROLE}.bam already present."
        fi
    done
    # --- 3. QUALITY CONTROL ---
    if [ ! -d "qc_child" ] || [ ! -f "${CASE}_multiqc_report.html" ]; then
        echo "[QC] Qualimap and FastQC..."
        fastqc *.bam
        qualimap bamqc -bam child.bam --feature-file "$BED_FILE" -outdir qc_child --java-mem-size=4G > /dev/null 2>&1
        qualimap bamqc -bam father.bam --feature-file "$BED_FILE" -outdir qc_father --java-mem-size=4G > /dev/null 2>&1
        qualimap bamqc -bam mother.bam --feature-file "$BED_FILE" -outdir qc_mother --java-mem-size=4G > /dev/null 2>&1
        multiqc . -f -n "${CASE}_multiqc_report"
    else
        echo "[SKIP] Report QC already present."
    fi

    echo "child" > samples.txt
    echo "father" >> samples.txt
    echo "mother" >> samples.txt


    # ---------------------------------------------------------
    # PART IV: VARIANT CALLING
    # ---------------------------------------------------------
    echo "[STEP] Variant Calling."
    freebayes -f "$REF_FASTA" -t "$BED_FILE" -m 20 -C 5 -Q 10 --min-coverage 10 \
    child.bam father.bam mother.bam > "${CASE}_raw.vcf"

    bcftools sort "${CASE}_raw.vcf" -Oz -o "${CASE}_sorted.vcf.gz"
    bcftools index "${CASE}_sorted.vcf.gz"

    # ---------------------------------------------------------
    # PART V: FILTERING (Inheritance)
    # ---------------------------------------------------------


    # --- INHERITANCE FILTERING  ---
    echo "[STEP] Inheritance filtering..."
    case "$CASE" in
        "trio_2") FILT='GT[0]="RA" && GT[1]="RA" && GT[2]="RR"' ;; # AD paternal
        "trio_3") FILT='GT[0]="RA" && GT[1]="RR" && GT[2]="RR"' ;; # de novo AD
        *)        FILT='GT[0]="AA" && GT[1]="RA" && GT[2]="RA"' ;; # AR
    esac

    bcftools view -R "$BED_FILE" "${CASE}_sorted.vcf.gz" | \
    bcftools view -i "$FILT" | \
    bcftools filter -i 'QUAL>10' -Ov -o "${CASE}_candidates.vcf"


    # --- ANNOTAZIONE CON VEP ---
    echo "[STEP] VEP annotation..."
    
    vep -i "${CASE}_candidates.vcf" -o "${CASE}_annotated.vcf" \
        --vcf --cache --offline --dir_cache "$VEP_CACHE" \
        --assembly GRCh38 --use_given_ref --mane --pick_allele \
        --af --af_1kg --af_gnomade --max_af --sift b --polyphen b --no_fasta --force_overwrite

    # --- FILTRO FINALE (Rare & High Impact) ---
    echo "[STEP] Final report generation, IMPACT HIGH"
    filter_vep -i "${CASE}_annotated.vcf" -o "${CASE}_HIGH.vcf" \
        --filter "IMPACT is HIGH and (not MAX_AF or MAX_AF < 0.0001)" --force_overwrite
        
   echo "[STEP] Final report generation, IMPACT MODERATE"
    filter_vep -i "${CASE}_annotated.vcf" -o "${CASE}_MODERATE.vcf" \
        --filter "IMPACT is MODERATE and (not MAX_AF or MAX_AF < 0.0001)" --force_overwrite



    # ---------------------------------------------------------
    # PART VII: TRACKS PER UCSC
    # ---------------------------------------------------------
    echo "[STEP] BedGraph generation for UCSC..."
    bedtools genomecov -ibam child.bam -bg -trackline -trackopts "name=\"${CASE}_child\"" -max 100 > "${CASE}_child.bg"
    bedtools genomecov -ibam father.bam -bg -trackline -trackopts "name=\"${CASE}_father\"" -max 100 > "${CASE}_father.bg"
    bedtools genomecov -ibam mother.bam -bg -trackline -trackopts "name=\"${CASE}_mother\"" -max 100 > "${CASE}_mother.bg"
    cd "$ROOT_DIR"
    echo "### COMPLETE: $CASE ###"
done
