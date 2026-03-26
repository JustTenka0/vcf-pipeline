# Trio Genomic Analysis Pipeline (vcf-pipeline)

An automated Bash pipeline for processing Trio Genomic Data (Child, Father, Mother). This tool handles the workflow from raw FASTQ reads to filtered and annotated variant calls (VCF).

---

## 📋 Prerequisites

Ensure the following tools are installed and accessible in your `$PATH`:
* **Alignment**: `bowtie2`, `samtools`
* **Quality Control**: `fastqc`, `qualimap`, `multiqc`
* **Variant Calling**: `freebayes`, `bcftools`, `bgzip`
* **Annotation**: `vep` (Ensembl Variant Effect Predictor) & `filter_vep`
* **Coverage**: `bedtools`

---

## 🛠 Setup and Input Requirements

The script is interactive, but it requires specific file naming and directory structures to function correctly.

### 1. Naming Convention
Before starting, ensure your FASTQ files are named using the following pattern (where `${CASE_ID}` is your unique identifier, e.g., `FAM01`):
* `${CASE_ID}_child.fq.gz`
* `${CASE_ID}_father.fq.gz`
* `${CASE_ID}_mother.fq.gz`

### 2. Reference Files
You must have the following reference data ready:
* **Reference Genome**: A FASTA file (e.g., `universe.fasta`) and its index (`.fai`).
* **Bowtie2 Index**: A set of `.bt2` files sharing the same prefix.
* **Target regions**: A `.bed` file defining your capture regions (e.g., Exome targets).
* **VEP Cache**: A local directory containing the Ensembl cache.

### 3. Folder Structure
The script will prompt for the **Reads Directory**. For best results, keep your raw data organized:
```text
project_folder/
├── pipeline.sh
└── raw_data/
    ├── FAM01_child.fq.gz
    ├── FAM01_father.fq.gz
    └── FAM01_mother.fq.gz

### 4. How to Use

1. **Prepare your environment**: Ensure all dependencies (Samtools, FreeBayes, etc.) are in your PATH.
2. **Execute the script**:
   ```bash
   bash pipeline.sh