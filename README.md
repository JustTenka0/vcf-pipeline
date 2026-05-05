# Trio Genomic Analysis Pipeline (vcf-pipeline)

An automated Bash pipeline for processing Trio Genomic Data (Child, Father, Mother) paired-ends simulated data. This tool handles the workflow from raw FASTQ reads to filtered and annotated variant calls (VCF).

## ❗NOTES
1. Make sure to change the inheritance and directory names of the `inheritance.txt` file.
2. Files are heavy, download them only if you really need them. 
    * Target and Indexing files are reperible in `ref/` directory. 
    * Original file used are in `data/` directory.
3. For this case the files name are named in the same way for each trios, but are different.

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

The script requires specific file naming and directory structures to function correctly.

### 1. Naming Convention
Before starting, ensure your FASTQ files are named using the following pattern. Since the files are paired ends simulated data, there are 2 fq.gz for each individual (where `${CASE_ID}` is your unique identifier, e.g., `HG00427`):
* `HG00427.target_R1.fq.gz`
* `HG00427.target_R2.fq.gz`

### 2. Reference Files
You must have the following reference data ready:
* **Reference Genome**: A FASTA file (e.g., `chr20.fasta`) and its index (`.fai`).
* **Bowtie2 Index**: A set of `.bt2` files sharing the same prefix.
* **Target regions**: A `.bed` file defining your capture regions (e.g., Exome targets).

### 3. Folder Structure
The script will prompt for the **Reads Directory**. For best results, keep your raw data organized:
```text
working_directory/
├── pipeline.sh
└── trio_1/
    ├── HG00427.target_R1.fq.gz
    ├── HG00427.target_R2.fq.gz
    ├── HG00428.target_R1.fq.gz
    ├── HG00428.target_R2.fq.gz
    ├── HG00429.target_R1.fq.gz
    └── HG00429.target_R2.fq.gz
└── trio_2/
    ...
```
### 4. How to Use

1. **Prepare your environment**: Ensure all dependencies (Samtools, FreeBayes, etc.) are in your PATH.
2. **Execute the script**:
   ```bash
   bash pipeline.sh