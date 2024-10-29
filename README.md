# RNAseq, ATACseq, and CUT&RUN Analysis for MSH2 and MLH1 in Tumor Metastasis

This repository contains scripts and files related to the analysis of RNAseq, ATACseq, and CUT&RUN data to investigate the contrasting roles of MSH2 and MLH1 in tumor metastasis.

## Overview

The code provided here is fully reproducible. Except for the original `fastq.gz` files, all other files and results can be recreated using the scripts included in this repository. Due to size constraints and data ownership, large files and non public files are omitted. For access to these files, please see the `.gitignore` and email [hkates@ufl.edu](mailto:hkates@ufl.edu).

## Directory Structure

# RNAseq Pipeline Summary

This directory contains the RNA-seq pipeline files and scripts to analyze the roles of MSH2 and MLH1 in tumor metastasis. Below is an overview of the main folders and the contents within each directory.

**ALL SCRIPTS SHOULD BE RUN FROM THE scripts/ DIRECTORY.**

**MAKE SURE THAT YOU CHANGE THE HIPERGATOR RESOURSE REQUESTS IN ANY .SBATCH FILES**

Nothing else needs to be changed

**Run *sbatch scripts by typing `sbatch FILENAME`**
**Run *R scripts by typing `module load R` and then `Rscript FILENAME`**

## Directory Structure

* **results/**: This directory contains subdirectories corresponding to the main stages of the RNA-seq analysis pipeline.
  * **1_fastqc/**: Quality control outputs from FASTQ files.
  * **2_fastp/**: Contains results from FASTQ trimming and preprocessing with Fastp.
  * **3_salmon/**: Outputs from transcript quantification using Salmon.
  * **4-6_deseq2/**: DESeq2 results for differential expression analysis, including MSH2 and MLH1 comparisons.
  * **7_GSEA/**: Results from Gene Set Enrichment Analysis (GSEA).
  * **8_GO/**: GO enrichment analysis outputs.

* **scripts/**: Contains all scripts used to perform each step of the pipeline, named sequentially to indicate the recommended order of execution.
  * **1a_fastqc.sbatch**: Script to perform quality control using FastQC.
  * **1b_summarize_fastqc.bash**: Summarizes FastQC results for review.
  * **2_fastp.sbatch**: Trims and preprocesses FASTQ files with Fastp.
  * **3a_salmon_human_index.sbatch**: Creates a human transcriptome index for Salmon.
  * **3b_salmon_human.sbatch**: Runs Salmon for transcript quantification.
  * **4a_deseq2.MLH1.R**: DESeq2 analysis for MLH1-related samples.
  * **4b_deseq2.MSH2.R**: DESeq2 analysis for MSH2-related samples.
  * **5_summarize_deseq2.R**: Summarizes DESeq2 results for both MSH2 and MLH1 comparisons.
  * **6_deseq2_visualizations.R**: Generates visualizations of DESeq2 differential expression results.
  * **7a_gsea.MLH1.R**: Gene Set Enrichment Analysis for MLH1.
  * **7b_gsea.MSH2.R**: Gene Set Enrichment Analysis for MSH2.
  * **8a_go_enrichment.MLH1.R**: GO enrichment analysis for MLH1-related genes.
  * **8b_go_enrichment.MSH2.R**: GO enrichment analysis for MSH2-related genes.
  * **logs/**: Directory for log files generated by each step of the pipeline.

Each script is designed to facilitate reproducibility and can be run sequentially to perform the complete RNA-seq analysis for the study.

# ATAC-seq Pipeline Summary

This pipeline processes ATAC-seq data starting from raw FASTQ files through trimming, alignment, peak calling, and differential binding analysis. Below is a summary of each step and the associated input/output files.

## Step 1: Trim Galore

Runs Trim Galore to remove adapter sequences from the raw FASTQ files.

**Input Files:**
  - Raw FASTQ files (`R1` and `R2`) for each lane:
    - Example: `231-MLH1KO-1_S1_L001_R1_001.fastq.gz`, `231-MLH1KO-1_S1_L001_R2_001.fastq.gz`

**Output Files:**
  - Trimmed FASTQ files for each lane:
    - Example: `MLH1KO-1_L001_val_1.fq.gz`, `MLH1KO-1_L001_val_2.fq.gz`
  - Trimming report files:
    - Example: `231-MLH1KO-1_S1_L001_R1_001.fastq.gz_trimming_report.txt`

## Step 2: FastQC (Pre- and Post-trimming)

Runs FastQC on both raw and trimmed FASTQ files to assess read quality.

**Input Files:**
  - Raw FASTQ files (`R1`, `R2`) for each lane
  - Trimmed FASTQ files (`val_1.fq.gz`, `val_2.fq.gz`)

**Output Files:**
  - FastQC HTML reports:
    - Example: `231-MLH1KO-1_S1_L001_R1_001_fastqc.html`
  - FastQC data files (`fastqc_data.txt` inside the `.zip`):
    - Example: `231-MLH1KO-1_S1_L001_R1_001_fastqc.zip`

## Step 3: Bowtie2 Alignment

Aligns trimmed reads to the reference genome using Bowtie2.

**Input Files:**
  - Trimmed FASTQ files from Trim Galore:
    - Example: `MLH1KO-1_L001_val_1.fq.gz`, `MLH1KO-1_L001_val_2.fq.gz`

**Output Files:**
  - SAM alignment files:
    - Example: `MLH1KO-1.sam`

## Step 4: Samtools Processing

Converts SAM to BAM, applies filtering, and sorts BAM files for downstream analysis.

**Input Files:**
  - SAM alignment files:
    - Example: `MLH1KO-1.sam`

**Output Files:**
  - Filtered and sorted BAM files:
    - Example: `MLH1KO-1.sorted.bam`
  - Flagstat statistics:
    - Example: `MLH1KO-1.bam.flagstat.log`

## Step 5: Peak Calling with Genrich

Calls ATAC-seq peaks using Genrich.

**Input Files:**
  - Filtered and sorted BAM files from Samtools:
    - Example: `MLH1KO-1.sorted.bam`

**Output Files:**
  - NarrowPeak files (peak regions):
    - Example: `MLH1KO-1.narrowPeak`

## Step 6: BigWig Generation

Creates BigWig files from BAM files for visualization in genome browsers.

**Input Files:**
  - Filtered and sorted BAM files from Samtools:
    - Example: `MLH1KO-1.sorted.bam`

**Output Files:**
  - BigWig files for visualization:
    - Example: `MLH1KO-1.bw`

## Step 7: Differential Binding Analysis with DiffBind

Performs differential binding analysis between groups (e.g., KO vs WT) using DiffBind.

**Input Files:**
  - NarrowPeak files for each sample:
    - Example: `MLH1KO-1.narrowPeak`
  - BAM files for each sample:
    - Example: `MLH1KO-1.sorted.bam`

**Output Files:**
  - CSV file with differential binding results:
    - Example: `MLH1_differential_binding_results.csv`
  - Filtered differential binding results for KO and R4 samples:
    - Example: `MLH1KO_filtered_differential_binding_results.csv`
    - Example: `MLH1R4_filtered_differential_binding_results.csv`
  - BED file with peaks data (optional):
    - Example: `MLH1_filtered_differential_binding_results.bed`

# Logging

Warnings and issues encountered during the pipeline run are logged in logs/:
  - Example: 
  
## Reproducibility

To ensure full reproducibility, follow the steps below:

1. **Download Original Data**: Obtain the original `fastq.gz` files or run scripts on hipergator where they will access data in /orange (persmission is for users in zhangw only)
2. **Run Scripts**: Use the scripts provided in the `scripts` directory to process the data and generate results.

## Contact

For access to large files or any other inquiries, please contact:

Heather Kates  
Email: [hkates@ufl.edu](mailto:hkates@ufl.edu)

## License

None

## Acknowledgements

This research was supported by [unknown]
