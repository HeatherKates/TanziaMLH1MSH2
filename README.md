# RNAseq, ATACseq, and CUT&RUN Analysis for MSH2 and MLH1 in Tumor Metastasis

This repository contains scripts and files related to the analysis of RNAseq, ATACseq, and CUT&RUN data to investigate the contrasting roles of MSH2 and MLH1 in tumor metastasis.

## Overview

The code provided here is fully reproducible. Except for the original `fastq.gz` files, all other files and results can be recreated using the scripts included in this repository. Due to size constraints and data ownership, large files and non public files are omitted. For access to these files, please see the `.gitignore` and email [hkates@ufl.edu](mailto:hkates@ufl.edu).

## Directory Structure

- **RNAseq**
  * **scripts/**: Contains all the scripts necessary to reproduce the analysis.  
    * script files are named beginning with a sequential order that they can be run in for full reproducibility. (A comes before B).
  * **results/**: Directory for storing analysis results.

- **ATACseq**
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

Warnings and issues encountered during the pipeline run are logged in:
  - Example: `QC_warnings.log`
  
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
