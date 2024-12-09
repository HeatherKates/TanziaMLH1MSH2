#!/bin/bash

# Define output CSV file
output_csv="diffbind_samplesheet.csv"

# Write header to the output file
echo "SampleID,Tissue,Factor,Condition,Treatment,Replicate,Tissue,bamReads,Peaks,PeakCaller" > $output_csv

# Process each directory (KO and R4)
for dir in ../KO ../R4; do
  # Read samplesheet.csv line by line (skip header)
  tail -n +2 "$dir/samplesheet.csv" | while IFS=',' read -r group replicate fastq1 fastq2 control; do
    # Generate SampleID, Factor, and Condition from directory and group
    SampleID="${dir}_${group}_${replicate}"
    Factor="$group"
    Condition="${dir}"

    # Define Treatment, Tissue, and PeakCaller
    Treatment="none"
    Tissue="MSH2"
    PeakCaller="narrow"

    # Locate BAM file and Peak file paths
    bamReads="$dir/02_alignment/bowtie2/target/markdup/${group}_R${replicate}.target.markdup.sorted.bam"
    Peaks="$dir/03_peak_calling/04_called_peaks/macs2/${group}_R${replicate}.macs2_peaks.narrowPeak"

    # Check if files exist
    if [[ ! -f $bamReads ]]; then
      echo "Warning: BAM file not found: $bamReads" >&2
    fi
    if [[ ! -f $Peaks ]]; then
      echo "Warning: Peak file not found: $Peaks" >&2
    fi

    # Append row to output CSV
    echo "$SampleID,$Tissue,$Factor,$Condition,$Treatment,$replicate,$Tissue,$bamReads,$Peaks,$PeakCaller" >> $output_csv
  done
done

echo "Samplesheet generated: $output_csv"

