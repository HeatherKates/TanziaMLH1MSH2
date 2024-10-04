#!/bin/bash
module load multiqc
cd ../
# Run multiqc to summarize the fastqc results
multiqc 2_fastQC/ -o 2_fastQC/
