#!/bin/bash
#SBATCH --job-name=fastp.%j
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=8gb
#SBATCH --qos=timgarrett
#SBATCH --account=timgarrett
#SBATCH -t 01:00:00
#SBATCH --output=/blue/zhangw/hkates/Tanzia_RNAseq/logs/fastp.%j.out
#SBATCH --error=/blue/zhangw/hkates/Tanzia_RNAseq/logs/fastp.%j.err
#SBATCH --array=12

echo This is task $SLURM_ARRAY_TASK_ID

echo -e "\nInfo: Starting a job on $(date) on $(hostname) in $(pwd).\n"

###################################################################
# This file is to perform quality filtering and trimming with fastp
###################################################################

module load fastp/0.23.4

# Define the base directory
BASE_DIR="/orange/zhangw/GE7607/GE-7607-WZhang-10B-22MNM5LT3-Lane2"

# Create an array of subdirectories dynamically
DIRS=$(ls -d ${BASE_DIR}/*/)

# Get the directory for the current task
DIR=$(echo $DIRS | tr ' ' '\n' | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Define the R1 and R2 files
Sample1=$(ls ${DIR}/*_R1_001.fastq.gz)
Sample2=$(ls ${DIR}/*_R2_001.fastq.gz)

# Define the output directory
OUT_DIR="/blue/zhangw/hkates/Tanzia_RNAseq/results/fastp"

# Extract sample name from the input file path
SAMPLE_NAME=$(basename ${Sample1} | sed 's/_R1_001.fastq.gz//')

# Run fastp with best practices
fastp -i ${Sample1} -I ${Sample2} \
      -o ${OUT_DIR}/${SAMPLE_NAME}_R1_trimmed.fastq.gz \
      -O ${OUT_DIR}/${SAMPLE_NAME}_R2_trimmed.fastq.gz \
      -h ${OUT_DIR}/${SAMPLE_NAME}_fastp_report.html \
      -j ${OUT_DIR}/${SAMPLE_NAME}_fastp_report.json \
      -q 20 -u 30 -n 5 \
      --cut_front --cut_tail --cut_right -w 3

echo -e "\nInfo: fastp job completed on $(date) on $(hostname).\n"

