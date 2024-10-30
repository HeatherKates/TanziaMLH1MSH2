#!/bin/bash
#SBATCH --job-name=fastqc.%j
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --qos=timgarrett-b
#SBATCH --account=timgarrett
#SBATCH -t 02:00:00
#SBATCH --output=logs/fastqc.%j.out
#SBATCH --error=logs/fastqc.%j.err
#SBATCH --array=1-1

echo "This is task $SLURM_ARRAY_TASK_ID"

echo -e "\nInfo: Starting a job on $(date) on $(hostname) in $(pwd).\n"

###################################################################
# This file is to perform FastQC quality control for two lanes per sample
###################################################################

module load fastqc/0.12.1

# Define the base directory
BASE_DIR="/orange/zhangw/GE7607/GE-7607-WZhang-10B-22MNM5LT3-Lane2"

# Create an array of subdirectories dynamically (one per sample)
DIRS=$(ls -d ${BASE_DIR}/*/)

# Get the directory for the current task
DIR=$(echo $DIRS | tr ' ' '\n' | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Define the R1 and R2 files
Sample1=$(ls ${DIR}/*_R1_001.fastq.gz)
Sample2=$(ls ${DIR}/*_R2_001.fastq.gz)

# Define the output directory
OUT_DIR="../results/1_fastqc"
mkdir -p ${OUT_DIR}

# Extract sample name from the input file path
SAMPLE_NAME=$(basename ${Sample1} | sed 's/_R1_001.fastq.gz//')

# Run FastQC on R1 and R2 files
fastqc ${Sample1} ${Sample2} -o ${OUT_DIR} -t 2

echo -e "\nInfo: FastQC job completed for ${SAMPLE_NAME} on $(date) on $(hostname).\n"
