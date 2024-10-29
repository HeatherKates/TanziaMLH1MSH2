#!/bin/bash
#SBATCH --job-name=fastqc.%j
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --qos=REPLACE_WITH_YOUR_QOS #Replace with your hipergator qos
#SBATCH --account=REPLACE_WITH_YOUR_ACCOUNT #Replace with your hipergator account
#SBATCH -t 02:00:00
#SBATCH --output=logs/fastqc.%j.out
#SBATCH --error=logs/fastqc.%j.err
#SBATCH --array=1-12

echo "This is task $SLURM_ARRAY_TASK_ID"

echo -e "\nInfo: Starting a job on $(date) on $(hostname) in $(pwd).\n"

###################################################################
# This file is to perform FastQC quality control for two lanes per sample
###################################################################

module load fastqc/0.12.1

# Define the base directory
BASE_DIR="/orange/zhangw/NS3842/NS-3842-WZhang-15B-223353LT1-Lane1-2"

# Create an array of subdirectories dynamically (one per sample)
DIRS=$(ls -d ${BASE_DIR}/*_L001/)

# Get the directory for the current task
DIR=$(echo $DIRS | tr ' ' '\n' | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Extract the base sample name (without the lane information)
SAMPLE_BASE=$(basename $DIR | sed 's/_L001//')

# Define the R1 and R2 files for both lanes
Sample1_L001=$(ls ${BASE_DIR}/${SAMPLE_BASE}_L001/*_R1_001.fastq.gz)
Sample2_L001=$(ls ${BASE_DIR}/${SAMPLE_BASE}_L001/*_R2_001.fastq.gz)

Sample1_L002=$(ls ${BASE_DIR}/${SAMPLE_BASE}_L002/*_R1_001.fastq.gz)
Sample2_L002=$(ls ${BASE_DIR}/${SAMPLE_BASE}_L002/*_R2_001.fastq.gz)

# Define the output directory
OUT_DIR="../results/1_fastqc/${SAMPLE_BASE}"
mkdir -p ${OUT_DIR}

# Run FastQC on both lanes
fastqc ${Sample1_L001} ${Sample2_L001} -o ${OUT_DIR} -t 2 
fastqc ${Sample1_L002} ${Sample2_L002} -o ${OUT_DIR} -t 2

echo -e "\nInfo: FastQC job completed for ${SAMPLE_BASE} on $(date) on $(hostname).\n"
