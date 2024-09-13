#!/bin/bash
#SBATCH --job-name=fastqc.%j
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --qos=timgarrett-b
#SBATCH --account=timgarrett
#SBATCH -t 02:00:00
#SBATCH --output=../logs/fastqc.%j.out
#SBATCH --error=../logs/fastqc.%j.err
#SBATCH --array=1-12

echo "This is task $SLURM_ARRAY_TASK_ID"

echo -e "\nInfo: Starting a job on $(date) on $(hostname) in $(pwd).\n"

###################################################################
# This file is to perform FastQC quality control per lane
###################################################################

#!/bin/bash

echo "This is task $SLURM_ARRAY_TASK_ID"

echo -e "\nInfo: Starting a job on $(date) on $(hostname) in $(pwd).\n"

###################################################################
# This file is to perform FastQC quality control on per-lane samples
###################################################################

module load fastqc/0.12.1

# Define the base directory
BASE_DIR="/orange/zhangw/NS3842/NS-3842-WZhang-15B-223353LT1-Lane1-2"

# Get the list of all directories (with lane info)
SAMPLES=($(ls -d ${BASE_DIR}/*))

# Convert the list of directories into an array
SAMPLES_ARRAY=("${SAMPLES[@]}")

# Get the directory for the current task
SAMPLE_DIR=${SAMPLES_ARRAY[$SLURM_ARRAY_TASK_ID]}

# Define the output directory
OUT_DIR="../results/fastqc/$(basename ${SAMPLE_DIR})"
mkdir -p ${OUT_DIR}

# Define the R1 and R2 files
Sample1_R1=$(ls ${SAMPLE_DIR}/*_R1_001.fastq.gz)
Sample2_R2=$(ls ${SAMPLE_DIR}/*_R2_001.fastq.gz)

# Run FastQC for the current sample
fastqc ${Sample1_R1} ${Sample2_R2} -o ${OUT_DIR}

echo -e "\nInfo: FastQC job completed for ${SAMPLE_DIR} on $(date) on $(hostname).\n"
