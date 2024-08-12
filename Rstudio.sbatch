#!/bin/bash
#SBATCH --job-name=rserver
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=96:00:00
#SBATCH --output=rserver_%j.log
#SBATCH --error=rserver_%j.error
#SBATCH --account=timgarrett
#SBATCH --qos=timgarrett-b

module purge; module load R/4.2
rserver
