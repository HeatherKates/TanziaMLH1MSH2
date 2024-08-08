#!/bin/bash
#SBATCH --job-name=rserver
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --time=10:00:00
#SBATCH --output=rserver_%j.log
#SBATCH --error=rserver_%j.error
#SBATCH --account=timgarrett
#SBATCH --qos=timgarrett

module purge; module load R/4.2
rserver
