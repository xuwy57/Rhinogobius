#!/usr/bin/bash
#SBATCH -p blade,gpu
#SBATCH -c 64
#SBATCH -J genespace
#SBATCH --mem=200G

source /public/apps/miniconda3/etc/profile.d/conda.sh

conda activate orthofinder
module load mcscanx
module load R/4.2.0

Rscript run.R
