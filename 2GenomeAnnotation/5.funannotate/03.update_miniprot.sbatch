#!/usr/bin/bash
#SBATCH -p long
#SBATCH -c 32
#SBATCH --mem=80G

source genome_miniprot.sh

source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate funannotate1.8.15

CPU=$SLURM_CPUS_PER_TASK

funannotate update -i traineid_miniprot --cpus $CPU --pasa_db mysql --alt_transcripts 0.3
