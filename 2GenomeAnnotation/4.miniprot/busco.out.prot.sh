#!/usr/bin/bash
#SBATCH -p blade,long,himem,hugemem,gpu
#SBATCH -c 16
source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate busco5.3.2

busco -c $SLURM_CPUS_PER_TASK -i out.prot.fa -m prot -o outpro -l actinopterygii_odb10 \
--offline --download_path /public2/shareddatabase/busco/download/ > busco.out.prot.log 2>&1
