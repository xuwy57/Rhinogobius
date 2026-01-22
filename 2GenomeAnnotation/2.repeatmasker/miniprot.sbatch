#!/usr/bin/bash
#SBATCH -p long,gpu,himem,hugemem,blade
#SBATCH -c 24
#SBATCH --mem=100G

hostname

CPU=$SLURM_CPUS_PER_TASK

module load miniprot

miniprot -I -t $CPU --gff ref.mpi ./all_prot.fa > out.gff

