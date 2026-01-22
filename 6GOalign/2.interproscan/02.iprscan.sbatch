#!/usr/bin/bash
#SBATCH -p blade,gpu,himem,hugemem
#SBATCH -c 48

module load python3
module load java

/public4/software/interproscan/interproscan-5.71-102.0/interproscan.sh -cpu $SLURM_CPUS_PER_TASK  --goterms  -i Rform.prot_no_stop.fasta -f tsv -dp
