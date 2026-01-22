#!/usr/bin/bash
#SBATCH -p blade
#SBATCH -c 20
#SBATCH --mem=40G

module load samtools

for i in ./*_pseudogenome2.fasta
do
	samtools faidx $i 
done
wait
