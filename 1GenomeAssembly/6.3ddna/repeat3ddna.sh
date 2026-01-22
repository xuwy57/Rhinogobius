#!/bin/bash
#SBATCH -p himem
#SBATCH -c 1

Contigs=/public3/group_crf/home/g21mazhy26/R.formosanus_data/11_3ddna+juicebox/references/ref.fa
MND=/public3/group_crf/home/g21mazhy26/R.formosanus_data/11_3ddna+juicebox/aligned/merged_nodups.txt
MAPQ=30
GAPSIZE=1000

module load java
/public3/group_crf/software/3d-dna-201008/run-asm-pipeline-post-review.sh -r genome.review.assembly --sort-output -q $MAPQ $Contigs $MND > genome_3ddna.log 2>&1
