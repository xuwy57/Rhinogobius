#!/usr/bin/bash
#SBATCH -p long
#SBATCH -c 12
#SBATCH --mem=50G

ASSEMBLY=/data/projects/zma/R.formosanus_data/11_3ddna+juicebox/3ddna/Rform_sorted.fa # the assembly to look at
SP=Rform_hic
REPMODEL="singularity run docker://dfam/tetools:latest " #command of repeatmodeler
THREADS=48

ln -sf $ASSEMBLY scf_sorted.fa

$REPMODEL BuildDatabase -name $SP -engine ncbi scf_sorted.fa

$REPMODEL RepeatModeler -LTRStruct -engine ncbi -threads $THREADS -database $SP
