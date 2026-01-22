#!/usr/bin/bash
#SBATCH -p long
#SBATCH -c 12
#SBATCH --mem=50G

ASSEMBLY=/data/projects/zma/R.formosanus_data/11_3ddna+juicebox/3ddna/Rform.fasta # the assembly to look at
SP=Rform
REPMODEL="singularity run docker://dfam/tetools:latest " #command of repeatmodeler
THREADS=48

ln -sf $ASSEMBLY scf.fa

$REPMODEL BuildDatabase -name $SP -engine ncbi scf.fa

$REPMODEL RepeatModeler -LTRStruct -engine ncbi -threads $THREADS -database $SP
