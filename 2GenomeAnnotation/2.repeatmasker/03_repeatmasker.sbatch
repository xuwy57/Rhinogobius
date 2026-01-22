#!/usr/bin/bash
#SBATCH -c 40
#SBATCH -p long
#SBATCH --mem=100G

ASSEMBLY=/data/projects/zma/R.formosanus_data/11_3ddna+juicebox/3ddna/Rform_sorted.fa #hic assembly
REPBASELIB=`realpath ./all_hic_te.fa` #repeat library, with killifish repeats
THREADS=40

ln -sf $ASSEMBLY scf.fa


singularity run docker://dfam/tetools:latest RepeatMasker -gff -html -pa $THREADS -lib $REPBASELIB  scf.fa

