#!/usr/bin/bash
#SBATCH -p blade,gpu,himem,hugemem
#SBATCH -c 64
#SBATCH -J cafe
#SBATCH --mem=100G

CAFE=/public3/group_crf/software/CAFE5/bin/cafe5
CPU=64

$CAFE -c $CPU -t tree.tre -i orthogroups.genecount.txt -p -e > err.est.log 2>&1 
