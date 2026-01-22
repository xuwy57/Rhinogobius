#!/usr/bin/bash
#SBATCH -p blade,gpu,himem,hugemem
#SBATCH -c 64
#SBATCH -J cafe
#SBATCH --mem=100G

CAFE=/public3/group_crf/software/CAFE5/bin/cafe5

CPU=64

for i in `seq 1 20`; do
	mkdir -p runresults_k6/$i
$CAFE -c $CPU -t tree.tre -i orthogroups.genecount.txt -p -eresults/Base_error_model.txt -k 6 -o  runresults_k6/$i >  runresults_k6/$i/est.log 2>&1 
done
