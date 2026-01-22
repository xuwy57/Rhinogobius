#!/usr/bin/bash
#SBATCH -c 64
#SBATCH -J cafe
#SBATCH --mem=100G

CAFE=/public3/group_crf/software/CAFE5/bin/cafe5
CPU=64

for i in `seq 1 20`; do
	mkdir -p runresults_k4/$i
$CAFE -c $CPU -t tree.tre -i orthogroups.genecount.txt -p -eresults/Base_error_model.txt -k 4 -o  runresults_k4/$i >  runresults_k4/$i/est.log 2>&1 
done
