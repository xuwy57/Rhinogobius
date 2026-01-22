
CAFE=/data/software/CAFE5/bin/cafe5
CPU=32

for i in `seq 18 20`; do
	mkdir -p runresults_k5/$i
$CAFE -c $CPU -t tree.tre -i orthogroups.genecount.txt -p -eresults/Base_error_model.txt -k 5 -o  runresults_k5/$i >  runresults_k5/$i/est.log 2>&1 
done
