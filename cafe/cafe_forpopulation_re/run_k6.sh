
CAFE=/data/software/CAFE5/bin/cafe5
CPU=84

for i in `seq 12 20`; do
	mkdir -p runresults_k6/$i
$CAFE -c $CPU -t tree.tre -i orthogroups.genecount.txt -p -eresults/Base_error_model.txt -k 6 -o  runresults_k6/$i >  runresults_k6/$i/est.log 2>&1 
done
