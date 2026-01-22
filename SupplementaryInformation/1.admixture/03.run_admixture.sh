Name=Pomcan
admixture=/data/software/admixture_linux-1.3.0/admixture
for K in 2 3 4 5 6 7 8 9 10; do
	$admixture --cv $Name.bed $K -j12 | tee log${K}.out 
done
