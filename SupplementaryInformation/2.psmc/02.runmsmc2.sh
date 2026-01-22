MSMC2=/data/software/msmc2-2.1.3/build/release/msmc2

sSampleList="samples.txt"
arrPops=( `cut -d" " -f1 $sSampleList` )
sIn=formsmc2_in
sOut=msmc2ret
recOverMu=1

mkdir -p $sOut


for sPop in "${arrPops[@]}"; do
	#sFiles="${sIn}/${sPop}*/*/formsmc2.multihetsep.txt"
	sFiles=`realpath ${sIn}/${sPop}*/*/formsmc2.multihetsep.txt` #do not include sex chromosome
	sCmd="$MSMC2 -o $sOut/$sPop -r $recOverMu -i 50 -t 24 $sFiles > /dev/null "
	echo $sCmd
	eval $sCmd > /dev/null 2>&1 & 
done

wait

