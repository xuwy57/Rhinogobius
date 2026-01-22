MSMC2=/data/software/msmc2-2.1.3/build/release/msmc2

nReps=30
nChunkSize=5000000


arrPops=( `ls formsmc2_in` )

sOutDir=bootstrapped
nCPUs=128
nCPUPerProcess=2
nUsedCPU=0


for sPop in "${arrPops[@]}"; do
	for i in $(seq 1 $nReps); do 

		sBSFolder="$sOutDir/$sPop/_${i}";
		[ ! -s $sBSFolder/out.final.txt ] && $MSMC2 -o $sBSFolder/out -i 50 -t $nCPUPerProcess -r 1 $sBSFolder/bootstrap_multihetsep.*.txt 2> /dev/null &
		nUsedCPU=$(( nUsedCPU + nCPUPerProcess ));
                #echo $nUsedCPU
                if [ $nUsedCPU -ge $nCPUs ]; then
                        nUsedCPU=0
                        echo Wait...
                        wait;
                fi


	done

done

wait




