vcfcaller=/data/software/msmc-tools/vcfAllSiteParser.py

CONVERT=/data/software/msmc-tools/generate_multihetsep.py
BCFTOOLS="bcftools"
sVCFFolder=./all_gvcf/
sVCFSuffix=.genotyped.g.vcf.gz
sChrSuffix="Chromosome"
sVCFFolder=`realpath $sVCFFolder`

#sSampleList="samples.txt"
sSampleList="samples.txt"
arrSamples=( `cut -d" " -f1 $sSampleList` )
arrCovUpper=(  `cut -d" " -f4 $sSampleList` )
arrCovLower=(  `cut -d" " -f3 $sSampleList` )

arrChr=( `seq 1 22` )

sOutDir=formsmc2_in/


for nSample in "${!arrSamples[@]}"; do 
	sSample=${arrSamples[$nSample]};
	nCovUpper=${arrCovUpper[$nSample]};
	nCovLower=${arrCovLower[$nSample]};
	sFilterCmd=${arrFilterCmd[$nSample]};

	for nChr in "${arrChr[@]}"; do 
		sChr=$sChrSuffix$nChr
		sChrDir=$sOutDir/$sSample/$sChr
		mkdir -p $sChrDir
		( if [ ! -s $sChrDir/out_mask.bed.gz ];then $BCFTOOLS view -r "$sChr"  -e ' (GT=="RR" && RGQ<30) || GT=="mis" || ( GT!="RR" && (QUAL<30 || GQ<30 || INFO/DP>'$nCovUpper' || INFO/DP<'$nCovLower' || TYPE!="snp" )) ' $sVCFFolder/$sSample$sVCFSuffix | $vcfcaller $sChr $sChrDir/out_mask.bed.gz | gzip -c > $sChrDir/out.vcf.gz; fi; \
		if [ ! -e $sChrDir/done.txt ]; then $CONVERT --mask $sChrDir/out_mask.bed.gz $sChrDir/out.vcf.gz > $sChrDir/formsmc2.multihetsep.txt 2> $sChrDir/formsmc2.multihetsep.log && touch $sChrDir/done.txt; fi; ) 2>$sChrDir/log.txt &
	done

done

wait

