trimjar=/data/software/Trimmomatic-0.39/trimmomatic-0.39.jar
#adapterfa=/data/software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa
adapterfa=/data/software/Trimmomatic-0.39/adapters/NexteraPE-PE.fa

for sR1 in /data/projects/shareddata/rhinogobius_spp/raw_data_2022_11_06/Rawdata/*/*_R1.fq.gz
do
	sStem=`basename $sR1`;
	sStem=${sStem/%_R1.fq.gz}
	sPath=`dirname $sR1`
	sR2=$sPath/${sStem}_R2.fq.gz
	java -jar $trimjar PE -threads 24 -phred33 $sR1 $sR2 \
		$sStem.paired_1.fq.gz $sStem.unpaired_1.fq.gz \
	       	$sStem.paired_2.fq.gz $sStem.unpaired_2.fq.gz \
		ILLUMINACLIP:$adapterfa:2:30:10 \
		LEADING:3 TRAILING:3 SLIDINGWINDOW:5:20 MINLEN:50 > trim.$sStem.log 2>&1 &
done
wait;


