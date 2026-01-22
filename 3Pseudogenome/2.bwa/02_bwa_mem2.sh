CPU=4
REF=./ref_Rform.fasta
for sR1 in /public3/group_crf/home/g21mazhy26/Rhinogobius_spp/DNA/GATK/trim_out2/*_1.fq.gz; do
        sDir=`dirname $sR1`
	sBase=`basename $sR1`
	sSample=${sBase/.paired_1.fq.gz/}
	sR2=${sDir}/${sSample}.paired_2.fq.gz
	( bwa mem -t $CPU $REF $sR1 $sR2 \
        | samtools view  -u - | samtools sort - -m 20g -o $sSample.bam ) > $sSample.log 2>&1 &

#	echo $sSample
#	echo $sR1
#	echo $sR2

done
wait


#bwa比对 生成bam文件，并排序

