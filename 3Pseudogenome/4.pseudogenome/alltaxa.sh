echo -n "" > alltaxa.txt
for i in /public3/group_crf/home/g21mazhy26/new_pseudogenome/gvcf/*.bam.genotyped.g.vcf.gz; do
	sDir=`dirname $i`
	sBase=`basename $i`
	sSample=${sBase/.bam.genotyped.g.vcf.gz/}
	echo -ne "$sSample\t$i\n" >> alltaxa.txt

#	echo $sSample
#	echo $sBase
done
