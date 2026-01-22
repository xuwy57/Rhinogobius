for sR1 in /data2/projects/zma/Rhinogobius_spp/RNA/02_cleandata/trim_out/*.paired_2.fq.gz; do
	sDir=`dirname $sR1`
	sBase=`basename $sR1`
	sSample=${sBase/.paired_2.fq.gz/}
#	echo $sSample
#	echo $sR1
        mkdir -p $sSample
done
