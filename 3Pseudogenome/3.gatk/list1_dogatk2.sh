#!/bin/bash
#do gatk2.sh for each 7 samples----sample.bam file
for sSample in /data/projects/zma/Rhinogobius_spp/DNA/bam/A1_R-sp_Anhuidongzhixiannixizhen.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/A3_R-sp_Anhuiqimendatanxiang.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/A6_R-lentiginis_Anhuishitai.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/A7_R-sp_Anhuishitai.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/A8_R-cf-wuyiensis_Anhuiqimenlikouzhen.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/F1_R-cf-rubromaculatus_Fujianputianshi.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/F2_R-xianshuiensis_Fujianputianshi.bam
do
#	echo $sSample
        sBase=`basename $sSample`
        sName=${sBase/.bam/}
#	echo $sName
	cd /data/projects/zma/Rhinogobius_spp/DNA/bam/gatk/$sName
        source /data/projects/zma/Rhinogobius_spp/DNA/bam/gatk/03_gatk2.sh $sSample &

done;
wait
