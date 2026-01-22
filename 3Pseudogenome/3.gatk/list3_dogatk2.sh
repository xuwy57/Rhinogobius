#!/bin/bash
#do gatk2.sh for each 7 samples----sample.bam file
for sSample in /data/projects/zma/Rhinogobius_spp/DNA/bam/GX1_R-sp_Guangxinanning.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/GX2_R-sp_Guangxifangchenggangmaluzhen.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/GX3_R-sp_Guangxifangchenggangmaluzhen.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/H2_R-sp_Hunanloudilouxingqu.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/HB1_R-houheensis_Hubeixuandongxian.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/HN1_R-nandujiangensis_Hainanhaikoushijiuzhouzhen.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/HN2_R-wanchuangensis_Hainanqionghaishihuishanzhen.bam
do
#	echo $sSample
        sBase=`basename $sSample`
        sName=${sBase/.bam/}
#	echo $sName
	cd /data/projects/zma/Rhinogobius_spp/DNA/bam/gatk/$sName
        source /data/projects/zma/Rhinogobius_spp/DNA/bam/gatk/03_gatk2.sh $sSample &

done;
wait
