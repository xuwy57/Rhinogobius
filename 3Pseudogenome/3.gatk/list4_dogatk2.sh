#!/bin/bash
#do gatk2.sh for each 7 samples----sample.bam file
for sSample in /data/projects/zma/Rhinogobius_spp/DNA/bam/HN3_R-leavelli_Hainanzhanzhoudachengzhen.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/HN4_R-leavelli_Hainanbaotingshilingzhen.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/HN5_R-linshuiensis_Hainanbaotingshilingzhen.bam\
	/data/projects/zma/Rhinogobius_spp/DNA/bam/HN6_R-sp_Hainanhaikou.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/J1_R-sp_Jiangxipingxiangshianyuanqu.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/J2_R-sp_Jiangxiwuyuan.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/L1_R-cf-brunneus_Liaoningkuandian.bam
do
#	echo $sSample
        sBase=`basename $sSample`
        sName=${sBase/.bam/}
#	echo $sName
	cd /data/projects/zma/Rhinogobius_spp/DNA/bam/gatk/$sName
        source /data/projects/zma/Rhinogobius_spp/DNA/bam/gatk/03_gatk2.sh $sSample &

done;
wait
