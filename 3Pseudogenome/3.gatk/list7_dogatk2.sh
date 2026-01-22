#!/bin/bash
#do gatk2.sh for each 9 samples----sample.bam file
for sSample in /data/projects/zma/Rhinogobius_spp/DNA/bam/Z3_R-wuyiensis_Zhejianglongquanshijinxizhen.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/Z4_R-cf-changtingensis_Zhejianglongquanshizhulongz.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/Z5_R-niger_Zhejianglongquanshizhulongzhen.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/Z6_R-giurinus_Zhejianglongquanshizhulongzhen.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/Z7_R-davidi_Zhejianganjixiantianhuangpingzhen.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/Z8_R-immaculatus_Zhejiangtongluxianlucixiang.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/Z9_R-sp_Zhejiangkaihuaxiansuzhuangzhen.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/H1_R-sp_Hunanlianyuanshimaotangzhen.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/GX2_R-sp_Guangxifangchenggangmaluzhen.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/HN2_R-wanchuangensis_Hainanqionghaishihuishanzhen.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/Z15_R-leavelli_Zhejiangchangshanxian.bam
do
#	echo $sSample
        sBase=`basename $sSample`
        sName=${sBase/.bam/}
#	echo $sName
	cd /data/projects/zma/Rhinogobius_spp/DNA/bam/gatk/$sName
        source /data/projects/zma/Rhinogobius_spp/DNA/bam/gatk/03_gatk2.sh $sSample &

done;
wait
