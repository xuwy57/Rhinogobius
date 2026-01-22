#!/bin/bash
#do gatk2.sh for each 8 samples----sample.bam file

for sSample in /data/projects/zma/Rhinogobius_spp/DNA/bam/Z13_R-sp_Zhejiangkaihuaxiansuzhuangzhen.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/Z16_R-cf-lentiginis_Zhejiangchangshanxian.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/Z17_R-cliffordpopei_Zhejiangchangshanxian.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/Z18_R-cliffordpopeiR-leavelli_Zhejiangchangshanxia.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/Z20_R-aporus_Zhejiangwuyixianxilianxiang.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/Z1_R-sp_Zhejianglongquanshizhulongzhen.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/Z4_R-cf-changtingensis_Zhejianglongquanshizhulongz.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/Z6_R-giurinus_Zhejianglongquanshizhulongzhen.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/Z7_R-davidi_Zhejianganjixiantianhuangpingzhen.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/Z8_R-immaculatus_Zhejiangtongluxianlucixiang.bam
do
#	echo $sSample
        sBase=`basename $sSample`
        sName=${sBase/.bam/}
#	echo $sName
	cd /data/projects/zma/Rhinogobius_spp/DNA/bam/gatk/$sName
        source /data/projects/zma/Rhinogobius_spp/DNA/bam/gatk/03_gatk2.sh $sSample &
done;
wait
