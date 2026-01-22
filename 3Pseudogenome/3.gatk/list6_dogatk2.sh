#!/bin/bash
#do gatk2.sh for each 7 samples----sample.bam file
for sSample in /data/projects/zma/Rhinogobius_spp/DNA/bam/Z17_R-cliffordpopei_Zhejiangchangshanxian.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/Z18_R-cliffordpopeiR-leavelli_Zhejiangchangshanxia.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/Z19_R-aporus_Zhejiangwuyixianchengqu.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/Z20_R-aporus_Zhejiangwuyixianxilianxiang.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/Z21_R-sp_Zhejiangchangshanxiansianzhen.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/Z1_R-sp_Zhejianglongquanshizhulongzhen.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/Z2_R-aporus_Zhejiangpananxianshuangxizhen.bam
do
#	echo $sSample
        sBase=`basename $sSample`
        sName=${sBase/.bam/}
#	echo $sName
	cd /data/projects/zma/Rhinogobius_spp/DNA/bam/gatk/$sName
        source /data/projects/zma/Rhinogobius_spp/DNA/bam/gatk/03_gatk2.sh $sSample &

done;
wait
