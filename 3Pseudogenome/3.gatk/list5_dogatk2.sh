#!/bin/bash
#do gatk2.sh for each 7 samples----sample.bam file
for sSample in /data/projects/zma/Rhinogobius_spp/DNA/bam/S1_R-szechuanensis_Sichuan.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/Z10_R-sp_Zhejianglongquanshizhulongzhen.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/Z11_R-sp_Zhejiangchangshanxian.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/Z12_R-davidi_Zhejiangchangshanxian.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/Z13_R-sp_Zhejiangkaihuaxiansuzhuangzhen.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/Z14_R-lentiginis_Zhejiangyuyaoliangnongzhen.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/Z15_R-leavelli_Zhejiangchangshanxian.bam\
	/data/projects/zma/Rhinogobius_spp/DNA/bam/Z16_R-cf-lentiginis_Zhejiangchangshanxian.bam
do
#	echo $sSample
        sBase=`basename $sSample`
        sName=${sBase/.bam/}
#	echo $sName
	cd /data/projects/zma/Rhinogobius_spp/DNA/bam/gatk/$sName
        source /data/projects/zma/Rhinogobius_spp/DNA/bam/gatk/03_gatk2.sh $sSample &

done;
wait
