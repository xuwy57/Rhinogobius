#!/bin/bash
#do gatk2.sh for each 8 samples----sample.bam file
for sSample in /data/projects/zma/Rhinogobius_spp/DNA/bam/F3_R-reticulatus_Fujianpuchengshifulingzhen.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/F4_R-reticulatus_Fujianfuzhou.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/F7b_R-sp_Fujianyunxiao.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/GD1_R-duospilus_Guangdonglongmen.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/GX2_R-sp_Guangxifangchenggangmaluzhen.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/HN2_R-wanchuangensis_Hainanqionghaishihuishanzhen.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/S1_R-szechuanensis_Sichuan.bam \
	/data/projects/zma/Rhinogobius_spp/DNA/bam/Z10_R-sp_Zhejianglongquanshizhulongzhen.bam 
do
#	echo $sSample
        sBase=`basename $sSample`
        sName=${sBase/.bam/}
#	echo $sName
	cd /data/projects/zma/Rhinogobius_spp/DNA/bam/gatk/$sName
        source /data/projects/zma/Rhinogobius_spp/DNA/bam/gatk/03_gatk2.sh $sSample &

done;
wait
