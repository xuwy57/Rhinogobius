#!/usr/bin/bash
#SBATCH -c 40
#SBATCH --mem=64G

for sR1 in /data2/projects/zma/Rhinogobius_spp/RNA/02_cleandata/trim_all/*_all.paired_1.fq.gz
do
	sStem=`basename $sR1`;
	sStem=${sStem/_all.paired_1.fq.gz}
	sPath=`dirname $sR1`
	sR2=$sPath/${sStem}_all.paired_2.fq.gz

	cd /data2/projects/zma/Rhinogobius_spp/RNA/04_STAR/${sStem}

#	STAR --runThreadN 40 --runMode alignReads --readFilesCommand zcat --twopassMode Basic --outSAMtype BAM Unsorted --outSAMunmapped None --genomeDir /data2/projects/zma/Rhinogobius_spp/RNA/04_STAR/ref --readFilesIn $sR1 --outFileNamePrefix DiffExpress

	STAR --runThreadN 40 --runMode alignReads --readFilesCommand zcat --twopassMode Basic --outSAMtype BAM Unsorted --outSAMunmapped None --genomeDir /data2/projects/zma/Rhinogobius_spp/RNA/04_STAR/ref --readFilesIn $sR1 $sR2 --outFileNamePrefix $sStem

done
wait;
