#!/bin/bash
#SBATCH -p hugemem
#SBATCH -o slurmlog/addslurm-%a.out
#SBATCH -e slurmlog/addslurm-%a.err
#SBATCH -a 0-50
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1


module load bcftools
CPU=${SLURM_CPUS_PER_TASK} #per task
nCount=0
nTotalJobs=${SLURM_ARRAY_TASK_COUNT}
nThisJob=${SLURM_ARRAY_TASK_ID}
nSampleLines=921612192 #total length of full chromosomes

for i in /public3/group_crf/home/g21mazhy26/new_pseudogenome/gvcf/*.bam.genotyped.g.vcf.gz; do
        nModulo=$(( $nCount % $nTotalJobs ))
        if (( $nModulo == $nThisJob )); then

        sID=`basename $i`
        sID=${sID/.bam.genotyped.g.vcf.gz/}
        echo $sID

        bash get_mean_dp.sh  $i $sID $nSampleLines
        fi
        nCount=$(( $nCount + 1 ))

done



