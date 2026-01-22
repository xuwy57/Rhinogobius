#!/bin/bash
#SBATCH -p hugemem,blade,gpu,himem
#SBATCH -o slurmlog/norm-%a.out
#SBATCH -e slurmlog/norm-%a.err
#SBATCH -a 0-50
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1


module load bcftools
CPU=${SLURM_CPUS_PER_TASK} #per task
nCount=0
nTotalJobs=${SLURM_ARRAY_TASK_COUNT}
nThisJob=${SLURM_ARRAY_TASK_ID}

for i in /public4/group_crf/shareddata/Rformosanus_gatk/Rformosanus.genotyped.g.vcf.gz; do
        nModulo=$(( $nCount % $nTotalJobs ))
        if (( $nModulo == $nThisJob )); then

        sID=`basename $i`
        sID=${sID/.genotyped.g.vcf.gz/}
        echo $sID

        bash normalize_site_dp.sh $i $sID
        fi
        nCount=$(( $nCount + 1 ))

done



