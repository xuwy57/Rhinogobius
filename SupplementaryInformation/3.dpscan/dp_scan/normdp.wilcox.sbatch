#!/bin/bash
#SBATCH -p blade,gpu,himem,hugemem
#SBATCH -o slurmlog/wk-%a.out
#SBATCH -e slurmlog/wk-%a.err
#SBATCH -a 0-99
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16


module load clusterbasics
module load R
nTotalJobs=${SLURM_ARRAY_TASK_COUNT}
nThisJob=${SLURM_ARRAY_TASK_ID}

php compute_depth_wilcoxp.php $nTotalJobs $nThisJob
