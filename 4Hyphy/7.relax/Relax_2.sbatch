#!/bin/bash
#SBATCH -p blade,himem,gpu
#SBATCH -o slurmlog/RelaxNothos_slurm-%a.out
#SBATCH -e slurmlog/RelaxNothos_slurm-%a.err
#SBATCH -a 0-63
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

module load hyphy/2.5.32
module load clusterbasics

#SET THIS TO THE CORRECT NUMBER
TreeFolder="/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/hyphy_1016/exportAln_genespace/Relax_Clade2"
OutDIR=`basename $TreeFolder`;


mkdir -p $OutDIR

php hyphyrelax.php -N ${SLURM_ARRAY_TASK_COUNT} -f ${SLURM_ARRAY_TASK_ID} -T $TreeFolder -o $OutDIR 

exit

