#!/usr/bin/bash
#SBATCH -p gpu
#SBATCH -c 64
#SBATCH --mem=150G
source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate funannotate1.8.15

source genomedef.sh

CPU=$SLURM_CPUS_PER_TASK

funannotate train -i $GENOME -o trained  \
--left /public3/group_crf/home/g21mazhy26/R.formosanus_data/02_clean.data/trimmomatic/RNA-2021-11-16/R2120716_paired_raw_1.fq.gz /public3/group_crf/home/g21mazhy26/R.formosanus_data/02_clean.data/trimmomatic/RNA-2022-01-11/R2124108_paired_raw_1.fq.gz \
--right /public3/group_crf/home/g21mazhy26/R.formosanus_data/02_clean.data/trimmomatic/RNA-2021-11-16/R2120716_paired_raw_2.fq.gz /public3/group_crf/home/g21mazhy26/R.formosanus_data/02_clean.data/trimmomatic/RNA-2022-01-11/R2124108_paired_raw_2.fq.gz  \
--no_trimmomatic \
--no_normalize_reads \
--pasa_min_avg_per_id 85 \
--jaccard_clip \
--memory 100G \
--trinity "$TRINITY_TRANSCRIPT" \
--species "$SPECIES" \
--pasa_db mysql \
--max_intronlen 200000 \
--aligners minimap2 minimap2 \
--strain "$STRAIN" --cpus $CPU 
