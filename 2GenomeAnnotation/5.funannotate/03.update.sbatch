#!/usr/bin/bash
#SBATCH -p long,blade,hugemem,gpu
#SBATCH -c 32
#SBATCH --mem=50G

source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate funannotate1.8.15

funannotate update -i traineid_miniprot --cpus 32 --pasa_db mysql --alt_transcripts 0.3 \
--left /public3/group_crf/home/g21mazhy26/R.formosanus_data/02_clean.data/trimmomatic/RNA-2021-11-16/R2120716_paired_raw_1.fq.gz /public3/group_crf/home/g21mazhy26/R.formosanus_data/02_clean.data/trimmomatic/RNA-2022-01-11/R2124108_paired_raw_1.fq.gz \
--right /public3/group_crf/home/g21mazhy26/R.formosanus_data/02_clean.data/trimmomatic/RNA-2021-11-16/R2120716_paired_raw_2.fq.gz /public3/group_crf/home/g21mazhy26/R.formosanus_data/02_clean.data/trimmomatic/RNA-2022-01-11/R2124108_paired_raw_2.fq.gz  \
--no_trimmomatic \
--no_normalize_reads \
