#!/bin/bash
#SBATCH --job-name=bcf_proc
#SBATCH --cpus-per-task=2        # 每个任务的线程数
#SBATCH --mem=8G                 # 每个任务内存
#SBATCH --array=1-51             # 样本数量
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err

module load bcftools
genome=/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/ref/ref/Rform.fasta

# 样本列表（对应 array index）
sample_list=(./gvcf/*.g.vcf.gz)
vcf_file=${sample_list[$SLURM_ARRAY_TASK_ID-1]}
sample_name=$(basename "$vcf_file" .g.vcf.gz)

# 创建输出目录
mkdir -p tmp_vcf norm_vcf2

# Step 1: 过滤 SNP + indel，并限制 ploidy
#bcftools view \
#    -O z \
#    -M 2 -m 2 \
#    -i 'TYPE=="snp" || TYPE=="indel"' \
#    --threads $SLURM_CPUS_PER_TASK \
#    "$vcf_file" > tmp_vcf/${sample_name}.filtered.vcf.gz

# 建立索引
#bcftools index tmp_vcf/${sample_name}.filtered.vcf.gz

# Step 2: 归一化
bcftools norm \
    -O z \
    -f "$genome" \
    -w 99999999 \                   # 可以调小，减少耗时
    --threads $SLURM_CPUS_PER_TASK \
    tmp_vcf/${sample_name}.filtered.vcf.gz \
    > norm_vcf2/${sample_name}.norm.vcf.gz

# 索引结果文件
bcftools index norm_vcf2/${sample_name}.norm.vcf.gz

