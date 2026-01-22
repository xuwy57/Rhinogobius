#!/bin/bash
#SBATCH --job-name=vcf_norm
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=1-49  # 这里改成 tasks.txt 的行数

module load bcftools
module load htslib

TASK_LIST=tasks.txt
GENOME="/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/ref/ref/Rform.fasta"
TMP_VCF_DIR="tmp_vcf"
NORM_VCF2_DIR="norm_vcf_array"

# 读取任务
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASK_LIST")
sample=$(echo "$line" | awk '{print $1}')
chrom=$(echo "$line" | awk '{print $2}')

# 样本名处理（可按需调整）
sample_name="${sample#norm_vcf/}"
sample_name="${sample_name%.g.snps.indels.vcf.gz}"
vcf="${TMP_VCF_DIR}/${sample_name}.filtered.vcf.gz"

# 输出目录
out_dir="${NORM_VCF2_DIR}/${sample_name}"
out_vcf="${out_dir}/${chrom}.norm.vcf.gz"

echo "[$(date)] 处理样本 $sample_name 染色体 $chrom"

bcftools view -r "$chrom" "$vcf" \
	    | bcftools norm -O z -f "$GENOME" -w 99999999 \
	        > "$out_vcf"

tabix -p vcf "$out_vcf"

echo "[$(date)] 完成样本 $sample_name 染色体 $chrom"

