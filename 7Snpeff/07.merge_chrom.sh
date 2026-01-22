#!/bin/bash
#SBATCH --job-name=merge_chrom
#SBATCH --output=logs/merge_%A_%a.out
#SBATCH --error=logs/merge_%A_%a.err
#SBATCH --mem=16G
#SBATCH --array=1-22

set -euo pipefail
module load bcftools

merge_dir="merged_by_chrom"
out_dir="final_merged"
mkdir -p "$out_dir"

# 获取任务编号对应的染色体目录
chrom_dir=$(ls -d ${merge_dir}/Chromosome* | sed -n "${SLURM_ARRAY_TASK_ID}p")
chrom=$(basename "$chrom_dir")

echo "[$(date)] Start merging $chrom"

out_vcf="${out_dir}/${chrom}.merged.vcf.gz"
extra_vcf1="minimap_outgroup_vcf/out.Rform.Tridentiger_barbatus.snps.indels.vcf.gz"
extra_vcf2="minimap_outgroup_vcf/out.Rform.Rhinogobius_giurinus.snps.indels.vcf.gz"

# 拆分大 VCF 的当前染色体
extra_chrom_vcf1="${out_dir}/${chrom}.extra1.vcf.gz"
bcftools view -r "$chrom" -Oz -o "$extra_chrom_vcf1" "$extra_vcf1"
bcftools index -t "$extra_chrom_vcf1"

extra_chrom_vcf2="${out_dir}/${chrom}.extra2.vcf.gz"
bcftools view -r "$chrom" -Oz -o "$extra_chrom_vcf2" "$extra_vcf2"
bcftools index -t "$extra_chrom_vcf2"

# 找到所有输入 VCF（该染色体目录下 + extra）
vcfs=(${chrom_dir}/*.vcf.gz "$extra_chrom_vcf1" "$extra_chrom_vcf2")

if [[ ${#vcfs[@]} -eq 0 ]]; then
    echo "[$(date)] No VCFs found for $chrom"
    exit 1
fi

# 合并
bcftools merge -Oz -o "$out_vcf" "${vcfs[@]}"

# 建索引
bcftools index -t "$out_vcf"

echo "[$(date)] Finished $chrom -> $out_vcf"
