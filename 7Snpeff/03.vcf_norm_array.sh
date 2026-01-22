#!/bin/bash
#SBATCH --job-name=vcf_norm
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

module load bcftools
module load htslib

SAMPLE_LIST=samples.txt
CHROM_LIST=chromosomes.txt
GENOME="/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/ref/ref/Rform.fasta"
TMP_VCF_DIR="tmp_vcf"
EXIST_VCF_BASEDIR="output_split2"
NORM_VCF2_DIR="norm_vcf_array"

mkdir -p "$NORM_VCF2_DIR"

# 样本和染色体数量
NUM_SAMPLES=$(wc -l < "$SAMPLE_LIST")
NUM_CHROMS=$(wc -l < "$CHROM_LIST")

# 任务 ID → 样本/染色体索引
sample_idx=$(( (SLURM_ARRAY_TASK_ID - 1) / NUM_CHROMS + 1 ))
chrom_idx=$(( (SLURM_ARRAY_TASK_ID - 1) % NUM_CHROMS + 1 ))

sample=$(sed -n "${sample_idx}p" "$SAMPLE_LIST")
chrom=$(sed -n "${chrom_idx}p" "$CHROM_LIST")

# 修改sample名字
sample_name="${sample#norm_vcf/}"
sample_out1_name="${sample_name%.vcf.gz}"
sample_name="${sample_name%.g.snps.indels.vcf.gz}"
vcf="${TMP_VCF_DIR}/${sample_name}.filtered.vcf.gz"
out1_dir="${EXIST_VCF_BASEDIR}/${sample_out1_name}"
out1_vcf="${out1_dir}/${chrom}.vcf"
out2_dir="${NORM_VCF2_DIR}/${sample_name}"
out2_vcf="${out2_dir}/${chrom}.norm.vcf.gz"

# 判断是否已存在norm结果
if [[ -f "$out1_vcf" ]]; then
    echo "✅ 已存在: $out1_vcf，跳过"
    exit 0
fi

mkdir -p "$out2_dir"

echo "[$(date)] 处理样本 $sample_name 染色体 $chrom"

bcftools view -r "$chrom" "$vcf" \
    | bcftools norm -O z -f "$GENOME" -w 99999999 \
    > "$out2_vcf"

tabix -p vcf "$out2_vcf"

echo "[$(date)] 完成样本 $sample_name 染色体 $chrom"
