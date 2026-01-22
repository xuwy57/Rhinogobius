#!/bin/bash
#SBATCH --job-name=vcf_split
#SBATCH --output=logs/vcf_split_%A_%a.out
#SBATCH --error=logs/vcf_split_%A_%a.err
#SBATCH --array=1-52
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00

set -euo pipefail
module load samtools

INPUT_DIR="norm_vcf"
OUTPUT_DIR="output_split2"

# 从样本列表中读取当前任务对应的文件名
SAMPLE_LIST="samples.txt"

# 获取当前任务对应的文件路径
vcf_path=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
sample_name=$(basename "$vcf_path" .vcf.gz)
sample_dir="$OUTPUT_DIR/$sample_name"
mkdir -p "$sample_dir"

echo "Processing $sample_name ..."

# 先保存表头
zcat "$vcf_path" 2>/dev/null | awk '/^#/ {print > hdr; next}
    { if(NF<10) exit; print > tmp }' hdr="$sample_dir/${sample_name}.header" tmp="$sample_dir/${sample_name}.tmp" || true

# 如果有数据文件则分染色体
if [[ -s "$sample_dir/${sample_name}.tmp" ]]; then
    awk -v dir="$sample_dir" '{print > dir"/"$1".vcf"}' "$sample_dir/${sample_name}.tmp"
    rm "$sample_dir/${sample_name}.tmp"
    echo "✅ Finished $sample_name"
else
    echo "⚠ $sample_name: No valid variant data found before corruption."
fi
