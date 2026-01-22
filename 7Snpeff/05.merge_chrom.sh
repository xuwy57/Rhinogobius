#!/bin/bash
#SBATCH --job-name=header_merge
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --array=1-22
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -euo pipefail
module load bcftools
module load samtools
# 目录设置
dir_raw="output_split2"
dir_norm="norm_vcf_array"
merge_dir="merged_by_chrom"
mkdir -p "$merge_dir"

# 当前染色体
chrom=$(sed -n "${SLURM_ARRAY_TASK_ID}p" chromosomes.txt)

files_to_merge=()

echo "[$(date)] Start processing $chrom"

# 处理 output_split2 文件：补 header + 压缩 + 索引
for s in "$dir_raw"/*; do
    raw_vcf="$s/${chrom}.vcf"
    header_file=$(ls "$s"/*header 2>/dev/null | head -n1)
    if [[ -f "$raw_vcf" && -f "$header_file" ]]; then
        tmp_vcf="$merge_dir/$(basename "$s")_${chrom}.vcf.gz"
        if [[ ! -f "$tmp_vcf" ]]; then
            echo "Adding header for $raw_vcf → $tmp_vcf"
            cat "$header_file" "$raw_vcf" | bgzip -c > "$tmp_vcf"
            bcftools index -t "$tmp_vcf"
        fi
        files_to_merge+=("$tmp_vcf")
    fi
done

# 处理 norm_vcf_array 文件：补 header + 压缩 + 索引
for s in "$dir_norm"/*; do
    norm_vcf="$s/${chrom}.norm.vcf.gz"
    sample_name=$(basename "$s")
    header_file=$(find "$dir_raw" -type f -name "${sample_name}*header" | head -n1)
    if [[ -f "$norm_vcf" && -f "$header_file" ]]; then
        tmp_vcf="$merge_dir/${sample_name}_${chrom}.norm.fixed.vcf.gz"
        if [[ ! -f "$tmp_vcf" ]]; then
            echo "Adding header for $norm_vcf → $tmp_vcf"
            (cat "$header_file"; zcat "$norm_vcf") | bgzip -c > "$tmp_vcf"
            bcftools index -t "$tmp_vcf"
        fi
        files_to_merge+=("$tmp_vcf")
    fi
done

# Merge 所有文件
if [[ ${#files_to_merge[@]} -gt 0 ]]; then
    echo "[$chrom] Merging ${#files_to_merge[@]} files..."
    bcftools merge -Oz -o "$merge_dir/${chrom}.merged.vcf.gz" "${files_to_merge[@]}"
    bcftools index -t "$merge_dir/${chrom}.merged.vcf.gz"
    echo "[$(date)] Finished merging $chrom"
else
    echo "[$chrom] No files to merge."
fi
