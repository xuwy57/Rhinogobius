#!/bin/bash
#SBATCH -p blade,himem,gpu,hugemem
#SBATCH --job-name=merge_filter_chr
#SBATCH --output=logs/merge_filter_%A_%a.out
#SBATCH --error=logs/merge_filter_%A_%a.err
#SBATCH --array=1-22
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00

set -eo pipefail

# ============================================
# 配置区域
# ============================================
INPUT_BASE="norm"                    # 每个样本的norm VCF目录
SAMPLE_LIST="samples.txt"            # 样本列表文件
CHROM_LIST="chromosomes.txt"         # 染色体列表文件
FILTER_DIR="filtered_per_chr"        # 输出目录

# ============================================
# 加载模块
# ============================================
module load bcftools

# ============================================
# 获取当前任务对应的染色体
# ============================================
if [[ ! -f "$CHROM_LIST" ]]; then
    echo "错误: 染色体列表文件不存在: $CHROM_LIST"
    exit 1
fi

chrom=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CHROM_LIST" | tr -d '\r ' | xargs)
if [[ -z "$chrom" ]]; then
    echo "警告: 任务 ${SLURM_ARRAY_TASK_ID} 无对应染色体，跳过"
    exit 0
fi

echo "[$(date)] 开始处理染色体: $chrom (任务 ${SLURM_ARRAY_TASK_ID})"

# ============================================
# 读取样本列表并收集VCF文件（严格使用你提供的正则）
# ============================================
if [[ ! -f "$SAMPLE_LIST" ]]; then
    echo "错误: 样本列表文件不存在: $SAMPLE_LIST"
    exit 1
fi

vcf_list=()

while IFS= read -r raw_line || [[ -n "$raw_line" ]]; do
    # === 严格使用你提供的表达式获取 sample name ===
    sample_name=$(basename "$raw_line" | sed 's/\.[^.]*(gz|vcf|bcf)*$//')
    
    if [[ -z "$sample_name" ]]; then
        echo "  警告: 无法提取样本名，跳过行: $raw_line"
        continue
    fi

    vcf="${INPUT_BASE}/${sample_name}/${chrom}.norm.vcf.gz"
    
    if [[ -f "$vcf" ]]; then
        vcf_list+=("$vcf")
        echo "  添加样本: ${sample_name}  (${vcf})"
    else
        echo "  警告: 文件不存在 -> ${sample_name} : ${vcf}"
    fi
done < <(sed 's/\r//g' "$SAMPLE_LIST")

if [[ ${#vcf_list[@]} -eq 0 ]]; then
    echo "错误: 染色体 ${chrom} 没有任何可用VCF文件"
    exit 1
fi

echo "  共收集到 ${#vcf_list[@]} 个样本的VCF，开始合并..."

# ============================================
# 执行合并 + 过滤
# ============================================
mkdir -p "$FILTER_DIR"
output="${FILTER_DIR}/${chrom}.filtered.vcf.gz"

if bcftools merge -m none -0 --threads $SLURM_CPUS_PER_TASK -Ou "${vcf_list[@]}" | \
   bcftools view -Ou -i 'TYPE=="snp" || TYPE=="indel"' --threads $SLURM_CPUS_PER_TASK | \
   bcftools view -Oz --threads $SLURM_CPUS_PER_TASK -o "$output"; then
    
    bcftools index --threads $SLURM_CPUS_PER_TASK "$output"
    echo "✓ 成功完成: ${output}"
    echo "  文件大小: $(du -h "$output" | cut -f1)"
else
    echo "✗ 合并或过滤失败: ${chrom}"
    rm -f "$output" "${output}.csi" "${output}.tbi"
    exit 1
fi

echo "[$(date)] 染色体 ${chrom} 处理完成"