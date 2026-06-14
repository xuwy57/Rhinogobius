#!/bin/bash
#SBATCH -p blade,himem,gpu,hugemem
#SBATCH --job-name=vcf_norm_array
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --array=1-53               # 根据 samples.txt 的行数设置（样本数）
#SBATCH --cpus-per-task=1          # 每个样本任务只用 1 个 CPU（串行处理染色体）
#SBATCH --mem=8G                   # 每个样本任务内存，可根据单个染色体大小调整

set -euo pipefail
module load bcftools 

SPLIT_NORM_DIR="norm"
SAMPLE_LIST="samples.txt"
CHROM_LIST="chromosomes2.txt"
GENOME="/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/ref/ref/Rform.fasta"

# ==================== 改进后的样本读取部分 ====================
raw_line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST" | tr -d '\r' | xargs)

if [[ -z "$raw_line" ]]; then
    echo "【错误】samples.txt 第 ${SLURM_ARRAY_TASK_ID} 行为空或不存在！"
    exit 1
fi

sample_name=$(basename "$raw_line" | sed 's/\.[^.]*(gz|vcf|bcf)*$//')
vcf_path="$raw_line"

echo "[$(date)] 处理样本: $sample_name (任务 ${SLURM_ARRAY_TASK_ID})"
echo "   VCF 路径 → $vcf_path"
# ============================================================

sample_dir="${SPLIT_NORM_DIR}/${sample_name}"
mkdir -p "$sample_dir"

echo "[$(date)] 处理样本: $sample_name (任务 ${SLURM_ARRAY_TASK_ID})"

# 检查输入文件及索引
if [[ ! -f "$vcf_path" ]]; then
    echo "错误: 输入文件不存在 $vcf_path"
    exit 1
fi
if [[ ! -f "${vcf_path}.tbi" && ! -f "${vcf_path}.csi" ]]; then
    echo "警告: 缺少索引，尝试创建..."
    tabix -p vcf "$vcf_path" || { echo "无法创建索引，退出"; exit 1; }
fi

# 遍历所有染色体，逐个标准化
while read -r chrom; do
    output="${sample_dir}/${chrom}.norm.vcf.gz"
    
    if [[ -f "$output" ]]; then
        echo "  ✓ $chrom 已存在，跳过"
        continue
    fi
    
    echo " → 处理染色体 $chrom"

    # 执行 norm 并捕获退出码
    if bcftools norm \
        -r "$chrom" \
        -f "$GENOME" \
        -Oz \
        -w 10000 \
        --threads 1 \
        "$vcf_path" > "$output.tmp"; then

        # 成功后再改名 + 建索引，避免半成品文件
        mv "$output.tmp" "$output" && \
        tabix -p vcf "$output" && \
        echo " ✓ 完成 $chrom"
    else
        echo " ✗ 失败 $chrom (bcftools norm 退出码 $?)"
        echo "   → 请检查 .err 文件或参考基因组是否匹配"
        rm -f "$output.tmp" "$output"
    fi
done < "$CHROM_LIST"

echo "[$(date)] 样本 $sample_name 处理完成"
