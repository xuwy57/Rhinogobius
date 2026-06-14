#!/bin/bash

set -eo pipefail

# ============================================
# 配置（请根据实际情况修改）
# ============================================
FILTER_DIR="filtered_per_chr"        # 过滤后的染色体VCF目录
CHROM_LIST="chromosomes.txt"          # 染色体列表（顺序）
FINAL_VCF="final_merged.vcf.gz"       # 最终输出文件名

# ============================================
# 加载模块（根据集群环境调整）
# ============================================
module load bcftools

# ============================================
# 按顺序收集过滤后的VCF
# ============================================
if [[ ! -f "$CHROM_LIST" ]]; then
    echo "错误: 染色体列表文件不存在: $CHROM_LIST"
    exit 1
fi

concat_list=()
while read -r chrom; do
    vcf="${FILTER_DIR}/${chrom}.filtered.vcf.gz"
    if [[ -f "$vcf" ]]; then
        concat_list+=("$vcf")
    else
        echo "警告: 染色体 $chrom 的过滤文件不存在，跳过"
    fi
done < "$CHROM_LIST"

if [[ ${#concat_list[@]} -eq 0 ]]; then
    echo "错误: 没有找到任何过滤后的染色体文件"
    exit 1
fi

echo "找到 ${#concat_list[@]} 个染色体文件，开始合并..."

# ============================================
# 合并所有染色体
# ============================================
bcftools concat -Oz -o "$FINAL_VCF" "${concat_list[@]}"
bcftools index "$FINAL_VCF"

echo "========================================"
echo "合并完成！"
echo "最终文件: $FINAL_VCF"
echo "包含染色体数: ${#concat_list[@]}"
echo "========================================"
