#!/bin/bash

# 输出文件
out_file="vcf_last_position.gz.tsv"

# 获取所有样本目录（只取目录名）
samples=($(ls -d norm_vcf_array/*/ | xargs -n1 basename))

# 获取所有出现过的 chromosome 名（去重并自然排序，去掉 .vcf.gz 后缀）
chroms=($(find norm_vcf_array -type f -name "*.vcf.gz" \
    | xargs -n1 basename \
    | sed 's/\.vcf\.gz$//' \
    | sort -Vu))  # -V 让数字按自然顺序排列

# 表头
echo -e "样本\t${chroms[*]}" | tr ' ' '\t' > "$out_file"

# 遍历样本
for s in "${samples[@]}"; do
    line="$s"
    for c in "${chroms[@]}"; do
        file="norm_vcf_array/$s/${c}.vcf.gz"
        if [[ -f "$file" ]]; then
            # 获取最后一行的 POS（忽略 # 开头的 header）
            pos=$(zgrep -v '^#' "$file" | tail -n 1 | cut -f2)
            line="$line\t${pos:-NA}"
        else
            line="$line\tNA"
        fi
    done
    echo -e "$line" >> "$out_file"
done
