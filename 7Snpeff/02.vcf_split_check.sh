#!/bin/bash

# 输出文件
out_file="vcf_last_position.tsv"

# 获取所有样本目录（只取目录名）
samples=($(ls -d output_split2/*/ | xargs -n1 basename))

# 获取所有出现过的 chromosome 名（去重并自然排序，去掉 .vcf 后缀）
chroms=($(find output_split2 -type f -name "*.vcf" \
    | xargs -n1 basename \
    | sed 's/\.vcf$//' \
    | sort -Vu))  # -V 让数字按自然顺序排列

# 表头
echo -e "样本\t${chroms[*]}" | tr ' ' '\t' > "$out_file"

# 遍历样本
for s in "${samples[@]}"; do
    line="$s"
    for c in "${chroms[@]}"; do
        file="output_split2/$s/${c}.vcf"
        if [[ -f "$file" ]]; then
            # 获取最后一行的 POS（忽略 # 开头的 header）
            pos=$(grep -v '^#' "$file" | tail -n 1 | cut -f2)
            line="$line\t${pos:-NA}"
        else
            line="$line\tNA"
        fi
    done
    echo -e "$line" >> "$out_file"
done

