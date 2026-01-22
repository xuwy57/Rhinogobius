#!/bin/bash

# 定义源目录和目标目录
src_dir="all_pep"
dest_dir="all_filter_pep"

# 确保目标目录存在
mkdir -p "$dest_dir"

# 遍历源目录下的所有文件
for file in "$src_dir"/*; do
    # 检查是否是文件（不是目录）
    if [ -f "$file" ]; then
        # 统计文件中 '>' 的数量
        count=$(grep -o '>' "$file" | wc -l)
        
        # 如果恰好有 58 个 '>'
        if [ "$count" -eq 58 ]; then
            # 移动文件到目标目录
            mv "$file" "$dest_dir/"
            echo "已移动: $file"
        fi
    fi
done

echo "筛选完成！"
