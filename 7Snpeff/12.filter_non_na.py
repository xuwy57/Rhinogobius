#!/usr/bin/env python3
import pandas as pd
import sys

if len(sys.argv) != 2:
    print("用法: python filter_non_na.py <变异文件>")
    sys.exit(1)

file_path = sys.argv[1]
df = pd.read_csv(file_path, sep='\t')

# 确定物种列的范围（假设前9列是固定信息，后面是物种列）
species_cols = df.columns[8:]

# 筛选不全为NA的行
non_na_df = df[~df[species_cols].isna().all(axis=1)]

# 统计结果
num_sites = len(non_na_df)
num_genes = non_na_df.iloc[:, 4].nunique()  # 假设第5列是gene_id

print(f'不全为NA的位点数: {num_sites}')
print(f'不全为NA的基因数: {num_genes}')
print(f'这些位点涉及{len(species_cols)}个物种')

# 保存结果
output_file = 'non_na_variants.txt'
non_na_df.to_csv(output_file, sep='\t', index=False)
print(f'结果已保存到 {output_file}')
