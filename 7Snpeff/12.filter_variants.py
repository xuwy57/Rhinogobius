#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import sys
import os

def process_variants_for_clade(variant_file, species_file, output_prefix, target_clades):
    """
    处理变异数据，为指定的clade筛选符合条件的位点
    """
    
    # 读取物种信息 - 使用空格分隔
    species_df = pd.read_csv(species_file, delim_whitespace=True)
    print(f"读取到 {len(species_df)} 个物种信息")
    print("物种信息文件的列名:", species_df.columns.tolist())
    
    # 检查必要的列是否存在
    if 'index' not in species_df.columns:
        print("错误: 物种信息文件中缺少 'index' 列")
        return
        
    # 定义clade和lifestyle列（根据实际列名调整）
    clade_col = 'clade' if 'clade' in species_df.columns else 'Clade'
    
    # 查找lifestyle列（不区分大小写）
    lifestyle_col = None
    for col in species_df.columns:
        if col.lower() == 'lifestyle':
            lifestyle_col = col
            break
    
    if lifestyle_col is None:
        print("错误: 未找到 lifestyle 列")
        return
    
    print(f"使用 '{lifestyle_col}' 作为 lifestyle 列")
    
    # 打印所有lifestyle的值，帮助诊断
    print("所有lifestyle值:", species_df[lifestyle_col].unique())
    
    # 查找洄游物种的标记（不区分大小写）
    amphidromous_marker = None
    for value in species_df[lifestyle_col].unique():
        if 'amphidromous' in str(value).lower():
            amphidromous_marker = value
            break
    
    if amphidromous_marker is None:
        print("错误: 未找到 Amphidromous 标记")
        return
    
    print(f"使用 '{amphidromous_marker}' 作为洄游物种标记")
    
    # 读取变异文件的第一行来确定总列数
    with open(variant_file, 'r') as f:
        first_line = f.readline().strip().split('\t')
    
    n_fixed_cols = 8  # 前8列是固定信息
    total_columns = len(first_line)
    n_species_cols = total_columns - n_fixed_cols
    print(f"变异文件总列数: {total_columns}")
    print(f"物种列数量: {n_species_cols}")
    
    # 检查物种索引是否有效
    max_index = species_df['index'].max() - 1
    min_index = species_df['index'].min() - 1
    if max_index > total_columns or min_index < 1:
        print(f"错误: 物种索引必须在1到{total_columns}之间")
        print(f"当前索引范围: {min_index} 到 {max_index}")
        return
    
    # 创建列名列表
    fixed_cols = ['chr', 'pos', 'ref', 'alt', 'gene_id', 'effect', 'impact', 'aa_change']
    # 为所有列创建名称，包括固定列和物种列
    all_cols = fixed_cols + [f'species_{i+1}' for i in range(n_species_cols)]
    
    # 为每个目标clade创建输出
    for clade in target_clades:
        print(f"\n处理 {clade} 数据...")
        
        # 获取当前clade的物种
        clade_species = species_df[species_df[clade_col] == clade]
        if clade_species.empty:
            print(f"警告: 未找到 {clade} 的物种，跳过处理")
            continue
            
        clade_indices = clade_species['index'].tolist()
        amphidromous_species = species_df[species_df[lifestyle_col] == amphidromous_marker]['index'].tolist()
        
        print(f"{clade}物种数量: {len(clade_indices)}")
        print(f"洄游物种数量: {len(amphidromous_species)}")
        print(f"洄游物种索引: {amphidromous_species}")
        
        # 处理数据
        chunksize = 10000
        processed_chunks = []
        
        for chunk in pd.read_csv(variant_file, sep='\t', header=None, names=all_cols, chunksize=chunksize):
            # 筛选HIGH影响的变异，并创建副本以避免SettingWithCopyWarning
            high_impact = chunk[chunk['impact'] == 'HIGH'].copy()
            
            if high_impact.empty:
                continue
                
            # 初始化计数列
            high_impact[f'{clade}_2_count'] = 0
            high_impact[f'{clade}_NA_count'] = 0
            high_impact['amphidromous_0_count'] = 0
            high_impact['amphidromous_NA_count'] = 0
            high_impact[f'{clade}_total'] = len(clade_indices)
            high_impact['amphidromous_total'] = len(amphidromous_species)
            
            # 对每个物种列进行检查和计数
            for idx in clade_indices:
                # index直接对应列位置，从1开始
                col_idx = idx - 2  # 转换为0-based索引
                if col_idx < 0 or col_idx >= len(all_cols):
                    print(f"警告: 列索引 {idx} 超出范围 (0-{len(all_cols)-1})，跳过")
                    continue
                    
                col_name = all_cols[col_idx]
                
                # 统计当前clade物种
                high_impact.loc[high_impact[col_name] == 2, f'{clade}_2_count'] += 1
                high_impact.loc[(high_impact[col_name].isna()) | (high_impact[col_name] == 'NA'), f'{clade}_NA_count'] += 1
            
            for idx in amphidromous_species:
                # index直接对应列位置，从1开始
                col_idx = idx - 2  # 转换为0-based索引
                if col_idx < 0 or col_idx >= len(all_cols):
                    print(f"警告: 列索引 {idx} 超出范围 (0-{len(all_cols)-1})，跳过")
                    continue
                    
                col_name = all_cols[col_idx]
                
                # 统计洄游物种
                high_impact.loc[high_impact[col_name] == 0, 'amphidromous_0_count'] += 1
                high_impact.loc[(high_impact[col_name].isna()) | (high_impact[col_name] == 'NA'), 'amphidromous_NA_count'] += 1
            
            # 筛选符合条件的行
            condition1 = (high_impact[f'{clade}_2_count'] + high_impact[f'{clade}_NA_count']) == high_impact[f'{clade}_total']
            condition2 = high_impact[f'{clade}_2_count'] > 0
            clade_ok = condition1 & condition2
            condition3 = (high_impact['amphidromous_0_count'] + high_impact['amphidromous_NA_count']) == high_impact['amphidromous_total']
            amphidromous_ok = condition3
            filtered = high_impact[clade_ok & amphidromous_ok]
            
            processed_chunks.append(filtered)
        
        # 合并所有处理过的块
        if processed_chunks:
            result = pd.concat(processed_chunks, ignore_index=True)
            
            # 选择需要输出的列
            output_cols = [
                'chr', 'pos', 'ref', 'alt', 'gene_id', 
                f'{clade}_2_count', f'{clade}_NA_count', 
                'amphidromous_0_count', 'amphidromous_NA_count'
            ]
            result = result[output_cols]
            
            # 保存结果
            output_file = f"{output_prefix}_{clade}.txt"
            result.to_csv(output_file, sep='\t', index=False)
            print(f"找到 {len(result)} 个符合条件的位点，结果已保存到 {output_file}")
            
            # 显示一些统计信息
            if len(result) > 0:
                print(f"统计摘要 ({clade}):")
                print(f"  平均{clade}物种中值为2的个数: {result[f'{clade}_2_count'].mean():.2f}")
                print(f"  平均{clade}物种中值为NA的个数: {result[f'{clade}_NA_count'].mean():.2f}")
                print(f"  平均洄游物种中值为0的个数: {result['amphidromous_0_count'].mean():.2f}")
                print(f"  平均洄游物种中值为NA的个数: {result['amphidromous_NA_count'].mean():.2f}")
        else:
            print(f"没有找到符合条件的位点 ({clade})")
            # 创建空结果文件
            output_cols = [
                'chr', 'pos', 'ref', 'alt', 'gene_id', 
                f'{clade}_2_count', f'{clade}_NA_count', 
                'amphidromous_0_count', 'amphidromous_NA_count'
            ]
            output_file = f"{output_prefix}_{clade}.txt"
            pd.DataFrame(columns=output_cols).to_csv(output_file, sep='\t', index=False)

def main():
    if len(sys.argv) != 4:
        print("用法: python filter_variants.py <变异文件> <物种信息文件> <输出文件前缀>")
        print("示例: python filter_variants.py variants.txt species_info.txt results")
        sys.exit(1)
    
    variant_file = sys.argv[1]
    species_file = sys.argv[2]
    output_prefix = sys.argv[3]
    
    # 检查文件是否存在
    if not os.path.exists(variant_file):
        print(f"错误: 变异文件 '{variant_file}' 不存在")
        sys.exit(1)
        
    if not os.path.exists(species_file):
        print(f"错误: 物种信息文件 '{species_file}' 不存在")
        sys.exit(1)
    
    # 定义要处理的clade列表
    target_clades = ['clade2', 'clade5', 'clade7']
    
    # 处理数据
    process_variants_for_clade(variant_file, species_file, output_prefix, target_clades)

if __name__ == "__main__":
    main()
