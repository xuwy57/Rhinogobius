#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import sys
import os
from collections import defaultdict

def load_transcript_to_gene_map(gff_file):
    """
    从 GFF3 文件读取 transcript_id -> gene_id 映射
    """
    transcript_to_gene = {}
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            feature_type = parts[2].lower()
            attributes = parts[8]

            attr_dict = {}
            for item in attributes.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    attr_dict[k.strip()] = v.strip()

            if feature_type in ["mrna", "transcript"]:
                transcript_id = attr_dict.get("ID")
                gene_id = attr_dict.get("Parent")
                if transcript_id and gene_id:
                    transcript_to_gene[transcript_id] = gene_id
    print(f"从 {gff_file} 中读取到 {len(transcript_to_gene)} 个转录本→基因映射")
    return transcript_to_gene

def load_gene_to_longest_transcript(gff_longest_file):
    """
    从最长转录本 GFF3 文件读取 gene_id -> longest_transcript_id 映射
    """
    gene_to_longest = {}
    with open(gff_longest_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            feature_type = parts[2].lower()
            attributes = parts[8]

            attr_dict = {}
            for item in attributes.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    attr_dict[k.strip()] = v.strip()

            if feature_type in ["mrna", "transcript"]:
                transcript_id = attr_dict.get("ID")
                gene_id = attr_dict.get("Parent")
                if transcript_id and gene_id:
                    gene_to_longest[gene_id] = transcript_id
    print(f"从 {gff_longest_file} 中读取到 {len(gene_to_longest)} 个基因→最长转录本映射")
    return gene_to_longest

def parse_sample_column_name(col_name):
    """
    从列名中提取纯净样本名：
    - 如果是文件路径，取 basename 并移除 .bam 后缀
    - 否则直接返回原字符串（如 Rformosanus）
    """
    if '/' in col_name:
        name = os.path.basename(col_name)
        if name.endswith('.bam'):
            name = name[:-4]
        return name
    else:
        return col_name

def process_genes_per_species_counting(variant_file, species_file, gff_file, gff_longest_file, output_prefix):
    # 1. 加载映射文件
    transcript_to_gene = load_transcript_to_gene_map(gff_file)
    gene_to_longest = load_gene_to_longest_transcript(gff_longest_file)

    # 2. 读取物种信息文件，获取所有有效物种名（species 列）
    species_df = pd.read_csv(species_file, delim_whitespace=True)
    if 'species' not in species_df.columns:
        raise ValueError("物种信息文件必须包含 'species' 列")
    valid_species_set = set(species_df['species'].str.strip())

    # 3. 解析变异文件头，获取固定列名和样本列
    with open(variant_file, 'r') as f:
        header_line = f.readline().strip().split('\t')
    
    # 实际固定列名（与文件完全一致）
    fixed_columns = ['CHROM', 'POS', 'ANC_ALLELE', 'DER_ALLELE', 'RNA_ID', 'ANNOTATION', 'IMPACT', 'PROTEIN_VARIANT']
    # 检查文件头是否匹配固定列
    if header_line[:len(fixed_columns)] != fixed_columns:
        print("警告：文件头固定列与预期不完全匹配，将按照实际列索引处理。")
        print(f"预期前{len(fixed_columns)}列：{fixed_columns}")
        print(f"实际前{len(fixed_columns)}列：{header_line[:len(fixed_columns)]}")
    
    # 固定列数
    n_fixed = len(fixed_columns)
    raw_sample_cols = header_line[n_fixed:]   # 样本列原始名称
    
    # 解析纯净样本名，并记录 (纯净名, 原始列名, 列索引)
    sample_info = []
    for idx, raw_name in enumerate(raw_sample_cols, start=n_fixed):
        clean_name = parse_sample_column_name(raw_name)
        if clean_name in valid_species_set:
            sample_info.append((clean_name, raw_name, idx))
        else:
            print(f"警告: 样本 '{clean_name}' 不在物种信息文件中，将跳过此列")
    
    if not sample_info:
        sys.exit("错误: 没有找到任何匹配的样本列，请检查物种信息文件与变异文件头的一致性")
    
    ordered_sample_names = [info[0] for info in sample_info]
    sample_col_index = {info[0]: info[2] for info in sample_info}

    # 4. 初始化统计字典
    gene_species_counts = defaultdict(lambda: defaultdict(int))
    gene_species_seen = defaultdict(lambda: defaultdict(bool))
    missing_gene_log = []
    missing_longest_log = set()

    # 5. 分块读取变异文件，仅处理 IMPACT == 'HIGH' 的行
    chunksize = 10000
    # 注意：读取时使用 header=0，pandas 会自动将第一行作为列名
    for chunk in pd.read_csv(variant_file, sep='\t', header=0, low_memory=False, chunksize=chunksize):
        # 检查是否存在 IMPACT 列（注意大写）
        if 'IMPACT' not in chunk.columns:
            print("错误：文件中没有 'IMPACT' 列，请检查文件格式")
            break
        high_chunk = chunk[chunk['IMPACT'] == 'HIGH'].copy()
        if high_chunk.empty:
            continue
        
        for _, row in high_chunk.iterrows():
            transcript_id = str(row['RNA_ID'])   # 注意列名 RNA_ID
            gene_id = transcript_to_gene.get(transcript_id)
            if gene_id is None:
                missing_gene_log.append(transcript_id)
                continue
            
            # 遍历所有有效样本
            for sample_name in ordered_sample_names:
                col_idx = sample_col_index[sample_name]
                # 获取该样本的基因型值（通过 iloc 按位置取）
                if col_idx >= len(row):
                    continue
                val = str(row.iloc[col_idx]).strip()
                if val in ['nan', 'NA', '']:
                    continue
                gene_species_seen[gene_id][sample_name] = True
                if val in ['2', '2.0', 2]:
                    gene_species_counts[gene_id][sample_name] += 1

    # 6. 汇总结果
    results = []
    all_genes = set(gene_species_counts.keys()) | set(gene_species_seen.keys())
    for gene_id in all_genes:
        longest_tx = gene_to_longest.get(gene_id)
        if longest_tx is None:
            missing_longest_log.add(gene_id)
            longest_tx = gene_id
        
        row = {'gene_id': longest_tx}
        for sample_name in ordered_sample_names:
            if not gene_species_seen[gene_id].get(sample_name, False):
                row[sample_name] = "NA"
            else:
                row[sample_name] = gene_species_counts[gene_id].get(sample_name, 0)
        results.append(row)
    
    # 7. 输出结果
    output_file = f"{output_prefix}_gene_species_counts.txt"
    log_file = f"{output_prefix}_missing_ids.log"
    
    if results:
        result_df = pd.DataFrame(results)
        result_df = result_df[['gene_id'] + ordered_sample_names]
        result_df.to_csv(output_file, sep='\t', index=False)
        
        with open(log_file, 'w') as log:
            log.write(f"# 未找到对应 gene_id 的转录本 ({len(set(missing_gene_log))} 个):\n")
            for tid in sorted(set(missing_gene_log)):
                log.write(f"{tid}\n")
            log.write(f"\n# 未找到最长转录本的基因 ({len(missing_longest_log)} 个):\n")
            for gid in sorted(missing_longest_log):
                log.write(f"{gid}\n")
        
        print(f"找到 {len(result_df)} 个基因的统计信息，结果已保存到 {output_file}")
        print(f"日志文件已保存到 {log_file}")
    else:
        print("⚠️ 没有找到符合条件的位点-物种对")

def main():
    if len(sys.argv) != 6:
        print("用法: python gene_species_counting.py <变异文件> <物种信息文件> <gff文件> <最长转录本gff文件> <输出前缀>")
        sys.exit(1)
    
    variant_file = sys.argv[1]
    species_file = sys.argv[2]
    gff_file = sys.argv[3]
    gff_longest_file = sys.argv[4]
    output_prefix = sys.argv[5]
    
    process_genes_per_species_counting(variant_file, species_file, gff_file, gff_longest_file, output_prefix)

if __name__ == "__main__":
    main()