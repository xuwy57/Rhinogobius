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

def load_gene_to_longest_transcript(gff_file):
    """
    从最长转录本 GFF3 文件读取 gene_id -> longest_transcript_id 映射
    """
    gene_to_longest = {}
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
                    gene_to_longest[gene_id] = transcript_id
    print(f"从 {gff_file} 中读取到 {len(gene_to_longest)} 个基因→最长转录本映射")
    return gene_to_longest

def process_genes_single_counting(variant_file, species_file, gff_file, gff_longest_file, output_prefix, target_clades):
    """
    处理基因数据，统计各个clade和lifestyle组中含frameshift_variant突变且基因型为2的物种数
    （每个物种在每个基因中最多计一次），最后输出基因对应的最长转录本 ID
    """
    # === Step 1. 载入 GFF 映射 ===
    transcript_to_gene = load_transcript_to_gene_map(gff_file)
    gene_to_longest = load_gene_to_longest_transcript(gff_longest_file)

    # === Step 2. 读取物种信息 ===
    species_df = pd.read_csv(species_file, delim_whitespace=True)
    print(f"读取到 {len(species_df)} 个物种信息")
    print("物种信息文件的列名:", species_df.columns.tolist())

    if 'index' not in species_df.columns:
        print("错误: 物种信息文件中缺少 'index' 列")
        return

    clade_col = 'clade' if 'clade' in species_df.columns else 'Clade'
    lifestyle_col = None
    for col in species_df.columns:
        if col.lower() == 'lifestyle':
            lifestyle_col = col
            break
    if lifestyle_col is None:
        print("错误: 未找到 lifestyle 列")
        return

    amphidromous_marker = None
    landlocked_marker = None
    for value in species_df[lifestyle_col].unique():
        if 'amphidromous' in str(value).lower():
            amphidromous_marker = value
        if 'landlocked' in str(value).lower():
            landlocked_marker = value

    species_df['index'] = species_df['index'] - 2
    n_fixed_cols = 8
    with open(variant_file, 'r') as f:
        first_line = f.readline().strip().split('\t')
    total_columns = len(first_line)
    n_species_cols = total_columns - n_fixed_cols

    fixed_cols = ['chr', 'pos', 'ref', 'alt', 'transcript_id', 'effect', 'impact', 'aa_change']
    all_cols = fixed_cols + [f'species_{i+1}' for i in range(n_species_cols)]

    group_indices = {}
    for clade in target_clades:
        clade_species = species_df[species_df[clade_col] == clade]
        if not clade_species.empty:
            group_indices[clade] = clade_species['index'].tolist()
    if amphidromous_marker:
        amphidromous_species = species_df[species_df[lifestyle_col] == amphidromous_marker]
        if not amphidromous_species.empty:
            group_indices['amphidromous'] = amphidromous_species['index'].tolist()
    if landlocked_marker:
        landlocked_species = species_df[species_df[lifestyle_col] == landlocked_marker]
        if not landlocked_species.empty:
            group_indices['landlocked'] = landlocked_species['index'].tolist()

    gene_stats = defaultdict(lambda: defaultdict(set))

    chunksize = 10000
    for chunk in pd.read_csv(variant_file, sep='\t', header=None, names=all_cols, chunksize=chunksize):
        frameshift_variants = chunk[chunk['effect'] == 'frameshift_variant'].copy()
        if frameshift_variants.empty:
            continue

        for _, row in frameshift_variants.iterrows():
            transcript_id = str(row['transcript_id'])
            gene_id = transcript_to_gene.get(transcript_id)
            if gene_id is None:
                continue

            for group_name, indices in group_indices.items():
                for idx in indices:
                    col_idx = idx - 1
                    if col_idx < 0 or col_idx >= len(all_cols):
                        continue
                    col_name = all_cols[col_idx]
                    if row[col_name] == 2:
                        gene_stats[gene_id][f'{group_name}_species'].add(idx)

    results = []
    for gene_id, stats in gene_stats.items():
        longest_tx = gene_to_longest.get(gene_id, gene_id)
        row = {'gene_id': longest_tx}
        for group_name in group_indices.keys():
            species_set = stats.get(f'{group_name}_species', set())
            row[f'{group_name}_count'] = len(species_set)
            row[f'{group_name}_total'] = len(group_indices[group_name])
        results.append(row)

    output_file = f"{output_prefix}_gene_stats.txt"
    if results:
        result_df = pd.DataFrame(results)
        result_df.to_csv(output_file, sep='\t', index=False)
        print(f"找到 {len(result_df)} 个基因的统计信息，结果已保存到 {output_file}")
    else:
        print("没有找到符合条件的基因")
        pd.DataFrame(columns=['gene_id'] + [f"{g}_count" for g in group_indices.keys()] + [f"{g}_total" for g in group_indices.keys()]).to_csv(output_file, sep='\t', index=False)

def main():
    if len(sys.argv) != 6:
        print("用法: python single_counting_frameshift.py <变异文件> <物种信息文件> <gff文件> <最长转录本gff文件> <输出前缀>")
        sys.exit(1)

    variant_file = sys.argv[1]
    species_file = sys.argv[2]
    gff_file = sys.argv[3]
    gff_longest_file = sys.argv[4]
    output_prefix = sys.argv[5]

    target_clades = ['clade2', 'clade5', 'clade7']
    process_genes_single_counting(variant_file, species_file, gff_file, gff_longest_file, output_prefix, target_clades)

if __name__ == "__main__":
    main()
