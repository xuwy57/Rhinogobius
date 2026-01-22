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

def process_genes_per_species_counting(variant_file, species_file, gff_file, gff_longest_file, output_prefix):
    transcript_to_gene = load_transcript_to_gene_map(gff_file)
    gene_to_longest = load_gene_to_longest_transcript(gff_longest_file)

    species_df = pd.read_csv(species_file, delim_whitespace=True)
    species_index_to_name = {row['index']: row['species'] for _, row in species_df.iterrows()}
    all_species_indices = species_df['index'].tolist()

    with open(variant_file, 'r') as f:
        first_line = f.readline().strip().split('\t')
    n_fixed_cols = 8
    n_species_cols = len(first_line) - n_fixed_cols
    fixed_cols = ['chr', 'pos', 'ref', 'alt', 'transcript_id', 'effect', 'impact', 'aa_change']
    all_cols = fixed_cols + [f'species_{i+1}' for i in range(n_species_cols)]

    # === 初始化 ===
    gene_species_counts = defaultdict(lambda: defaultdict(int))
    gene_species_seen = defaultdict(lambda: defaultdict(bool))  # 是否在此基因中出现过该物种
    missing_gene_log = []
    missing_longest_log = set()

    # === 统计 ===
    chunksize = 10000
    for chunk in pd.read_csv(variant_file, sep='\t', header=None, names=all_cols, chunksize=chunksize):
        high_impact = chunk[(chunk['impact'] == 'HIGH')].copy()
        if high_impact.empty:
            continue

        for _, row in high_impact.iterrows():
            transcript_id = str(row['transcript_id'])
            gene_id = transcript_to_gene.get(transcript_id)
            if gene_id is None:
                missing_gene_log.append(transcript_id)
                continue

            for species_index in all_species_indices:
                col_idx = n_fixed_cols + (species_index - 9)
                if col_idx >= len(all_cols):
                    continue
                col_name = all_cols[col_idx]
                val = str(row[col_name]).strip()
                species_name = species_index_to_name.get(species_index, f"species_{species_index}")

                if val == "nan" or val == "NA":
                    continue
                # 标记该物种在该基因出现过
                gene_species_seen[gene_id][species_name] = True
                # 若该位点有值且等于2，则计数
                if val in ["2", "2.0", 2]:
                    gene_species_counts[gene_id][species_name] += 1

    # === 汇总结果 ===
    results = []
    for gene_id in set(list(gene_species_counts.keys()) + list(gene_species_seen.keys())):
        longest_tx = gene_to_longest.get(gene_id)
        if longest_tx is None:
            missing_longest_log.add(gene_id)
            longest_tx = gene_id

        row = {'gene_id': longest_tx}
        for species_index in all_species_indices:
            species_name = species_index_to_name.get(species_index, f"species_{species_index}")
            # 若该物种在此基因从未出现过 → 输出 NA
            if not gene_species_seen[gene_id].get(species_name, False):
                row[species_name] = "NA"
            else:
                row[species_name] = gene_species_counts[gene_id].get(species_name, 0)
        results.append(row)

    output_file = f"{output_prefix}_gene_species_counts.txt"
    log_file = f"{output_prefix}_missing_ids.log"

    if results:
        result_df = pd.DataFrame(results)
        all_species_names = [species_index_to_name.get(idx, f"species_{idx}") for idx in all_species_indices]
        result_df = result_df[['gene_id'] + all_species_names]

        result_df.to_csv(output_file, sep='\t', index=False)

        with open(log_file, 'w') as log:
            log.write(f"# 未找到对应 gene_id 的转录本 ({len(set(missing_gene_log))} 个):\n")
            for tid in sorted(set(missing_gene_log)):
                log.write(f"{tid}\n")
            log.write(f"\n# 未找到最长转录本的基因 ({len(missing_longest_log)} 个):\n")
            for gid in sorted(missing_longest_log):
                log.write(f"{gid}\n")

        print(f" 找到 {len(result_df)} 个基因的统计信息，结果已保存到 {output_file}")
        print(f" 日志文件已保存到 {log_file}")
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
