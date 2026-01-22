import pandas as pd
import numpy as np
import argparse
from collections import defaultdict

def find_nearest_features(region_start, region_end, feature_list, downstream_range=1000, search_upstream=False):
    """
    找到距离目标区域最近的基因组特征
    :param region_start: 目标区域起始位置
    :param region_end: 目标区域终止位置
    :param feature_list: 候选特征列表 [(start, end, feature_id, strand)]
    :param downstream_range: 下游搜索范围
    :param search_upstream: 是否搜索上游特征
    :return: (关联特征ID列表, 关联类型, 最小距离)
    """
    min_distance = float('inf')
    nearest_features = []
    
    for feat_start, feat_end, feat_id, strand in feature_list:
        # 计算重叠情况
        if region_end >= feat_start and region_start <= feat_end:
            distance = 0
            feature_type = "overlap"
        # 计算下游距离
        elif feat_start > region_end and (feat_start - region_end) <= downstream_range:
            distance = feat_start - region_end
            # 根据链方向确定上下游关系
            if strand == "+":
                feature_type = "downstream"
            else:  # 负链基因的下游在左侧
                feature_type = "upstream" if feat_start < region_start else "downstream"
        # 计算上游距离（如果启用）
        elif search_upstream and feat_end < region_start and (region_start - feat_end) <= downstream_range:
            distance = region_start - feat_end
            if strand == "+":
                feature_type = "upstream"
            else:  # 负链基因的上游在右侧
                feature_type = "downstream" if feat_end > region_end else "upstream"
        else:
            continue
        
        # 更新最近特征
        if distance < min_distance:
            min_distance = distance
            nearest_features = [(feat_id, feature_type, distance)]
        elif distance == min_distance:
            nearest_features.append((feat_id, feature_type, distance))
    
    return nearest_features

def main():
    # 设置命令行参数
    parser = argparse.ArgumentParser(description='查找显著区域及其最近mRNA')
    parser.add_argument('-r', '--regions', required=True, help='查询区域文件路径')
    parser.add_argument('-g', '--gff', required=True, help='GFF注释文件路径')
    parser.add_argument('-o', '--output', required=True, help='输出文件路径')
    parser.add_argument('-p', '--pvalue', type=float, default=0.05, 
                        help='P值显著性阈值 (默认: 0.05)')
    parser.add_argument('-d', '--downstream', type=int, default=1000, 
                        help='下游搜索范围(bp) (默认: 1000)')
    parser.add_argument('-u', '--upstream', action='store_true', 
                        help='是否搜索上游mRNA (默认: 仅下游)')
    parser.add_argument('--mrna-feature', default='mRNA', 
                        help='GFF中mRNA特征类型 (默认: mRNA)')
    parser.add_argument('--mrna-attribute', default='ID', 
                        help='GFF中mRNA ID属性名 (默认: ID)')
    parser.add_argument('--region-chr', type=int, default=0, 
                        help='区域文件中染色体列索引 (默认: 0)')
    parser.add_argument('--region-start', type=int, default=1, 
                        help='区域文件中起始位置列索引 (默认: 1)')
    parser.add_argument('--region-end', type=int, default=2, 
                        help='区域文件中结束位置列索引 (默认: 2)')
    parser.add_argument('--region-migratory', type=int, default=3, 
                        help='区域文件中洄游型深度列索引 (默认: 3)')
    parser.add_argument('--region-land', type=int, default=4, 
                        help='区域文件中陆峰型深度列索引 (默认: 4)')
    parser.add_argument('--region-pvalue', type=int, default=7, 
                        help='区域文件中P值列索引 (默认: 7)')
    
    args = parser.parse_args()
    
    # 读取查询区域文件
    try:
        region_cols = {
            'chr': args.region_chr,
            'start': args.region_start,
            'end': args.region_end,
            'migratory': args.region_migratory,
            'land': args.region_land,
            'p_value': args.region_pvalue
        }
        
        # 读取区域文件
        region_df = pd.read_csv(
            args.regions, 
            sep='\t', 
            header=None,
            usecols=list(region_cols.values()),
            names=list(region_cols.keys())
        )
        print(f"成功读取区域文件: {args.regions}, 共 {len(region_df)} 行")
    except Exception as e:
        print(f"读取区域文件错误: {e}")
        return
    
    # 筛选显著区域
    sig_df = region_df[
        (region_df['migratory'] > region_df['land']) & 
        (region_df['p_value'] < args.pvalue)
    ].copy()
    print(f"找到 {len(sig_df)} 个显著区域 (P < {args.pvalue})")
    
    # 读取GFF注释文件
    try:
        gff_df = pd.read_csv(
            args.gff, 
            sep='\t', 
            header=None,
            comment='#',
            names=['chr', 'source', 'feature', 'start', 'end', 
                   'score', 'strand', 'phase', 'attributes']
        )
        print(f"成功读取GFF文件: {args.gff}, 共 {len(gff_df)} 行")
    except Exception as e:
        print(f"读取GFF文件错误: {e}")
        return
    
    # 创建mRNA位置索引 {chrom: [(start, end, mrna_id, strand)]}
    mrna_index = defaultdict(list)
    mrna_count = 0
    
    for _, row in gff_df.iterrows():
        if str(row['feature']).lower() == args.mrna_feature.lower():
            # 从属性字段提取mRNA ID
            attrs = str(row['attributes'])
            mrna_id = None
            
            # 尝试提取ID
            for attr in attrs.split(';'):
                if f"{args.mrna_attribute}=" in attr:
                    mrna_id = attr.split(f"{args.mrna_attribute}=")[-1].strip()
                    break
                elif f"{args.mrna_attribute}:" in attr:
                    mrna_id = attr.split(f"{args.mrna_attribute}:")[-1].strip()
                    break
            
            if mrna_id:
                strand = row['strand'] if pd.notna(row['strand']) else '+'
                mrna_index[row['chr']].append((
                    int(row['start']), 
                    int(row['end']), 
                    mrna_id, 
                    strand
                ))
                mrna_count += 1
    
    print(f"索引到 {mrna_count} 个mRNA (特征类型: '{args.mrna_feature}', ID属性: '{args.mrna_attribute}')")
    
    # 为每个显著区域查找最近的关联mRNA
    results = []
    no_mrna_count = 0
    
    for idx, row in sig_df.iterrows():
        chrom, start, end = row['chr'], row['start'], row['end']
        nearest_mrnas = []
        min_distance = float('inf')
        assoc_type = "none"
        
        if chrom in mrna_index:
            # 查找重叠mRNA
            overlapping_mrnas = []
            for mrna in mrna_index[chrom]:
                if end >= mrna[0] and start <= mrna[1]:
                    overlapping_mrnas.append(mrna)
            
            # 如果有重叠mRNA
            if overlapping_mrnas:
                nearest = find_nearest_features(
                    start, end, 
                    overlapping_mrnas,
                    downstream_range=args.downstream,
                    search_upstream=args.upstream
                )
                if nearest:
                    nearest_mrnas = nearest
                    min_distance = 0
                    assoc_type = "overlap"
            
            # 如果没有重叠mRNA，查找下游mRNA
            if not nearest_mrnas:
                candidate_mrnas = []
                for mrna in mrna_index[chrom]:
                    # 检查是否在下游范围内
                    if (mrna[0] > end and (mrna[0] - end) <= args.downstream) or \
                       (args.upstream and mrna[1] < start and (start - mrna[1]) <= args.downstream):
                        candidate_mrnas.append(mrna)
                
                if candidate_mrnas:
                    nearest = find_nearest_features(
                        start, end, 
                        candidate_mrnas,
                        downstream_range=args.downstream,
                        search_upstream=args.upstream
                    )
                    if nearest:
                        nearest_mrnas = nearest
                        min_distance = min(g[2] for g in nearest)  # 取最小距离
                        assoc_type = nearest[0][1]  # 取第一个mRNA的类型
        
        # 处理结果
        mrna_ids = ";".join([g[0] for g in nearest_mrnas]) if nearest_mrnas else "None"
        mrna_types = ";".join(set(g[1] for g in nearest_mrnas)) if nearest_mrnas else "none"
        
        results.append({
            'chr': chrom,
            'start': start,
            'end': end,
            'migratory_depth': row['migratory'],
            'land_depth': row['land'],
            'p_value': row['p_value'],
            'associated_mrnas': mrna_ids,
            'association_type': mrna_types,
            'distance_to_mrna': min_distance if min_distance != float('inf') else -1
        })
        
        if not nearest_mrnas:
            no_mrna_count += 1
    
    # 创建结果DataFrame并保存
    result_df = pd.DataFrame(results)
    result_df.to_csv(args.output, sep='\t', index=False)
    
    # 输出统计信息
    print(f"\n处理完成! 结果已保存至 {args.output}")
    print(f"显著区域总数: {len(sig_df)}")
    print(f"找到关联mRNA的区域: {len(sig_df) - no_mrna_count}")
    print(f"未找到关联mRNA的区域: {no_mrna_count}")
    print("\n关联类型分布:")
    print(result_df['association_type'].value_counts())
    
    # 提取所有去重mRNA ID
    all_mrnas = []
    for mrna_list in result_df['associated_mrnas']:
        if mrna_list != "None":
            all_mrnas.extend(mrna_list.split(';'))
    
    unique_mrnas = sorted(set(all_mrnas))
    print(f"\n去重后唯一mRNA ID数量: {len(unique_mrnas)}")
    
    # 保存mRNA ID列表
    mrna_list_file = args.output.replace('.tsv', '_mrna_ids.txt')
    with open(mrna_list_file, 'w') as f:
        f.write("\n".join(unique_mrnas))
    print(f"唯一mRNA ID列表已保存至: {mrna_list_file}")

if __name__ == "__main__":
    main()