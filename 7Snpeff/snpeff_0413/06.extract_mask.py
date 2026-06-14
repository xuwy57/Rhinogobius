#!/usr/bin/env python3
# extract_mask.py
import sys
import argparse
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description='从 pairwise.simp.ref 文件中提取 mask 区域（目标基因组上的可信区间）')
    parser.add_argument('-i', '--input', required=True, help='输入文件，列4=chr,5=start,6=end (1-based)')
    parser.add_argument('-o', '--output', default='mask_regions.bed', help='输出 BED 文件（0-based half-open）')
    parser.add_argument('--chr-col', type=int, default=5, help='...')
    parser.add_argument('--start-col', type=int, default=6, help='...')
    parser.add_argument('--end-col', type=int, default=7, help='...')
    return parser.parse_args()

def main():
    args = parse_args()
    mask_intervals = defaultdict(list)
    with open(args.input, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            if len(cols) < max(args.chr_col, args.start_col, args.end_col):
                continue
            chrom = cols[args.chr_col-1]
            start = int(cols[args.start_col-1]) - 1   # to 0‑based
            end   = int(cols[args.end_col-1])         # already 0‑based half‑open
            mask_intervals[chrom].append((start, end))

    with open(args.output, 'w') as out:
        for chrom in sorted(mask_intervals.keys()):
            for start, end in mask_intervals[chrom]:
                out.write(f"{chrom}\t{start}\t{end}\n")
    print(f"Mask 区域已保存至 {args.output}")

if __name__ == '__main__':
    main()
