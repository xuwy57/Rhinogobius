#!/bin/bash
set -euo pipefail
module load bcftools
module load samtools

dir_raw="output_split2"
dir_norm="norm_vcf_array"
merge_dir="merged_by_chrom"

echo "[$(date)] Start fixing missing tbi..."

# 遍历每个染色体目录
for chrom_dir in "$merge_dir"/Chromosome*; do
    chrom=$(basename "$chrom_dir")

    echo "[$(date)] Checking $chrom_dir"

    for f in "$chrom_dir"/*.vcf.gz; do
        [[ -f "$f" ]] || continue
        if [[ ! -f "${f}.tbi" ]]; then
            echo "[$chrom] Missing index for $f"
            sample=$(basename "$f")
            sample=${sample%_Chromosome*.vcf.gz}
            echo "[$chrom] Sample = $sample"
            raw_vcf="$dir_raw/$sample/${chrom}.vcf"
            header_file=$(ls "$dir_raw/$sample"/*header 2>/dev/null | head -n1)
            if [[ -f "$raw_vcf" && -f "$header_file" ]]; then
              echo "Rebuilding raw file for $sample $chrom"
              cat "$header_file" "$raw_vcf" | bgzip -c > "$f"
              bcftools index -t "$f"
            fi
        fi
    done
done

echo "[$(date)] Finished fixing."
