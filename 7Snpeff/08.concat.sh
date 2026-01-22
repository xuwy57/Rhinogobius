#!/bin/bash
#SBATCH --job-name=concat
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

set -euo pipefail
module load bcftools
module load samtools
merge_dir="final_merged"

merged_list=$(ls "$merge_dir"/*.merged.vcf.gz | sort -V)
bcftools concat -Oz -o merged.snps.indels.vcf.gz $merged_list
bcftools index -t merged.snps.indels.vcf.gz

echo "✅ Final VCF: merged.snps.indels.vcf.gz"

