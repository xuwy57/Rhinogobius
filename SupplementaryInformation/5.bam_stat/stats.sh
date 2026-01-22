#!/usr/bin/bash
#SBATCH --job-name=bam
#SBATCH -p blade
#SBATCH -c 8
#SBATCH --mem=100G

module load samtools

# 创建结果表格文件（带时间戳）
output="coverage_summary_$(date +%Y%m%d).tsv"
echo -e "Sample\tAverage_Coverage" > "$output"

# 使用find定位所有BAM文件（包括子目录）
find /public3/group_crf/home/g21mazhy26/Rhinogobius_spp/DNA/GATK -name "*.bam" -print0 | while IFS= read -r -d $'\0' bam; do
    # 获取文件名作为样本ID
    sample=$(basename "$bam")
    sample="${sample%.bam}"

    # 检查BAM文件是否存在
    if [[ ! -f "$bam" ]]; then
        echo "错误: 文件不存在 - $bam" >&2
        continue
    fi

    # 检查BAM索引是否存在，不存在则创建
    if [[ ! -f "${bam}.bai" ]]; then
        echo "为 $sample 创建索引..."
        samtools index -@ 4 "$bam"
    fi

    # 并行计算覆盖度
    (
        # 计算平均覆盖度（修正awk语法）
        cov=$(samtools depth -@ 2 "$bam" 2>/dev/null | \
              awk '{sum+=$3; total++} END {if(total>0) printf "%.2f", sum/total; else printf "0.00"}')

        # 加锁写入结果
        flock -x 3
        echo -e "$sample\t$cov" >> "$output"
    ) 3>>"$output" &

    # 控制并行进程数（最多8个）
    if (( $(jobs -r -p | wc -l) >= 8 )); then wait -n; fi
done

# 等待所有后台任务完成
wait
echo "计算完成！结果保存在 $output"
