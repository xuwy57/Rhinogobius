# 1. 批量处理所有BAM文件
parallel -j 8 "mosdepth -t 4 -n -Q 1 {.}_cov {}" ::: *.bam

# 2. 提取关键指标并生成表格
echo -e "Sample\tMean_Cov\tSD_Cov\t10X_Coverage\t20X_Coverage" > coverage_summary.tsv

for f in *_cov.mosdepth.summary.txt; do
    sample=$(basename $f _cov.mosdepth.summary.txt)
    
    # 提取指标
    mean_cov=$(awk '/total/ {print $4}' $f)
    sd_cov=$(awk '/total/ {print $5}' $f)
    cov_10x=$(awk -F'\t' '/10X/ {print $2*100}' ${f%.summary.txt}.mosdepth.global.dist.txt)
    cov_20x=$(awk -F'\t' '/20X/ {print $2*100}' ${f%.summary.txt}.mosdepth.global.dist.txt)
    
    # 写入表格
    echo -e "$sample\t$mean_cov\t$sd_cov\t$cov_10x\t$cov_20x" >> coverage_summary.tsv
done
