
#!/bin/bash
# filter_last_two_have_data.sh

input_file="Whole_Genome_Annotated.vcf.gz"  
output_file="filtered_last_two_samples.vcf"

zcat -f "$input_file" | awk '
BEGIN {FS="\t"; OFS="\t"}
{
    if ($1 ~ /^#/) {
        print
        next
    }
    
    last_sample  = $(NF)     # 最后一个样本
    last2_sample = $(NF-1)   # 倒数第二个样本
    
    # 判断是否有数据：只要不是以 ./.:. 开头就认为有数据
    has_data_last  = (last_sample  !~ /^\.\/\.:\./)
    has_data_last2 = (last2_sample !~ /^\.\/\.:\./)
    
    # 只要最后两个样本中**至少有一个有数据**就保留
    if (has_data_last && has_data_last2) {
        print
    }
}' > "$output_file"

echo "筛选完成！"
echo "输出文件 → $output_file"
echo "共筛选出 $(grep -v "^#" "$output_file" | wc -l) 行有数据的记录"