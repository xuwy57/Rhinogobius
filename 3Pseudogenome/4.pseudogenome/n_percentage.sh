#!/bin/bash

# 定义FASTA文件名数组
#fasta_files=("file1.fasta" "file2.fasta" "file3.fasta")

# 循环遍历FASTA文件名数组
for file in ./*_pseudogenome.fasta
do
    # 获取fasta文件中所有碱基数
    total=`grep -v ">" $file | wc -m`

    # 获取fasta文件中所有N的数目
    nCount=`grep -v ">" $file | tr -cd 'N' | wc -m`

    # 计算N所占比例
    percentage=`echo "scale=2;($nCount/$total)*100" | bc -l`

    # 输出结果
    echo "在 $file 中所有N占总碱基数的比例为: $percentage%" >> n_count_916.log
wait;
done

