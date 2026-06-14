#!/usr/bin/bash
#SBATCH -p blade,himem,gpu,hugemem
#SBATCH -c 1
#SBATCH --mem=40G
#SBATCH -J php_getSV
#SBATCH -o php_%j.log
#SBATCH -e php_%j.err

module load clusterbasics

# 2. 检查输入文件是否存在
VCF_FILE="ann.final_merged.vcf.gz"
CLOSE_OUTGROUP="/data/projects/zma/Rhinogobius_spp/DNA/bam/Z6_R-giurinus_Zhejianglongquanshizhulongzhen.bam"
FAR_OUTGROUP="Tridentiger_barbatus"

if [ ! -f "$VCF_FILE" ]; then
    echo "错误：未找到输入文件 $VCF_FILE"
    exit 1
fi

# 检查掩码文件（PHP 脚本中写死的命名格式）
MASK1="pairwise.simp.ref.${CLOSE_OUTGROUP}.txt"
MASK2="pairwise.simp.ref.${FAR_OUTGROUP}.txt"


if [ ! -f "$MASK1" ] || [ ! -f "$MASK2" ]; then
    echo "警告：掩码文件 $MASK1 或 $MASK2 不存在，脚本可能会报错或不应用过滤。"
fi

# 3. 运行 PHP 脚本
echo "开始运行 PHP 变异提取脚本..."
echo "开始时间: $(date)"

# 使用 php 命令运行脚本
# 因为 PHP 脚本内部使用了 zcat，所以不需要预先解压
php -d memory_limit=30G 05.getSV.RT.php

echo "处理完成时间: $(date)"

# 4. 检查输出结果
OUT_FILE="SV.polarizedRT.out.txt"
if [ -f "$OUT_FILE" ]; then
    echo "成功生成输出文件: $OUT_FILE"
    echo "总行数: $(wc -l < $OUT_FILE)"
else
    echo "处理可能失败，未找到输出文件。"
fi
