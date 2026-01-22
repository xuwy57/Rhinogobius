# busco /public3/group_crf/home/g21mazhy26/R.formosanus_data/10_funannotate/funannotate/merge_gff/busco_prot/run_actinopterygii_odb10/full_table.tsv
# 第三列，空/rna id，e.g. FUN_018780-T1

# bam /data/projects/zma/Rhinogobius_spp/DNA/bam
# 软链接
# *.bam

# bedtools计算bam中基因覆盖度
# bedtools会保留mapQ低的read
# 把每个基因家族的拷贝数计算出来即可
# 记录busco覆盖度
# 以busco校准，记录多拷贝基因覆盖度

# step1 
#mkdir -p cov
#gff=/data2/projects/wxu/Rhinogobius_spp/ref/ref/Rform.longest_isoform.gff3
#cat $gff | awk -F '[\t;]' '/^Chr/{if($3=="CDS"){print $1,$4-1,$5,$9}}' OFS="\t" | sed -e 's/ID=//' -e 's/gene_id=//' -e 's/transcript_id=//' -e 's/gene_name=//' | sort -k1,1 -k2,2n > cds.bed

#for bam in ../bam/*.bam; do
#       base=`basename $bam`
#       pop=${base/%.bam}
#       samtools view -F 256 -b $bam | samtools view -F 2048 -b - | bedtools coverage -a cds.bed -b - -d -sorted > cov/$pop.cov &
#done
#wait


# step2
#mkdir -p depth
#for cov in cov/*.cov; do
#       base=`basename $cov`
#       pop=${base/%.cov}
#       python get_depth.py $cov depth/$pop.depth.txt &
#done
#wait

# step3
#get_families_count.py usage:
#script, depth_file, busco_file, out_name = argv
#attention, 需要改脚本里面参考物种的列位置
busco_file=full_table.tsv
for depth in depth/*.txt; do
	base=`basename $depth`
       pop=${base/%.depth.txt}
       python get_families_count.py $depth $busco_file $pop.orthogroups_count.txt &
done
wait

# step4
python get_counttable.py ./ R.pop.orthogroups.genecount.txt
IN=Orthogroups.GeneCount.tsv
# 下面根据需要的外群留下之前orthofinder结果列

cat $IN | cut -f 1,2,3,4,5,6,7,8,9 |  awk '{if ($1=="Orthogroup") {print "Desc\t"$0 } else {print "(null)\t"$0 } }' > outgroups.orthogroups.genecount.txt
paste outgroups.orthogroups.genecount.txt R.pop.orthogroups.genecount.txt > orthogroups.genecount.txt

#awk 'BEGIN {OFS="\t"} NR==1 {print; next} {sum=0; for (i=1; i<=NF; i++) sum+=$i; printf "%s\t%d\n", $0, sum}' tmp.genecount.txt > orthogroups.genecount.txt
#sed -i 's/\r//g' orthogroups.genecount.txt



