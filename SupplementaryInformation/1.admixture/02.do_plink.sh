plink=/data/software/Plink/plink
Name=Rhinogobius
#STEP 1
#bcftools view -h filtered_all.samples.vcf.gz > header.txt
#while read old_name new_name; do
#    sed -i "s#$old_name#$new_name#g" header.txt
#done < sample_mapping.txt
#{
#	  cat header.txt  # 插入你的自定义表头
#	   bcftools view filtered_all.samples.vcf.gz | grep -v '^#'  # 删除原注释和旧表头，保留数据行
#    } > merged.g.vcf
bgzip merged.g.vcf
tabix -p vcf merged.g.vcf.gz

bcftools view -O z -o filtered.snp.merged.vcf.gz --type "snps" -m2 -M2  -i 'F_PASS(GT!="mis")>0.8' merged.g.vcf.gz &&

#STEP2
# Generate the input file in plink format
$plink --vcf filtered.snp.merged.vcf.gz  --make-bed --out $Name --allow-extra-chr --const-fid 0 --geno 0.999 &&

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
awk '{$1="0";print $0}' $Name.bim > $Name.bim.tmp &&
mv $Name.bim.tmp $Name.bim
wait

