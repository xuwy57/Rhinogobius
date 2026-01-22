#!/bin/bash
#SBATCH -p blade
#SBATCH -c 52
#SBATCH --mem=100G

module load bcftools
genome=/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/ref/ref/Rform.fasta

#bcftools view -O z -M 2 -m 2 -i 'TYPE=="snp" || TYPE=="indel"' outgroup_vcf/out.Carassius_auratus.Cyprinus_carpio_carpio.vcf.gz | \
#	bcftools norm -w 99999999 -O z -f $genome -  > outgroup_vcf/out.Carassius_auratus.Cyprinus_carpio_carpio.snps.indels.vcf.gz &
#bcftools index outgroup_vcf/out.Carassius_auratus.Cyprinus_carpio_carpio.snps.indels.vcf.gz

#bcftools view -O z -M 2 -m 2 -i 'TYPE=="snp" || TYPE=="indel"' outgroup_vcf/out.Carassius_auratus.Sinocyclocheilus_anshuiensis.vcf.gz | bcftools norm -w 99999999 -O z -f $genome -  > outgroup_vcf/out.Carassius_auratus.Sinocyclocheilus_anshuiensis.snps.indels.vcf.gz &
#bcftools index outgroup_vcf/out.Carassius_auratus.Sinocyclocheilus_anshuiensis.snps.indels.vcf.gz

for i in ./*.g.vcf.gz; do
	sBase=`basename $i`
	sBase=${sBase/%.g.vcf.gz/}
	bcftools view -O z -M 2 -m 2 -i 'TYPE=="snp" || TYPE=="indel"' $i | bcftools norm -w 99999999 -O z -f $genome - > norm_vcf/$sBase.g.snps.indels.vcf.gz &
#	bcftools index norm_vcf/$sBase.g.snps.indels.vcf.gz &
done

for i in ./gvcf/*.g.vcf.gz; do
	sBase=`basename $i`
	sBase=${sBase/%.g.vcf.gz/}
	bcftools view -O z -M 2 -m 2 -i 'TYPE=="snp" || TYPE=="indel"' $i | bcftools norm -w 99999999 -O z -f $genome - > norm_vcf/$sBase.g.snps.indels.vcf.gz &
#	bcftools index norm_vcf/$sBase.g.snps.indels.vcf.gz &
done

wait

echo "All tasks completed"  # 只有所有任务完成后才会执行
