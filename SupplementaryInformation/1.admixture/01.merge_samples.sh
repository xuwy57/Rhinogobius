#bcftools merge  -o all.samples.vcf.gz -O  z --threads 64 /data2/projects/wxu/Rhinogobius_spp/psmc/all_gvcf/*.bam.genotyped.g.vcf.gz /data2/projects/wxu/Rhinogobius_spp/psmc/psmc_formosanus/gatk/Rformosanus.genotyped.g.vcf.gz
#tabix all.samples.vcf.gz
#bcftools query -l all.samples.vcf.gz > all_samples.txt
bcftools view -S keep_samples.txt -O z -o filtered_all.samples.vcf.gz all.samples.vcf.gz
