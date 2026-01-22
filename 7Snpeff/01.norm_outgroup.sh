module load bcftools
genome=/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/ref/ref/Rform.fasta

#bcftools view -O z -M 2 -m 2 -i 'TYPE=="snp" || TYPE=="indel"' minimap_outgroup_vcf/out.Rform.Tridentiger_barbatus.vcf.gz > minimap_outgroup_vcf/filtered.vcf.gz
#bcftools norm -w 99999999 -O z -f $genome minimap_outgroup_vcf/filtered.vcf.gz  > minimap_outgroup_vcf/out..Rform.Tridentiger_barbatus.snps.indels.vcf.gz 
#bcftools index minimap_outgroup_vcf/out..Rform.Tridentiger_barbatus.snps.indels.vcf.gz

bcftools view -O z -M 2 -m 2 -i 'TYPE=="snp" || TYPE=="indel"' minimap_outgroup_vcf/out.Rform.Rhinogobius_giurinus.vcf.gz > minimap_outgroup_vcf/Rhinogobius_giurinus.filtered.vcf.gz
bcftools norm -w 99999999 -O z -f $genome minimap_outgroup_vcf/Rhinogobius_giurinus.filtered.vcf.gz  > minimap_outgroup_vcf/out.Rform.Rhinogobius_giurinus.snps.indels.vcf.gz
bcftools index minimap_outgroup_vcf/out.Rform.Rhinogobius_giurinus.snps.indels.vcf.gz

