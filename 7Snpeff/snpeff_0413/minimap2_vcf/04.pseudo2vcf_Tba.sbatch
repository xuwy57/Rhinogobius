module load samtools 
module load bcftools
module load clusterbasics

refid=rform
REF=/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/ref/ref/Rform.fasta
mkdir -p gvcf_expanded
sSampleID=Tridentiger_barbatus
ALT=/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/snpEff/snpeff_0316/outgroup_gvcf/rform_Tba/pseudogenome.Rform.Tridentiger_barbatus.fa
minimapVCF=/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/snpEff/minimap_outgroup_vcf_0318/out.Rform.Tridentiger_barbatus.vcf.gz
php 04.pseudogenomes2gvcf.php "$REF" "$ALT" "$ALT" "$sSampleID" "gvcf_expanded/$sSampleID.g.vcf.gz" "$minimapVCF" "$minimapVCF"

