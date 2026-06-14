module load samtools 
module load bcftools
module load clusterbasics

refid=rform
REF=/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/ref/ref/Rform.fasta
mkdir -p gvcf_expanded
sSampleID=Chaeturichtys_stigmatias
ALT=/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/snpEff/snpeff_0316/outgroup_gvcf/rform_Cstig/pseudogenome.Rform.Chaeturichtys_stigmatias.fa
minimapVCF=/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/snpEff/minimap_outgroup_vcf_0318/out.Rform.Chaeturichtys_stigmatias.vcf.gz
php 04.pseudogenomes2gvcf.php "$REF" "$ALT" "$ALT" "$sSampleID" "gvcf_expanded/$sSampleID.g.vcf.gz" "$minimapVCF" "$minimapVCF"
