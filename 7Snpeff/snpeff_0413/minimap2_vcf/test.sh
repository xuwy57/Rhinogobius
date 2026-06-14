REF=/public4/courses/ec3121/shareddata/Camellia_Sect_Chrysantha/partition_contigs_hap/purge/anc1/haphic/04.build/reviewed.chr.fa
ALT=/public4/courses/ec3121/shareddata/Camellia_Sect_Chrysantha/whole_genome_aln_cth/cth.vs.Camellia_chekiangoleosa/pseudogenome.cth.Camellia_chekiangoleosa.fa
minimapVCF=/public4/courses/ec3121/shareddata/Camellia_Sect_Chrysantha/whole_genome_aln_cth/cth.vs.Camellia_chekiangoleosa/out.cth.Camellia_chekiangoleosa.vcf.gz
php pseudogenomes2gvcf.php "$REF" "$ALT" "$ALT" Ccheck "testout.g.vcf.gz" "$minimapVCF" "$minimapVCF"
