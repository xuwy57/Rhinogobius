repeatmasker_gff=/data/projects/zma/R.formosanus_data/10_funannotate/repeatmasker/hic/scf_unsorted.fa.out.gff
#grep 'Chr' /data/projects/zma/R.formosanus_data/10_funannotate/repeatmasker/hic/scf_unsorted.fa.out.gff > repeatmasker.gff
#cat /data2/projects/wxu/Rhinogobius_spp/cafe5/ref/Rform_longest_isoform.gff3 | awk -F '[\t;]' '/^Chr/{if($3=="CDS"){print $1,$4-1,$5,$7,$9}}' OFS="\t" | sed -e 's/ID=//' -e 's/gene_id=//' -e 's/transcript_id=//' -e 's/gene_name=//' | sort -k1,1 -k2,2n > cds.bed

bedtools intersect -a cds.bed -b repeatmasker.gff -s -v > intersect.bed
