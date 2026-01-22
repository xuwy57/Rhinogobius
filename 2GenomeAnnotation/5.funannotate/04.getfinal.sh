sp=btpf0
genome=trained/predict_misc/genome.softmasked.fa
tooldir=/public/apps/miniconda3/envs/funannotate/opt/pasa-2.5.2/misc_utilities/
gff=`ls trained/update_results/Bahaha_taipingensis_HIFIASM*.gff3`
grep -vP "^#" $gff | \
sed 's/FUN_/'$sp'/g' | \
sed 's/novel_/'$sp'_/g' \
> $sp.full.gff3



grep -P "\ttRNA\t" $gff > $sp.tRNA.gff3

$tooldir/gff3_file_single_longest_isoform.pl $sp.full.gff3 > $sp.longest_isoform.gff3

$tooldir/gff3_file_to_proteins.pl $sp.full.gff3 $genome cDNA > $sp.full.mrna.fa
$tooldir/gff3_file_to_proteins.pl $sp.full.gff3 $genome prot > $sp.full.prot.fa
$tooldir/gff3_file_to_proteins.pl $sp.full.gff3 $genome CDS > $sp.full.cds.fa

$tooldir/gff3_file_to_proteins.pl $sp.longest_isoform.gff3 $genome cDNA > $sp.longest_isoform.mrna.fa
$tooldir/gff3_file_to_proteins.pl $sp.longest_isoform.gff3 $genome prot > $sp.longest_isoform.prot.fa
$tooldir/gff3_file_to_proteins.pl $sp.longest_isoform.gff3 $genome CDS > $sp.longest_isoform.cds.fa
