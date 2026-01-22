ref_gtf=../ref/Rform.gff3.with_geneid.gtf
featureCounts -a $ref_gtf -Q 1 -T 40 -p -t exon -o all.counts.txt /data2/projects/zma/Rhinogobius_spp/RNA/04_STAR/*A/*.out.bam
#for dir in ../STAR/*RNA; do
#	cd $dir
#	for bam in *bam; do
#		sample=${bam/%_gillRNAAligned.sortedByCoord.out.bam}
#		cd /data2/projects/wxu/apl_rnaseq/featurecount
#		featureCounts -T 1 -a $ref_gtf -t exon -g gene_id -o $sample.counts.txt $dir/$bam &
#	done
#done
cat all.counts.txt | cut -f1,7- > counts.txt
