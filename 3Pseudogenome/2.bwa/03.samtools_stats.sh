for r1 in ./*.bam.dedup.bam
do
	sp=`basename $r1`;
	sp=${sp/%.bam.dedup.bam}
	samtools flagstat $r1 > "$sp".stats
done
wait;
