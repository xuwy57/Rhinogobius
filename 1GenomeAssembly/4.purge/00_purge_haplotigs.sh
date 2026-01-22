#genome=/data/projects/zma/R.formosanus_data/07_nextpolish/01_rundir/genome.nextpolish.fasta
#reads=/data/projects/zma/R.formosanus_data/02_clean.data/ONT_BMK21110_clean.fq.gz
#sreads=/data/projects/zma/R.formosanus_data/02_clean.data/trimmomatic/DNA/R-for.paired_*.fq.gz
threads=20

#mkdir reads.fq sreads.fq

#ln -sf `realpath $genome` genome.fa
#ln -sf `realpath $reads` reads.fq
#ln -sf `realpath $sreads` sreads.fq

#minimap2 -t $threads -ax map-ont genome.fa reads.fq/*.fq.gz --secondary=no \
#	| samtools sort -m 40G -o laligned.bam -T tmp.lali

#minimap2 -t $threads -ax map-ont genome.fa sreads.fq/*.fq.gz --secondary=no \
#	| samtools sort -m 40G -o saligned.bam -T tmp.sali

#samtools index saligned.bam
#samtools index laligned.bam
samtools merge -@ $threads aligned.bam laligned.bam saligned.bam

/data/software/purge_haplotigs/bin/purge_haplotigs  hist  -b aligned.bam  -g genome.fa  -t $threads 

#exit;


#LOW=8
#MID=45
#HIGH=97

#/data/software/purge_haplotigs/bin/purge_haplotigs  cov   -i aligned.bam.gencov  -l $LOW  -m $MID  -h $HIGH  

#/data/software/purge_haplotigs/bin/purge_haplotigs  purge  -t $threads  -g genome.fa  -c coverage_stats.csv
