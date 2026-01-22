mkdir -p bams

for i in /data2/projects/zma/Rhinogobius_spp/DNA/bam/*/*.bam.dedup.bam; do
	sBase=`basename $i`
	sStem=${sBase/.bam.dedup.bam/}
	mkdir -p bams/$sStem
	ln -sf `realpath $i`  bams/$sStem/sorted.bam
done
