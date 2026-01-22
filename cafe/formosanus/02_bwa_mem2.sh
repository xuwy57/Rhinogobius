CPU=32
REF=./ref_Rform.fasta
sSample=Rformosanus
sR1=/data2/projects/wxu/Rhinogobius_spp/cafe5/formosanus/DNA/R-for.paired_1.fq.gz
sR2=/data2/projects/wxu/Rhinogobius_spp/cafe5/formosanus/DNA/R-for.paired_2.fq.gz
bwa mem -t $CPU $REF $sR1 $sR2 \
        | samtools view  -u - | samtools sort - -m 20g -o $sSample.bam  > $sSample.log 2>&1 


