CPU=32
REF=./ref_Rform.fasta
sSample=Rsimilis
sR1=Rsimilis_rename_pass_1.fastq.gz
sR2=Rsimilis_rename_pass_2.fastq.gz
bwa mem -t $CPU $REF $sR1 $sR2 \
        | samtools view  -u - | samtools sort - -m 20g -o $sSample.bam  > $sSample.log 2>&1 


