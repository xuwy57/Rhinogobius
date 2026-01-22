trimjar=/data/software/Trimmomatic-0.39/trimmomatic-0.39.jar
adapterfa=/data/software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa


sSample=Rsimilis
sR1=./SRR15013109.split/SRR15013109_pass_1.fastq.gz
sR2=./SRR15013109.split/SRR15013109_pass_2.fastq.gz
java -jar $trimjar PE -threads 12 -phred33 $sR1 $sR2 $sSample.paired_1.fq.gz $sSample.unpaired_1.fq.gz  $sSample.paired_2.fq.gz $sSample.unpaired_2.fq.gz  ILLUMINACLIP:$adapterfa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > trim.$sSample.log 2>&1 &

