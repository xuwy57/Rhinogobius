sga=/data/software/sga/src/SGA/sga
sp=rfor
r1=/data/projects/shareddata/tra/trimmed/R-for.paired_1.fq.gz
r2=/data/projects/shareddata/tra/trimmed/R-for.paired_2.fq.gz
nSubSampleLn=600000000 #about 45x coverage
CPU=64


mkdir -p $sp

(
mkfifo $sp/r1.fastq
mkfifo $sp/r2.fastq

zcat -f $r1 | head -n $nSubSampleLn > $sp/r1.fastq &
zcat -f $r2 | head -n $nSubSampleLn > $sp/r2.fastq &


$sga preprocess --pe-mode 1  $sp/r1.fastq $sp/r2.fastq > $sp/$sp.fastq
$sga index -a ropebwt --no-reverse -t $CPU  $sp/$sp.fastq
mv $sp.* $sp/
$sga preqc -v -t $CPU  $sp/$sp.fastq >  $sp/$sp.preqc
#sga-preqc-report.py mygenome.preqc sga/src/examples/*.preqc

) >> $sp/log.txt 2>&1
