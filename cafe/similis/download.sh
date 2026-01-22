#prefetch -X 200G SRR15013109 

fastq-dump -v --outdir $sSample"_split"  --readids --read-filter pass --dumpbase --split-3 --clip --gzip --skip-technical SRR15013109
