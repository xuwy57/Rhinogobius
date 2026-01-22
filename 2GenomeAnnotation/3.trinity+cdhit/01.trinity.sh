/data/software/trinityrnaseq-v2.12.0/Trinity\
        --CPU 64 --seqType fq --max_memory 100G --samples_file ./sampleall.txt \
	--jaccard_clip --output trinity_rna > out.log 2>&1


#trinity拼接清后的fastq文件，样本多可以用制表符制作一个sample.txt 存放reads
