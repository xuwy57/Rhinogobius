cd-hit-est\
	-i /data/projects/zma/R.formosanus_data/03_trinity.data/trinity_rna_all/Trinity.fasta\
       	-o R.formosanus_cdhit_all\
       	-c 0.95 -d 0 -n 10 -M 16000 -T 8 1>runcdhit_all.log 2>&1
