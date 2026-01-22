busco -c 24 -m genome -o buscoout --cpu 4\
    -i /data/projects/zma/R.formosanus_data/07_nextpolish/01_rundir/genome.nextpolish.fasta\
    -l /data/projects/zma/database/BUSCO/actinopterygii_odb10\
    --offline --download_path /data/projects/zma/database/BUSCO/\
    1>buscoout.log 2>&1
