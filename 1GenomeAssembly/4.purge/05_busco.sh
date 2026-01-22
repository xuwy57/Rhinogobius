busco -c 24 -m genome -o buscoout --cpu 8\
    -i 03_curated.fasta\
    -l /data/projects/zma/database/BUSCO/actinopterygii_odb10\
    --offline --download_path /data/projects/zma/database/BUSCO/\
    1>buscoout.log 2>&1
