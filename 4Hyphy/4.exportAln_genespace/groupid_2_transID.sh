cat genespace.orthogroups.txt | awk -F'\t' '{split($7, a, "|"); print $1"\t"a[3]}'  > groupid2transid.txt
