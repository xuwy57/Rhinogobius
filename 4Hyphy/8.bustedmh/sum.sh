cat Relax_LTY/* | awk -F'\t' '$9 < 0.05' | grep 'Success' > sumLTY_p0.05.txt
