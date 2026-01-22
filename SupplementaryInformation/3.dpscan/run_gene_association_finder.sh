module load python3
python find_gene_associations.py \
  -r dp_scan/sign.1e-5.txt \
  -g Rform.gff3 \
  -o results.tsv -u
