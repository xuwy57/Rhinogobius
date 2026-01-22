from sys import argv
from collections import defaultdict

script, cov, out_name= argv

inp = open(cov, 'r')
out_f = open(out_name, 'w')


cds_depth = []
rna_depth = defaultdict(list)
for line in inp:
    line_split = line.split('\t')
    chr = line_split[0]
    id = line_split[3].strip('cds.')
    start = line_split[1]
    end = line_split[2]
    raw_depth = int(line_split[5])
    rna_depth[id].append(raw_depth)
inp.close()


for rna in rna_depth:
    out_f.write(rna + '\t' + str(sum(rna_depth[rna]) / len(rna_depth[rna])) + '\n')
out_f.close()