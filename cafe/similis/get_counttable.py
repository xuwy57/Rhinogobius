from sys import argv
import os
from collections import defaultdict
import csv

script, count_dir, out_name = argv

out_f = open(out_name, 'w')
tsv_w = csv.writer(out_f, delimiter='\t')

inp_orthogroup = open('Orthogroups.tsv', 'r')

file_list = os.listdir(count_dir)
og_dic = defaultdict(list)
for file in file_list:
    if '.orthogroups_count.txt' in file:
        inp = open(file, 'r')
        for line in inp:
            line_split = line.split('\t')
            og_id = line_split[0].strip()
            og_dic[og_id].append(line_split[1].strip())
        inp.close()
for og in og_dic:
    tsv_w.writerow(og_dic[og])
out_f.close()


