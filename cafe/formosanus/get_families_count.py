from sys import argv
from collections import defaultdict

script, depth_file, busco_file, out_name = argv

inp = open(depth_file, 'r')
inp_busco = open(busco_file, 'r')
inp_orthogroup = open('Orthogroups.tsv', 'r')
out_f = open(out_name, 'w')

rna_depth = defaultdict(list)

busco_dic = {}
for line in inp_busco:
    line_split = line.split('\t')
    if not line.startswith('#') and len(line_split) > 3:
        busco_dic[line_split[2]] = ''
inp_busco.close()

# attention, line5 reference specie
og_dic = defaultdict(list)
for line in inp_orthogroup:
    line_split = line.split('\t')
    if line.startswith('OG'):
        og_id = line_split[0]
        if line_split[5] != '':
            og_dic[og_id] = line_split[5].split(', ')
        else:
            og_dic[og_id] = []
inp_orthogroup.close()

rna_depth_raw = defaultdict(list)
for line in inp:
    line_split = line.split('\t')
    id = line_split[0].split('.')[0]
    depth = float(line_split[1].strip())
    rna_depth_raw[id].append(float(depth))
inp.close()

rna_depth = {}
for id in rna_depth_raw:
    av_depth = sum(rna_depth_raw[id]) / len(rna_depth_raw[id])
    rna_depth[id] = av_depth

busco_depth = 0
busco_num = 0
for busco in busco_dic:
    if busco not in rna_depth:
        continue
    busco_num += 1
    busco_depth += float(rna_depth[busco])
busco_depth = busco_depth / busco_num
print(out_name + '\tbusco depth:\t' + str(busco_depth) + '\tbusco num:\t'  + str(busco_num))

out_f.write('Orthogroup' + '\t' + depth_file.split('/')[-1].strip('.depth.txt') + '\t\n')
for og in og_dic:
    if og_dic[og] == []:
        out_f.write(og + '\t0\n')
    else:
        out_dic = {}
        family_count = 0
        for rna in og_dic[og]:
            if rna not in rna_depth:
                out_dic[rna] = 'NA'
            else:
                out_dic[rna] = str(rna_depth[rna] / busco_depth)
                family_count += rna_depth[rna] / busco_depth
        out_f.write(og + '\t' + str(round(family_count)) + '\t')
        n_print = 0
        for rna in out_dic:
            out_f.write(rna + ':' + str(out_dic[rna]))
            n_print += 1
            if n_print == len(out_dic):
                out_f.write('\n')
            else:
                out_f.write(',')
out_f.close()
