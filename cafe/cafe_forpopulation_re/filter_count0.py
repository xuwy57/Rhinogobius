out_f = open('filtered.txt', 'w')

inp_orthogroup = open('orthogroups.genecount.txt', 'r')

for line in inp_orthogroup:
    if line.startswith('(null)'):
        line_split = line.split('\t')
        for i in line_split[2:]:
            i = i.strip()
            if int(i) != 0:
                out_f.write(line)
                break
    else:
        out_f.write(line)
inp_orthogroup.close()
out_f.close()
