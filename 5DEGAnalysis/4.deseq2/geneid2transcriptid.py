from sys import argv

script, file_gff, count = argv

inp_gff = open(file_gff, 'r')
inp_count = open(count, 'r')
out_name = count + '.renamed.txt'
out_f = open(out_name, 'w')


dic = {}

for line in inp_gff:
    if line.startswith('#') or "=" not in line:
        continue
    elif not line.startswith('#'):
        line_split = line.split('\t')
        feature = line_split[2]
        attributes = line_split[8]
        attributes_split = attributes.strip().split(';')
        if 'RNA' in feature:
            transcript = attributes_split[0]
            transcript_id = attributes_split[0].strip('ID=')
            gene_id = attributes_split[1].strip('Parent=')
            dic[gene_id] = transcript_id
        else:
            continue
inp_gff.close()

for line in inp_count:
    line_split = line.split('\t')
    if line_split[0] == 'Geneid':
        out_f.write(line)
    else:
        if line_split[0] in dic:
            out_f.write(dic[line_split[0]] + '\t' + '\t'.join(line_split[1:]))
        else:
            continue
inp_count.close()
out_f.close()