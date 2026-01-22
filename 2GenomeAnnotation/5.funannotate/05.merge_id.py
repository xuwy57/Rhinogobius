from sys import argv
from collections import defaultdict

script, SP, file_gff = argv

inp_f_gff = open(file_gff, 'r')
out_name = SP + '.merge_id.gff'
out_f = open(out_name, 'w')

dic_annotation = defaultdict(list)
dic_mRNA_parent = {}

locus = ''
gene = ''

for line in inp_f_gff:
    if line.startswith('#'):
        out_f.write(line)
    elif not line.startswith('#'):
        line_split = line.split('\t')
        feature = line_split[2]
        attributes = line_split[8]
        attributes_split = attributes.strip().split(';')
        if feature == 'locus':
            locus = attributes_split[0]
            locus_id = attributes_split[0].strip('ID=')
            out_f.write('\t'.join(line_split[0:8]) + '\t')
            out_f.write(locus + ';' + 'geneIDs=' + locus_id + ';' + attributes_split[2] + '\n')
            out_f.write('\t'.join(line_split[0:2]) + '\t' + 'gene' + '\t' + '\t'.join(line_split[3:8]) + '\t')
            out_f.write('ID=' + locus_id + ';' + 'locus=' + locus_id + '\n')
        elif feature == 'gene':
            continue
        elif feature == 'CDS' or feature == 'exon':
            out_f.write(line)
        else:
            out_f.write('\t'.join(line_split[0:8]) + '\t')
            parent_id = locus_id
            out_f.write(attributes_split[0] + ';' + 'Parent=' + parent_id + ';' + ';'.join(attributes_split[2:]) + '\n')
inp_f_gff.close()
out_f.close()