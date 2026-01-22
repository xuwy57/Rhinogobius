from sys import argv
from collections import defaultdict

script, file_gff = argv

inp_f_gff = open(file_gff, 'r')
out_name = file_gff + '.with_geneid.gtf'
out_f = open(out_name, 'w')

for line in inp_f_gff:
    if line.startswith('#'):
        out_f.write(line)
    elif not line.startswith('#'):
        line_split = line.split('\t')
        feature = line_split[2]
        attributes = line_split[8]
        attributes_split = attributes.strip().split(';')
        if feature == 'gene':
            locus = attributes_split[0]
            locus_id = attributes_split[0].strip('ID=')
            gene_id = 'gene_id "' + locus_id + '"'
        elif 'RNA' in feature:
            transcript = attributes_split[0]
            transcript_id = attributes_split[0].strip('ID=')
            transcript_id = 'transcript_id "' + transcript_id + '"'
            out_f.write('\t'.join(line_split[0:8]) + '\t')
            out_f.write(gene_id + '; ' + transcript_id + ';\n')
        elif feature == 'CDS' or feature == 'exon':
            out_f.write('\t'.join(line_split[0:8]) + '\t')
            out_f.write(gene_id + '; ' + transcript_id + ';\n')
        else:
            continue
inp_f_gff.close()
out_f.close()