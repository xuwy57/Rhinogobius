from sys import argv

script, paf = argv

inp_paf = open(paf, 'r')

querey_length = 0
target_length = 0
matching_bases = 0
align_length = 0

for line in inp_paf:
    line_split = line.split('\t')
    querey_length += int(line_split[1])
    target_length += int(line_split[6])
    matching_bases += int(line_split[9])
    align_length += int(line_split[10])

inp_paf.close()
print('querey_length = ' + str(querey_length))
print('target_length = ' + str(target_length))
print('matching_bases = '+ str(matching_bases))
print('align_length = ' + str(align_length))