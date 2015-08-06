'''
Created on Jun 25, 2014
@author: HSH
'''

f = open('/Users/HSH/Roche/Data/all14.txt')
f.readline()

cancer_type = set([])
map(lambda x: cancer_type.add(x.split('\t')[0]),f.readlines())
cancer_type = [c for c in cancer_type]
f.close()

f = open('/Users/HSH/Roche/Data/all14.txt')
f.readline()
data = {}

for line in f.readlines():
    line_split = line.split('\t')
    numStudy = line_split[34]
    if not numStudy:
        continue
    cancer = line_split[0]
    gene5p = line_split[2]
    gene3p = line_split[3]
    fusion = gene5p+'/'+gene3p
    if fusion not in data:
        data[fusion] = {}
        for c in cancer_type:
            data[fusion][c] = '0'
    data[fusion][cancer] = line_split[6]

ofile = open('/Users/HSH/Roche/20140625/data_1_all_documented_fusion','w')
header = ''+'\t'+'\t'.join(cancer_type)
ofile.write(header+'\n')

for d in data:
    line_out = '\t'.join([d]+[data[d][c] for c in cancer_type])
    ofile.write(line_out+'\n')

f.close()
ofile.close()