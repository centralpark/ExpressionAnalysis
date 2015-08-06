'''
Created on Jul 3, 2014
@author: HSH
'''

ifile_path='/Users/HSH/Roche/Data/expression/gene_annotation.txt'
f = open(ifile_path)
f.readline()
pseudogene_set = set([])
for line in f.readlines():
    line_split = line.rstrip().split('\t')
    if len(line_split) < 3:
        continue
    probe_id = line_split[0]
    symbol = line_split[1]
    title = line_split[2]
    if 'pseudogene' in title:
        pseudogene_set.add(probe_id)
f.close()

probe_set = set([])
f = open('/Users/HSH/Roche/Data/expression/THCA.gct')
f.readline()
map(lambda x: probe_set.add(x.split('\t')[0]),f.readlines())
f.close()

probe_set_2 = set([])
f = open('/Users/HSH/Roche/Data/expression/BRCA.gct')
f.readline()
map(lambda x: probe_set_2.add(x.split('\t')[0]),f.readlines())
f.close()

result_1 = probe_set & pseudogene_set
print 'Number of pseudogene probe used in THCA gene expression assay'
print str(len(result_1))

if probe_set == probe_set_2:
    print 'THCA and BRCA expression assays use the same set of probes'
else:
    print 'THCA and BRCA expression assays use the different set of probes'
