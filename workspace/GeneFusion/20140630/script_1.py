'''
Created on Jun 30, 2014

@author: HSH
'''

f = open('/Users/HSH/Roche/Data/expression/COAD_1598.gct')
lines = f.readlines()
print lines[2].split('\t')[2]
f.close()

f = open('/Users/HSH/Roche/Data/expression/COAD_1599.gct')
lines = f.readlines()
print lines[2].split('\t')[2]
f.close()

f = open('/Users/HSH/Roche/Data/expression/COAD_1865.gct')
lines = f.readlines()
print lines[2].split('\t')[2]
f.close()

# Test if genes are the same for multiple files of a certain cancer type
f = open('/Users/HSH/Roche/Data/expression/COAD_1598.gct')
f.readline()
f.readline()
f.readline()
genes_1598 = set([])

for line in f.readlines():
    g = line.split('\t')[0]
    genes_1598.add(g)

f.close()

f = open('/Users/HSH/Roche/Data/expression/COAD_1599.gct')
f.readline()
f.readline()
f.readline()
genes_1599 = set([])

for line in f.readlines():
    g = line.split('\t')[0]
    genes_1599.add(g)

f.close()

f = open('/Users/HSH/Roche/Data/expression/COAD_1865.gct')
f.readline()
f.readline()
f.readline()
genes_1865 = set([])

for line in f.readlines():
    g = line.split('\t')[0]
    genes_1865.add(g)

f.close()

print len(genes_1598),len(genes_1599),len(genes_1865)
print len(genes_1598 & genes_1599 & genes_1865)



########################################################################
import re
import rpy2.robjects as robjects

f = open('/Users/HSH/Roche/Data/expression/COAD_1598.gct')
header = f.readlines()[2].rstrip()
sample_fullname_list = header.split('\t')[2:]
sample_list = map(lambda x:re.match('TCGA-\w{2}-\w{4}-\w{3}-\w{3}-\w{4}-\w{2}',x).group(),sample_fullname_list)
sample_1598 = set(sample_list)
f.close()

f = open('/Users/HSH/Roche/Data/expression/COAD_1599.gct')
header = f.readlines()[2].rstrip()
sample_fullname_list = header.split('\t')[2:]
sample_list = map(lambda x:re.match('TCGA-\w{2}-\w{4}-\w{3}-\w{3}-\w{4}-\w{2}',x).group(),sample_fullname_list)
sample_1599 = set(sample_list)
f.close()

f = open('/Users/HSH/Roche/Data/expression/COAD_1865.gct')
header = f.readlines()[2].rstrip()
sample_fullname_list = header.split('\t')[2:]
sample_list = map(lambda x:re.match('TCGA-\w{2}-\w{4}-\w{3}-\w{3}-\w{4}-\w{2}',x).group(),sample_fullname_list)
sample_1865 = set(sample_list)
f.close()

samples = {'Sample_1598':list(sample_1598),'Sample_1599':list(sample_1599),'Sample_1865':list(sample_1865)}
robjects.globalenv['samples'] = robjects.ListVector(samples)
R = robjects.r
R_script = '''
library(VennDiagram)
venn.plot <- venn.diagram(samples,'~/Roche/20140630/Figure 1.png')
'''
R(R_script)
