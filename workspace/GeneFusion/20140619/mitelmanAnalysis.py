'''
Created on Jun 19, 2014

@author: HSH
'''

import mysql.connector
import functools
import string
import re

cnx = mysql.connector.connect(user='root',database='MITELMAN')
cursor = cnx.cursor()

filename = '/Users/HSH/Roche/20140619/mitelman/molbiolclinassoc.dat'
f = open(filename,'r')
f.readline()
f.readline()
outFile = open('/Users/HSH/Roche/20140619/topMorphGene','w')

for line in f.readlines():
    line = line.split('\t')
    morph = line[3]
    top = line[4]
    genes = line[7]
    if not genes:
        continue
    genes = genes.split(',')
    map(lambda x: outFile.write('\t'.join([top,morph,x])+'\n'),genes)

f.close()
outFile.close()

f = open('/Users/HSH/Roche/20140619/topMorphGene')
lines = map(functools.partial(string.split,sep='\t'),f.readlines())
f.close()
cancer = filter(lambda x:not x[0].isspace(),lines)
cancer = map(lambda x : x[0]+'/'+x[1],cancer)
cancer = set(cancer)
gene_pair = set(map(lambda x: x[2].rstrip(),lines))

data = [[0]*len(cancer) for i in range(len(gene_pair))]
cancer_dict = {}
cancer_row = []
i = 0
for s in cancer:
    cancer_dict[s] = i
    cancer_row.append(s)
    i = i+1
gene_pair_dict = {}
i = 0
for s in gene_pair:
    if not re.search('\w+/\w+',s):
        continue
    gene_pair_dict[s] = i
    i = i+1

for line in lines:
    cancer_key = line[0]+'/'+line[1]
    gene_pair_key = line[2].rstrip()
    if (cancer_key in cancer_dict) and (gene_pair_key in gene_pair_dict):
##        if row == 1601:
##            print line
##            print str(row)+'\t'+str(col)
        data[gene_pair_dict[gene_pair_key]][cancer_dict[cancer_key]] = 1
##        if data[1601][28]==1:
##            print line
##            print str(row)+'\t'+str(col)
##            break

outFile = open('/Users/HSH/Roche/20140619/cancerFusionMatrix','w')
for i in range(len(cancer_row)):
    top,morph = cancer_row[i].split('/')
    query = 'SELECT Benamning FROM Koder WHERE Kod=%s' % top
    cursor.execute(query)
    query_result = cursor.fetchall()
    cancer_row[i] = str(query_result[0][0])
    query = 'SELECT Benamning FROM Koder WHERE Kod=%s' % morph
    cursor.execute(query)
    query_result = cursor.fetchall()
    cancer_row[i] += ' '
    cancer_row[i] += str(query_result[0][0])
    cancer_row[i] = cancer_row[i].replace(' ','_')
line_out = ''+'\t'+'\t'.join(cancer_row)+'\n'
outFile.write(line_out)
for s in gene_pair_dict.keys():
    line_out = s+'\t'+'\t'.join(map(str,data[gene_pair_dict[s]]))+'\n'
    outFile.write(line_out)

outFile.close()
cnx.close()
