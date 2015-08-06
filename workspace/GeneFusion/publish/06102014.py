'''
Created on Jun 10, 2014

@author: HSH
'''

import string,re
import mysql.connector
from collections import Counter

filename = '/Users/HSH/Roche/20140610/all14.txt'
f = open(filename,'r')
# header line
f.readline()
outFile = open('/Users/HSH/Roche/20140610/all_fusion_location','w')

cnx = mysql.connector.connect(user='root',database='ROCHE')
cursor = cnx.cursor()

partners = []

for line in f:
    columns = line.split('\t')
    cancer = columns[0]
    gene5p = columns[2]
    gene3p = columns[3]
    isParalog = columns[8]
    # consider only GBM
    if cancer != 'GBM':
        continue
    if isParalog == 'B':
        continue
    query = 'SELECT chrom,txStart,txEnd FROM refGene WHERE name2=%s'
    cursor.execute(query,(gene5p,))
    query_result = cursor.fetchall()
    if len(query_result) < 0.5:
        continue
    (chr5p,start5p,end5p) = query_result[0]
    chr5p = string.lower(chr5p)
    if start5p > end5p:
        start5p,end5p = end5p,start5p
    query = 'SELECT chrom,txStart,txEnd FROM refGene WHERE name2=%s'
    cursor.execute(query,(gene3p,))
    query_result = cursor.fetchall()
    if len(query_result) < 0.5:
        continue
    (chr3p,start3p,end3p) = query_result[0]
    chr3p = string.lower(chr3p)
    if start3p > end3p:
        start3p,end3p = end3p,start3p
    pattern = re.compile('chr[\d+]_')
    if pattern.match(chr5p) or pattern.match(chr3p):
        continue
    # calculate histogram data file for Circos
    partners.extend([gene5p,gene3p])
    line_out = []
    line_out = string.join([chr5p,str(start5p),str(end5p),chr3p,str(start3p),str(end3p)], '\t')
    line_out = line_out.replace('chr','hs')
    outFile.write(line_out +'\n')

hist_data_file = open('/Users/HSH/Roche/20140610/hist_data','w')
cnt = Counter(partners)
for gene in cnt:
    query = 'SELECT chrom,txStart,txEnd FROM refGene WHERE name2=%s'
    cursor.execute(query,(gene,))
    query_result = cursor.fetchall()
    if len(query_result) < 0.5:
        continue
    (chromo,start,end) = query_result[0]
    chromo = chromo.replace('chr','hs')
    if start > end:
        start,end = end,start
    line_out = string.join([chromo,str(start),str(end),str(cnt[gene])],'\t')
    hist_data_file.write(line_out + '\n')


cnx.close() 
f.close()
outFile.close()
hist_data_file.close()