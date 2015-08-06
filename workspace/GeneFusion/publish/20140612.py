'''
Created on Jun 12, 2014

@author: HSH
'''

import string
import re
import functools
import mysql.connector
from collections import Counter


filename = '/Users/HSH/Roche/20140610/all14.txt'
f = open(filename,'r')
# header line
f.readline()
col_tumor = 0
col_gene5p = 2
col_gene3p = 3
col_isParalog = 8
col_numNormal = 10
col_avgRead = 23
col_polyA = 24
col_cnv3p = 28
col_cnv5p = 29
col_tumorAge = 31
col_numStudy = 34
col_numStudyFlex = 35
lines = f.readlines()
f.close()

lines = map(functools.partial(string.split,sep='\t'),lines)

# focus on a specific cancer type
lines = filter(lambda x:x[col_tumor]=='GBM',lines)

# should not occur in normal tissue
lines = filter(lambda x:not x[col_numNormal],lines)

# filter paralog
lines = filter(lambda x:x[col_isParalog]=='0',lines)

# filter number of study
lines = filter(lambda x:x[col_numStudyFlex],lines)

# filter copy number variation (CNV)
lines = filter(lambda x:-0.1<float(x[col_cnv3p])<0.1,lines)

# filter polyA
lines = filter(lambda x:int(float(x[col_polyA]))<=2,lines)

# filter average number of support reads
lines = filter(lambda x:float(x[col_avgRead])>=1,lines)

# prepare data file for circos plot
outFile = open('/Users/HSH/Roche/20140611/circos_data_GBM','w')
cnx = mysql.connector.connect(user='root',database='ROCHE')
cursor = cnx.cursor()
partners = []

for line in lines:
    gene5p = line[col_gene5p]
    gene3p = line[col_gene3p]
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

hist_data_file = open('/Users/HSH/Roche/20140611/circos_hist_data','w')
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
outFile.close()
hist_data_file.close()