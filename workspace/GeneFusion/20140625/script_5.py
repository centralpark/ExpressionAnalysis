'''
Created on Jun 25, 2014

@author: HSH
'''

import string
import re
import functools
import mysql.connector


filename = '/Users/HSH/Roche/Data/all14.txt'
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
lines = filter(lambda x:x[col_tumor]=='LUAD',lines)

# should not occur in normal tissue
lines = filter(lambda x:not x[col_numNormal],lines)

# filter paralog
lines = filter(lambda x:x[col_isParalog]=='0',lines)

# filter number of study
lines = filter(lambda x:x[col_numStudy],lines)

# prepare data file for circos plot
outFile = open('/Users/HSH/Roche/20140625/circos_data_LUAD','w')
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

    
cnx.close() 
outFile.close()