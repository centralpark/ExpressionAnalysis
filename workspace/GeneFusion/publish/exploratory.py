'''
Created on Jun 4, 2014

@author: Siheng He
'''

import string
import mysql.connector

filename = '/Users/HSH/Roche/all16.txt'
f = open(filename,'r')
# header line
f.readline()
outFile = open('/Users/HSH/Roche/all_fusion_location','a')

cnx = mysql.connector.connect(user='root',database='ROCHE')
cursor = cnx.cursor()

for line in f:
    columns = line.split('\t')
    gene5p = columns[2]
    gene3p = columns[3]
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
    line_out = []
    line_out = string.join([chr5p,str(start5p),str(end5p),chr3p,str(start3p),str(end3p)], '\t')
    line_out = line_out.replace('chr','hs')
    outFile.write(line_out +'\n')
    
cnx.close() 
f.close()
outFile.close()