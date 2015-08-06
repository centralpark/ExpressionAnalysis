'''
Prepare data for R DESeq2 package
Created on Jul 2, 2014
@author: HSH
'''

import os,re

message_flag = True

def prepareColdata(gct_path,ofile_path,fusion=None,cancer=None,map_path='/Users/HSH/Roche/Data/connect/ccle_tcga_mapping.txt',
                   all8_path='/Users/HSH/Roche/Data/all8.txt'):
    '''
    Prepare colData for constructing DESeqDataSet from scratch
    Input
        -- fusion: can be gene fusion pair, e.g. ('EML4','ALK'), or several fusions grouped by 3' gene, e.g. ('
    '''
    global message_flag
    barcode_exception = set([])
    
    if not cancer:
        filename = os.path.split(gct_path)[1]
        cancer = filename.split('.')[0]
        
    if not fusion:
        raise Exception('Provide a fusion to look at!')
    
    map_dict = {}
    f = open(map_path)
    f.readline()
    for line in f.readlines():
        line_split = line.split('\t')
        analysis_id = line_split[16].lower()
        barcode = line_split[1]
        map_dict[analysis_id] = barcode
    f.close()
    
    f = open(gct_path)
    samples = f.readline().rstrip().split('\t')[1:]
    f.close()
    
    col_data = {}
    pid2barcode = {}
    for s in samples:
        pid = s[0:12]
        col_data[s] = [pid,'T','N','','']
        pid2barcode[pid] = s
    
    f = open(all8_path)
    f.readline()
    # Test whether grouping is needed
    if type(fusion[0]) is str:
        for line in f.readlines():
            line_split = line.split('\t')
            tumor = line_split[3]
            analysis_id = line_split[0]
            if tumor != cancer or analysis_id not in map_dict:
                continue
            gene_5p = line_split[1]
            gene_3p = line_split[2]
            barcode = map_dict[analysis_id]
            # Some barcodes are not TCGA, e.g. 'CCLE-S-117-RNA-08'
            if not re.match('TCGA',barcode):
                continue
            # Some samples do not have gene expression data available, see example in COAD
            if barcode not in col_data:
                if message_flag:
                    print 'The following barcodes mismatches between gct and all8.txt:'
                    message_flag = False
                if barcode not in barcode_exception:
                    print barcode
                    barcode_exception.add(barcode)
                pid = barcode[0:12]
                if pid in pid2barcode:
                    barcode = pid2barcode[pid]
                else:
                    continue
                
            col_data[barcode][1] = line_split[16]
            col_data[barcode][3] = line_split[24]
            col_data[barcode][4] = line_split[26]
            if gene_5p == fusion[0] and gene_3p == fusion[1]:
                col_data[barcode][2] = 'Y'
    else:
        for line in f.readlines():
            line_split = line.split('\t')
            tumor = line_split[3]
            analysis_id = line_split[0]
            if tumor != cancer or analysis_id not in map_dict:
                continue
            gene_5p = line_split[1]
            gene_3p = line_split[2]
            barcode = map_dict[analysis_id]
            # Some barcodes are not TCGA, e.g. 'CCLE-S-117-RNA-08'
            if not re.match('TCGA',barcode):
                continue
            col_data[barcode][1] = line_split[16]
            col_data[barcode][3] = line_split[24]
            col_data[barcode][4] = line_split[26]
            if gene_5p in fusion[0] and gene_3p == fusion[1]:
                col_data[barcode][2] = 'Y'
    f.close()
    
    ofile = open(ofile_path,'w')
    header = '\t'.join(['','patient','tumor','harbor_fusion','age','gender'])
    ofile.write(header+'\n')
    for s in samples:
        line_out = '\t'.join([s]+col_data[s])
        ofile.write(line_out+'\n')    
    ofile.close()
    
    return [col_data,pid2barcode]

if __name__ == '__main__':
    pass