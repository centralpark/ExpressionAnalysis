'''
Prepare data for R DESeq2 package

Created on Jul 2, 2014
@author: HSH
'''

import os
import pandas as pd

message_flag = True

def listand(listA,listB):
    if len(listA) != len(listB):
        raise('The length of two lists must be the same!')
    N = len(listA)
    result = [None]*N
    for i in range(N):
        result[i] = listA[i] and listB[i]
    return result
        

def prepareColdata(gct_path,ofile_path,fusion=None,cancer=None,map_path='/Users/HSH/Roche/Data/connect/ccle_tcga_mapping.txt',
                   all8_path='/Users/HSH/Roche/Data/all8.txt'):
    '''
    Prepare colData for constructing DESeqDataSet from scratch
    Input
        -- fusion: can be gene fusion pair, e.g. ('EML4','ALK'), or several fusions grouped by 3' gene, e.g. (['EML4','AK7'],'ALK')
    '''    
    if not cancer:
        filename = os.path.split(gct_path)[1]
        cancer = filename.split('.')[0]
        
    if not fusion:
        raise Exception('Provide a fusion to look at!')
    
    datAll = pd.io.parsers.read_table(all8_path,low_memory=False)
    datCancer = datAll.ix[datAll.Tumor==cancer,['Tumor','gene5p','gene3p','PatientID','Sample','age','gender']]
    
    # get all barcodes in expression data file
    f = open(gct_path)
    samples = f.readline().rstrip().split('\t')[1:]
    f.close()
    
    col_data = {}
    g5p = fusion[0]
    g3p = fusion[1]
    for s in samples:
        pid = s[0:12]
        sample = int(s[13:15])
        if sample < 10:
            sampleType = 'T'
        elif sample < 20:
            sampleType = 'N'
        else:
            sampleType = 'C'
        select = listand((datCancer.Sample==sampleType).values,(datCancer.PatientID==pid).values)
        if any(select):
            dat_s = datCancer.ix[select,]
            # change the column location accordingly            
            age = dat_s.iloc[0,5]
            if pd.isnull(age):
                age = ''
            else:
                age = str(age)
            gender = dat_s.iloc[0,6]
            if pd.isnull(gender):
                gender = ''
            select_2 = listand((dat_s.gene5p==g5p).values,(dat_s.gene3p==g3p).values)
            if any(select_2):
                harbor_fusion = 'Y'
            else:
                harbor_fusion = 'N'
        else:
            print('No matched aliquot_id:'+s)
            continue
        col_data[s] = [pid,sampleType,harbor_fusion,age,gender]
    
    ofile = open(ofile_path,'w')
    header = '\t'.join(['','patient','tumor','harbor_fusion','age','gender'])
    ofile.write(header+'\n')
    for s in samples:
        if s not in col_data:
            continue
        line_out = '\t'.join([s]+col_data[s])
        ofile.write(line_out+'\n')    
    ofile.close()
    
    return col_data

if __name__ == '__main__':    
    pass
#    cancer = 'THCA'
#    gene = 'MCM9'
#    partner = 'ASF1A'
#    date = '20140829'
#    gct_file = '/Users/HSH/Roche/Data/expression/' + cancer + '.gct'
#    col_file = '/Users/HSH/Roche/' + date + '/' + '_'.join([cancer,gene,'colData'])
#    
#    result = prepareColdata(gct_file,col_file,fusion=(partner,gene),cancer=cancer,
#                            map_path='/Users/HSH/Roche/Data/connect/ccle_tcga_mapping_v1.txt')