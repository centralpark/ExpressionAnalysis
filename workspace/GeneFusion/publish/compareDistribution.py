'''
Created on Jul 24, 2014

@author: HSH
'''

import rpy2.robjects as robjects
import sys
sys.path.append('/Users/HSH/Roche/workspace/GeneFusion')
from all8_header import *



def getallfusion(filename='/Users/HSH/Roche/Data/all15.txt'):
    '''
    Find all validated gene fusions
    '''
    fusion = set([])
    if filename.endswith(('xls','xlsx')):
        import xlrd
        workbook = xlrd.open_workbook(filename)
        worksheet = workbook.sheet_by_index(0)
        num_rows = worksheet.nrows - 1
        curr_row = 0
        while curr_row < num_rows:
            curr_row += 1
            curr_fusion = [worksheet.cell_value(curr_row, i) for i in [0,4,5]]
            # convert unicode to str to avoid problems
            if isinstance(curr_fusion[0],unicode):
                curr_fusion = map(str,curr_fusion)
            if tuple(curr_fusion) in fusion:
                print curr_fusion
            fusion.add(tuple(curr_fusion))
    else:
        f = open(filename)
        f.readline()
        for line in f.readlines():
            line_split = line.split('\t')
            numStudy = line_split[34]
            if numStudy:
                fusion_this = [line_split[i] for i in [0,2,3]]
                fusion.add(tuple(fusion_this))
        f.close()
    return fusion


def compareage(gene5p,gene3p,cancer):
    '''Compare the age distribution of population with the fusion gene5p-gene3p or X-gene3p in cancer
and population without the fusion.
If gene5p has more than one gene symbol, all X-gene3p fusion will be grouped and analyzed together.'''
    R = robjects.r
    f = open('/Users/HSH/Roche/Data/all8.txt')
    f.readline()
    file_out = open('/Users/HSH/Roche/temp','w')
    # write header
    line_out = '\t'.join(['patient_id','age','harbor_fusion']) + '\n'
    file_out.write(line_out)
    pid = set([])
    pid_age = {}
    for line in f.readlines():
        line_split = line.split('\t')
        pid_this = line_split[col_PatientID]
        skip = False
        try:
            age_num = int(line_split[col_age])
        except ValueError:
            skip = True
        if skip:
            continue
        if line_split[col_Tumor]==cancer:
            pid_age[pid_this] = line_split[col_age]
            # if only one gene5p is provided, paired fusion is analyzed, otherwise grouped
            if type(gene5p) is str:
                if [line_split[i] for i in [col_gene5p,col_gene3p,col_Sample]] == [gene5p,gene3p,'T']:
                    pid.add(pid_this)
            elif type(gene5p) is list or tuple:
                for g in gene5p:
                    if [line_split[i] for i in [col_gene5p,col_gene3p,col_Sample]] == [g,gene3p,'T']:
                        pid.add(pid_this)
                        continue
            else:
                print gene5p
                raise Exception('Gene5p type unknown!')
    
    n_Y = 0     # number of sample harboring the fusion
    n_N = 0     # number of sample not harboring the fusion
    
    for pid_this in pid_age.keys():
        if pid_this in pid:
            line_out = pid_this+'\t'+pid_age[pid_this]+'\t'+'Y'+'\n'
            n_Y += 1
        else:
            line_out = pid_this+'\t'+pid_age[pid_this]+'\t'+'N'+'\n'
            n_N += 1
        file_out.write(line_out)
    
    if (n_Y == 0) or (n_N == 0):
        return None    
        
    f.close()
    file_out.close()

    R_script = '''
    dat <- read.table('/Users/HSH/Roche/temp',sep='\t',header=T)
    population_1 <- dat[dat[,3]=='Y',2]
    population_2 <- dat[dat[,3]=='N',2]
    result <- wilcox.test(age~harbor_fusion,data=dat)
    '''
    R(R_script)
    R_result = robjects.globalenv['result']
    W = R_result.rx('statistic')[0][0]
    p_value = R_result.rx('p.value')[0][0]
    return [n_Y,n_N,W,p_value]


def compareallage(ifilename,ofilename,group=True):
    all_fusion = getallfusion(ifilename)
    
    # Analyze age influence for all documented gene fusions
    file_result = open(ofilename,'w')
    if group:
        file_result.write('\t'.join(['Cancer','Gene5p','Gene3p','SampleSizeWithFusion',\
                                 'SampleSizeWithoutFusion','W','pValue']) + '\n')
        from collections import defaultdict
        grouped_fusion = defaultdict(list)
        for cancer,gene5p,gene3p in all_fusion:
            grouped_fusion[(cancer,gene3p)].append(gene5p)
        for key in grouped_fusion:
            gene5p = grouped_fusion[key]
            gene3p = key[1]
            cancer = key[0]
            result = compareage(gene5p,gene3p,cancer)
            if not result:
                continue
            line_out = '\t'.join([cancer,','.join(gene5p),gene3p]+map(str,result))
            file_result.write(line_out+'\n')
    else:
        file_result.write('\t'.join(['Cancer','Gene5p','Gene3p','SampleSizeWithFusion',\
                                 'SampleSizeWithoutFusion','W','pValue']) + '\n')
        for fusion in all_fusion:
            gene5p = fusion[1]
            gene3p = fusion[2]
            cancer = fusion[0]
            result = compareage(gene5p,gene3p,cancer)
            if not result:
                continue
            line_out = '\t'.join([cancer,gene5p,gene3p]+map(str,result))
            file_result.write(line_out+'\n')
    file_result.close()


def compareOS(gene5p,gene3p,cancer):
    '''
    Compare overall survival distribution.
    '''
    R = robjects.r
    f = open('/Users/HSH/Roche/Data/all8.txt')
    f.readline()
    file_out = open('/Users/HSH/Roche/temp','w')
    # write header
    line_out = '\t'.join(['patient_id','OS','harbor_fusion']) + '\n'
    file_out.write(line_out)

    pid = set([])
    pid_OS = {}
    for line in f.readlines():
        line_split = line.split('\t')
        pid_this = line_split[col_PatientID]
        skip = False
        try:
            int(line_split[col_OS])
        except ValueError:
            skip = True
        if skip:
            continue
        if line_split[col_Tumor]==cancer:
            pid_OS[pid_this] = line_split[col_OS]
            # if only one gene5p is provided, paired fusion is analyzed, otherwise grouped
            if type(gene5p) is str:
                if [line_split[i] for i in [col_gene5p,col_gene3p,col_Sample]] == [gene5p,gene3p,'T']:
                    pid.add(pid_this)
            elif type(gene5p) is list:
                for g in gene5p:
                    if [line_split[i] for i in [col_gene5p,col_gene3p,col_Sample]] == [g,gene3p,'T']:
                        pid.add(pid_this)
                        continue
            else:
                raise Exception('Gene5p type unknown!')
    
    n_Y = 0     # number of sample harboring the fusion
    n_N = 0     # number of sample not harboring the fusion
    
    for pid_this in pid_OS.keys():
        if pid_this in pid:
            line_out = pid_this+'\t'+pid_OS[pid_this]+'\t'+'Y'+'\n'
            n_Y += 1
        else:
            line_out = pid_this+'\t'+pid_OS[pid_this]+'\t'+'N'+'\n'
            n_N += 1
        file_out.write(line_out)
    
    if (n_Y == 0) or (n_N == 0):
        return None    
        
    f.close()
    file_out.close()

    R_script = '''
    dat <- read.table('/Users/HSH/Roche/temp',sep='\t',header=T)
    population_1 <- dat[dat[,3]=='Y',2]
    population_2 <- dat[dat[,3]=='N',2]
    result <- wilcox.test(OS~harbor_fusion,data=dat,conf.int=T)
    '''
    R(R_script)
    R_result = robjects.globalenv['result']
    W = R_result.rx('statistic')[0][0]
    p_value = R_result.rx('p.value')[0][0]
    estimate = R_result.rx('estimate')[0][0]
    return [n_Y,n_N,W,p_value,estimate]


def compareallOS(ifilename,ofilename,group=True):
    all_fusion = getallfusion(ifilename)
    
    # Analyze OS influence for all documented gene fusions
    file_result = open(ofilename,'w')
    if group:
        file_result.write('\t'.join(['Cancer','Gene5p','Gene3p','SampleSizeWithFusion',\
                                 'SampleSizeWithoutFusion','W','pValue','shift']) + '\n')
        from collections import defaultdict
        grouped_fusion = defaultdict(list)
        for cancer,gene5p,gene3p in all_fusion:
            grouped_fusion[(cancer,gene3p)].append(gene5p)
        for key in grouped_fusion:
            gene5p = grouped_fusion[key]
            gene3p = key[1]
            cancer = key[0]
            result = compareOS(gene5p,gene3p,cancer)
            if not result:
                continue
            line_out = '\t'.join([cancer,','.join(gene5p),gene3p]+map(str,result))
            file_result.write(line_out+'\n')
    else:
        file_result.write('\t'.join(['Cancer','Gene5p','Gene3p','SampleSizeWithFusion',\
                                 'SampleSizeWithoutFusion','W','pValue','shift']) + '\n')
        for fusion in all_fusion:
            gene5p = fusion[1]
            gene3p = fusion[2]
            cancer = fusion[0]
            result = compareOS(gene5p,gene3p,cancer)
            if not result:
                continue
            line_out = '\t'.join([cancer,gene5p,gene3p]+map(str,result))
            file_result.write(line_out+'\n')
    file_result.close()

if __name__ == '__main__':
#     all_fusion = getallfusion('/Users/HSH/Roche/Data/filteredFusions_v2.xls')
#     compareallage('/Users/HSH/Roche/Data/filteredFusions_v2.xls','/Users/HSH/Roche/20140805/age_statisticalAnalysis_v1.txt',group=False)
#     compareallOS('/Users/HSH/Roche/Data/filteredFusions_v2.xls','/Users/HSH/Roche/20140805/OS_statisticalAnalysis_v1.txt',group=False)
    pass