'''
Created on Jun 25, 2014
@author: HSH
'''

import rpy2.robjects as robjects

# map the column header name to numeric index
# file all8.txt
col_analysis_id = 0
col_gene5p = 1
col_gene3p = 2
col_Tumor = 3
col_Tool = 4
col_gene5p_start = 5
col_gene5p_end = 6
col_gene3p_start = 7
col_gene3p_end = 8
col_gene5p_break = 9
col_gene3p_break = 10
col_spanning_reads = 11
col_polyA = 12
col_coor_break3p = 13
col_coor_break5p = 14
col_PatientID = 15
col_Sample = 16
col_FalsePositive = 17
col_isKinase = 18
col_DriverMut = 19
col_CbioportalFusion = 20
col_5p_cnv = 21
col_3p_cnv = 22
col_diagnosis = 23
col_age = 24
col_averga_age_tumor_type = 25
col_gender = 26
col_exon_jump_fold = 27
col_rpkm_before_avr = 28
col_rpkm_after_avr = 29
col_exon_fusion = 30
col_perc_kinase_exon = 31
col_frame_orientation = 32
col_fusion_type = 33

def getdocumentedfusion(filename='/Users/HSH/Roche/Data/all14.txt'):
    '''Find all documented gene fusions'''
    fusion = set([])
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
    '''Compare the age distribution of population with the fusion gene5p-gene3p in cancer
and population without the fusion'''
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
        pid_this = line_split[15]
        skip = False
        try:
            age_num = int(line_split[24])
        except ValueError:
            skip = True
        if skip:
            continue
        if line_split[3]==cancer:
            pid_age[pid_this] = line_split[24]
            if [line_split[i] for i in [1,2,16]] == [gene5p,gene3p,'T']:
                pid.add(pid_this)
    
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


all_fusion = getdocumentedfusion()

# Analyze age influence for all documented gene fusions
file_result = open('/Users/HSH/Roche/20140625/age_statisticalAnalysis.txt','w')
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

