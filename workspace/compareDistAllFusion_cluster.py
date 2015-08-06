'''
Created on Aug 8, 2014
@author: HSH
'''

import re,sys,os
import multiprocessing as mp

def compareage(arg):
    '''
    Compare the age distribution of population with the fusion gene5p-gene3p or X-gene3p in cancer
    and population without the fusion.
    If gene5p has more than one gene symbol, all X-gene3p fusion will be grouped and analyzed together.
    ''' 
    gene5p = arg[1]
    gene3p = arg[2]
    cancer = arg[0]
    folder = arg[3]
    group = arg[4]
    
    col_gene5p = 1
    col_gene3p = 2
    col_Tumor = 3
    col_PatientID = 15
    col_Sample = 16
    col_age = 24
    col_OS = 27
    
    f = open('/pred_cello1/analysis/geneFusion/tcga/integrate/all8.txt')
    f.readline()
    if group:
        ofilename = os.path.join(folder, '_'.join([cancer,gene3p]))
    else:
        ofilename = os.path.join(folder, '_'.join([cancer,gene5p,gene3p]))
    file_out = open(ofilename,'w')
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
            elif isinstance(gene5p,list) or isinstance(gene5p,tuple):
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
    
    file_out.close()
    
    if (n_Y == 0) or (n_N == 0):
        os.remove(ofilename)    
        
    f.close()


def compareallage(fusions,folder,group=False):
    '''
    Input:
        -- fusions: a list of fusions
        -- folder: temperory folder to store the files
    '''
    pool = mp.Pool(processes=12)
    if not os.path.exists(folder):
        os.mkdir(folder)
    arg_list = []
    for fusion in fusions:
        arg_list.append(list(fusion) + [folder] + [group])
    pool.map(compareage,arg_list)


def compareOS(arg):
    '''
    Write file for the fusion in arg
    '''
    gene5p = arg[1]
    gene3p = arg[2]
    cancer = arg[0]
    folder = arg[3]
    group = arg[4]
    
    col_gene5p = 1
    col_gene3p = 2
    col_Tumor = 3
    col_PatientID = 15
    col_Sample = 16
    col_age = 24
    col_OS = 27
    
    f = open('/pred_cello1/analysis/geneFusion/tcga/integrate/all8.txt')
    f.readline()
    if group:
        ofilename = os.path.join(folder, '_'.join([cancer,gene3p]))
    else:
        ofilename = os.path.join(folder, '_'.join([cancer,gene5p,gene3p]))
    file_out = open(ofilename,'w')
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
            elif isinstance(gene5p,list) or isinstance(gene5p,tuple):
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
    
    file_out.close()
    
    if (n_Y == 0) or (n_N == 0):
        os.remove(ofilename)    
        
    f.close()


def compareallOS(fusions,folder,group=False):
    '''
    Input:
        -- fusions: a list of fusions
        -- folder: temperory folder to store the files
    '''
    pool = mp.Pool(processes=12)
    if not os.path.exists(folder):
        os.mkdir(folder)
    arg_list = []
    for fusion in fusions:
        arg_list.append(list(fusion) + [folder] + [group])
    pool.map(compareOS,arg_list)


if __name__ == '__main__':
    # Parsing command line arguments
    file_input = sys.argv[1]
    temp_folder_age = sys.argv[2]
    temp_folder_OS = sys.argv[3]
    
    f = open(file_input)
    fusions = []
    for line in f.readlines():
        line_split = line.rstrip().split('\t')
        fusion_this = [line_split[i] for i in [0,1,2]]
        # split multiple gene5p
        gene5p = fusion_this[1]
        gene5p_split = gene5p.split(',')
        if len(gene5p_split) == 1:
            fusion_this[1] = gene5p
        else:
            fusion_this[1] = gene5p_split
        fusions.append(fusion_this)
    f.close()
    compareallage(fusions,temp_folder_age,group=False)
    compareallOS(fusions,temp_folder_OS,group=True)