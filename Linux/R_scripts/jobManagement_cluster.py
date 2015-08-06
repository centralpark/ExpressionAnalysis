'''
Statistical analysis for all gene fusion pairs

Command line parameter:
argv[1]: temporary folder to store job data
argv[2]: number of compute nodes to use
argv[3]: temporary folder to store age file
argv[4]: temporary folder to store OS file

Created on Aug 8, 2014
@author: HSH
'''

import re,sys,os,random,subprocess
from subprocess import Popen,PIPE,STDOUT
sys.path.append('/pred_cello1/analysis/geneFusion/tcga/integrate/R_scripts/python_library/GeneFusion')
from misc import randompartition
from all8_header import *

def sanityCheck(x):
    '''
    Check if the gene names or tumor names follow normal standard.
    '''
    if type(x) is str:
        matchObj = re.search('[^A-Za-z0-9]',x)
        if matchObj:
            return False
        else:
            return True
    elif type(x) is list:
        for y in x:
            if not sanityCheck(y):
                return False
        return True
    else:
        return False


def getallfusion(filename='/pred_cello1/analysis/geneFusion/tcga/integrate/all15.txt'):
    '''
    Find all validated gene fusions
    Output:
        -- fusion:
            collection of fusion
            Data Type: Set
    '''
    fusion = set([])
    f = open(filename)
    f.readline()
    for line in f.readlines():
        line_split = line.split('\t')
        fusion_this = [line_split[i] for i in [0,2,3]]
        # check for white space (null) values
        check_val = map(lambda x: x.strip()=='',fusion_this)
        if any(check_val):
            continue
        if not sanityCheck(fusion_this):
            continue
        fusion.add(tuple(fusion_this))
    f.close()
    return fusion
    
    
if __name__ == '__main__':
    # Parsing command line arguments
    dirJob = sys.argv[1]
    N = int(sys.argv[2])
    temp_age_folder = sys.argv[3]
    temp_OS_folder = sys.argv[4]
    
    from collections import defaultdict
    
    # Running on rnuuspr66 server
    # temporary folder to store qsub job related files
    if not os.path.exists(dirJob):
        os.mkdir(dirJob)
    fusion = getallfusion()
    fusion = list(fusion)
    
    # randomly partition the fusion list into N nearly equal parts
    fusionPartition = randompartition(fusion,N)
    for i in range(N):
        name = 'Job' + str(i) + '_datafile'
        filename = os.path.join(dirJob,name)
        ofile = open(filename,'w')
        fp = fusionPartition[i]
        for f in fp:
            line = '\t'.join([f[0],f[1],f[2]])
            ofile.write(line+'\n')
        ofile.close()
    
    for i in range(N):
        name = 'Job' + str(i) + '_datafile'
        filename = os.path.join(dirJob,name)
        qsubString = '''#!/bin/sh
#PBS -N AgeOSJob%d
#PBS -m e
#PBS -V
#PBS -q long
#PBS -j oe
#PBS -l select=1:ncpus=12:mem=12gb

python /pred_cello1/analysis/geneFusion/tcga/integrate/R_scripts/compareDistAllFusion_cluster.py %s %s %s
        ''' % (i,filename,temp_age_folder,temp_OS_folder)
        p = Popen(['qsub'],stdin=PIPE,stdout=PIPE)
        qsub_stdout = p.communicate(input=qsubString)[0]
        print(qsub_stdout)
