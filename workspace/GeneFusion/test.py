import rpy2.robjects as robjects
import sys
sys.path.append('/Users/HSH/Roche/workspace/GeneFusion')
from all8_header import *

cancer = 'THCA'
gene5p = 'CCDC6'
gene3p = 'RET'

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
        elif type(gene5p) is list:
            for g in gene5p:
                if [line_split[i] for i in [col_gene5p,col_gene3p,col_Sample]] == [g,gene3p,'T']:
                    pid.add(pid_this)
                    continue
        else:
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

f.close()
file_out.close()