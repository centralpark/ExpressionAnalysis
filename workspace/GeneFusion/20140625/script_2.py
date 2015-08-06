'''
Created on Jun 25, 2014
@author: HSH
'''

'''Create data for CCDC6-RET fusion in THCA'''
gene5p = 'CCDC6'
gene3p = 'RET'
cancer = 'THCA'

f = open('/Users/HSH/Roche/Data/all8.txt')
f.readline()
file_out = open('/Users/HSH/Roche/20140625/CCDC6-RET_THCA','w')
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

f.close()
file_out.close()