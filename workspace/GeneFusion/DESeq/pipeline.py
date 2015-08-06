'''
Pipeline work (Python part) for differential expression analysis

Created on Jul 28, 2014
@author: HSH
'''

import sys
sys.path.append('/Users/HSH/Roche/workspace/GeneFusion/DESeq')

from getExpressionData import writeTCGAExpressionMatrix
from prepareDESeqData import prepareColdata

cancer = 'BLCA'
gene = 'TACC3'
partner = 'FGFR3'
date = '20140908'
gct_file = '/Users/HSH/Roche/Data/expression/' + cancer + '.gct'
col_file = '/Users/HSH/Roche/' + date + '/' + '_'.join([cancer,gene,'colData'])

writeTCGAExpressionMatrix([cancer], folder='/Users/HSH/Roche/Data/expression/')
result = prepareColdata(gct_file,col_file,fusion=(partner,gene),cancer=cancer,
                        map_path='/Users/HSH/Roche/Data/connect/ccle_tcga_mapping_v1.txt')