"""
Created on Tue Sep  9 10:41:49 2014

@author: HSH
"""

import os
import sys
sys.path.append('/home/staff/hes14/Roche/workspace/GeneFusion/DESeq')

from getExpressionData import writeTCGAExpressionMatrix
from prepareDESeqData import prepareColdata

def preparedata(cancer,gene5p,gene3p,folder=None):
    if not folder:
        folder = os.getcwd()
    gct_file = os.path.join(folder,cancer + '.gct')
    if not os.path.isfile(gct_file):
        writeTCGAExpressionMatrix([cancer], folder=folder)
    col_file = os.path.join(folder,'_'.join([cancer,gene5p,gene3p,'colData']))
    prepareColdata(gct_file,col_file,fusion=(gene5p,gene3p),cancer=cancer,
                   map_path='/home/staff/hes14/Roche/Data/connect/ccle_tcga_mapping_v1.txt')

if __name__ == '__main__':
    cancer = sys.argv[1]
    gene5p = sys.argv[2]
    gene3p = sys.argv[3]
    folder = sys.argv[4]
    preparedata(cancer,gene5p,gene3p,folder)