"""
Created on Thu Sep 25 11:25:45 2014
@author: HSH
"""

import sys
import pandas as pd
import numpy as np

if __name__ == '__main__':
    # Parsing command line arguments
    fileAge = sys.argv[1]
    fileOS = sys.argv[2]
    file_in = sys.argv[3]    
    file_out = sys.argv[4]    
    
    datMain = pd.io.parsers.read_table(file_in,
                                       low_memory=False)
    datAge = pd.io.parsers.read_table(fileAge)
    datOS = pd.io.parsers.read_table(fileOS)
    nrow = datMain.shape[0]
    
    datMain['Age_Shift_Grouped'] = pd.Series([np.nan]*nrow, index=datMain.index)
    datMain['Age_pValue_Grouped'] = pd.Series([np.nan]*nrow, index=datMain.index)
    datMain['OS_Shift_Grouped'] = pd.Series([np.nan]*nrow, index=datMain.index)
    datMain['OS_pValue_Grouped'] = pd.Series([np.nan]*nrow, index=datMain.index)
    for i in range(nrow):
        sel_row = datAge.ix[(datAge.Cancer==datMain.ix[i,'Tumor']) & (datAge.Gene3p==datMain.ix[i,'gene3p']),['pValue','Shift']]
        if len(sel_row.index) < 1:
            continue
        ind = sel_row.index[0]
        datMain.ix[i,'Age_Shift_Grouped'] = datAge.ix[ind,'Shift']
        datMain.ix[i,'Age_pValue_Grouped'] = datAge.ix[ind,'pValue']
        
        sel_row = datOS.ix[(datOS.Cancer==datMain.ix[i,'Tumor']) & (datOS.Gene3p==datMain.ix[i,'gene3p']),['pValue','Shift']]
        if len(sel_row.index) < 1:
            continue
        ind = sel_row.index[0]
        datMain.ix[i,'OS_Shift_Grouped'] = datOS.ix[ind,'Shift']
        datMain.ix[i,'OS_pValue_Grouped'] = datOS.ix[ind,'pValue']
    datMain.to_csv(file_out,sep='\t',index=False)