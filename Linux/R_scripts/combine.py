"""
Combine the statistical result as new columns to all15.txt

Command line parameter:
argv[1]: file of age statistical result
argv[2]: file of OS statistical result
argv[3]: all15.txt file
argv[4]: output file

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
    
    datMain['Age_Shift'] = pd.Series([np.nan]*nrow, index=datMain.index)
    datMain['Age_pValue'] = pd.Series([np.nan]*nrow, index=datMain.index)
    datMain['OS_Shift'] = pd.Series([np.nan]*nrow, index=datMain.index)
    datMain['OS_pValue'] = pd.Series([np.nan]*nrow, index=datMain.index)    
    for i in range(nrow):
        sel_row = datAge.ix[(datAge.Cancer==datMain.ix[i,'Tumor']) & (datAge.Gene5p==datMain.ix[i,'gene5p']) & (datAge.Gene3p==datMain.ix[i,'gene3p']),['pValue','Shift']]
        if len(sel_row.index) < 1:
            continue
        ind = sel_row.index[0]
        datMain.ix[i,'Age_Shift'] = datAge.ix[ind,'Shift']
        datMain.ix[i,'Age_pValue'] = datAge.ix[ind,'pValue']
        
        sel_row = datOS.ix[(datOS.Cancer==datMain.ix[i,'Tumor']) & (datOS.Gene5p==datMain.ix[i,'gene5p']) & (datOS.Gene3p==datMain.ix[i,'gene3p']),['pValue','Shift']]
        if len(sel_row.index) < 1:
            continue
        ind = sel_row.index[0]
        datMain.ix[i,'OS_Shift'] = datOS.ix[ind,'Shift']
        datMain.ix[i,'OS_pValue'] = datOS.ix[ind,'pValue']
    datMain.to_csv(file_out,sep='\t',index=False)