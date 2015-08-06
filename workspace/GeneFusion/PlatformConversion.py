'''
Created on Jul 10, 2014
@author: HSH
'''

import os

def fileConversion(ifile_path):
    '''
    Mac OS to Windows
    '''
    f = open(ifile_path)
    content_mac = f.read()
    content_win = content_mac.replace('~/Roche','U:/My Documents')
    name = f.name
    f.close()
    os.remove(name)
    ofile_path = name
    ofile = open(ofile_path,'w')
    ofile.write(content_win)
    ofile.close()


def mac2linux(path, cluster_name = False, ofile_path = None):
    '''
    Mac OS to Roche Linux
    Input:
        -- cluster_name: change the file name to cluster
    '''
    if os.path.isdir(path):
        file_paths = [os.path.join(path,f) for f in os.listdir(path)]
        file_paths = [f for f in file_paths if os.path.isfile(f)]
    else:
        file_paths = [path]
    for file_path in file_paths:
        f = open(file_path)
        content_mac = f.read()
        ext = os.path.splitext(f.name)[1]
        content_linux = content_mac.replace('/Users/HSH','/home/staff/hes14')
        if ext=='.py':
            if cluster_name:
                name = os.path.splitext(f.name)[0] + '_cluster.py' 
            else:
                name = os.path.splitext(f.name)[0] + '.py'
            f.close()
        elif ext=='.R':
            content_linux = content_linux.replace('~/Roche','/home/staff/hes14/Roche')
            if cluster_name:
                name = os.path.splitext(f.name)[0] + '_cluster.R' 
            else:
                name = os.path.splitext(f.name)[0] + '.R'
            f.close()
        else:
            f.close()
            continue  
        ofile_path = name
        ofile = open(ofile_path,'w')
        ofile.write(content_linux)
        ofile.close()
    return


if __name__ == '__main__':
    folderList = ['/home/staff/hes14/Roche/R_workspace','/home/staff/hes14/Roche/workspace/GeneFusion/DESeq']    
    for folder in folderList:
        print folder
        mac2linux(folder)
