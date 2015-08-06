'''
Compile a GENE-SAMPLE expression matrix from TCGA expression files (*.rsem.genes.results)

Created on Jul 1, 2014
@author: HSH
'''

import os,re,sys
sys.path.append('/Users/HSH/Roche/workspace/GeneFusion/DESeq')

# define frequently used mappings
# probe_id -> gene symbol, gene type
global id_symbol_map
id_symbol_map = {}
f = open('/Users/HSH/Roche/Data/expression/gene_annotation_v1.txt')
f.readline()
for line in f.readlines():
    line_split = line.rstrip().split('\t')
    if len(line_split) < 2:
        continue
    probe_id,symbol,gene_type = line_split
    id_symbol_map[probe_id] = [symbol,gene_type]
f.close()

global data_type
data_type = 'raw'

def getPseudogene(ifile_path='/Users/HSH/Roche/Data/expression/gene_annotation.txt',
                  ofile_path='/Users/HSH/Roche/Data/expression/pseudogene_list',mode='data'):
    '''Get a list of pseudogenes
    Input:
        -- ifile: input file where to look for pseudogenes
    '''
    f = open(ifile_path)
    f.readline()
    if mode.lower() == 'file':
        ofile = open(ofile_path,'w')
        line_out = 'Probe_ID' + '\t' + 'Gene_Symbol'
        ofile.write(line_out)
        for line in f.readlines():
            line_split = line.rstrip().split('\t')
            if len(line_split) < 3:
                continue
            probe_id = line_split[0]
            symbol = line_split[1]
            title = line_split[2]
            if 'pseudogene' in title:
                continue
            line_out = probe_id + '\t' + symbol
            ofile.write(line_out)
        ofile.close()
        f.close()
        
        return
    
    elif mode.lower() == 'data':
        result = []
        for line in f.readlines():
            line_split = line.rstrip().split('\t')
            probe_id = line_split[0]
            if len(line_split) < 3:
                continue
            title = line_split[2]
            if 'pseudogene' in title:
                result.append(probe_id)
        f.close()
        
        return result
    
    else:
        f.close()
        raise Exception('Mode not supported!')


def readRSEM(file_path):
    '''Read individual *.rsem.gene.results file from TCGA'''
    
    f = open(file_path)
    f.readline()
    expression = {}

    for line in f.readlines():
        line_split = line.rstrip().split('\t')
        gene_id = re.search('.*\|(\d+)',line_split[0]).group(1)
        raw_count = line_split[1]
        expression[gene_id] = raw_count
        
    f.close()
    
    return expression


def readCancer(path,map_dictionary=None,data_type='raw'):
    '''Process all *.rsem.gene.results files for a single cancer type'''
    if not map_dictionary:
        raise Exception('No sample ID mapping provided!')
    
    data_files = os.listdir(path)
    cancer_expression = {}
    
    if data_type.lower() == 'raw':
        data_files = filter(lambda x: x.endswith('rsem.genes.results'),data_files)
    elif data_type.lower() == 'normalized':
        data_files = filter(lambda x: x.endswith('rsem.genes.normalized_results'),data_files)
    else:
        raise Exception('Define the correct file type!')
        
    
    for df in data_files:
        analysis_id = re.search('\w{8}-\w{4}-\w{4}-\w{4}-\w{12}',df).group()
        if analysis_id not in map_dictionary:
#            print(analysis_id)            
            continue
        sample_id = map_dictionary[analysis_id]
        fullname = os.path.join(path,df)
        sample_expression = readRSEM(fullname)
        cancer_expression[sample_id] = sample_expression
        
    return cancer_expression


def writeExpressionMatrix(path,ofile_path,map_dict,filtering=True):
    '''Write the expression matrix for a single cancer type
    Input
        -- filtering: only keep protein-coding gene
    '''
    expression = readCancer(path,map_dict)
    expression_copy = expression.copy()
    ofile = open(ofile_path,'w')
    # Make sure the expression data is for the same set of genes, which is usually the case
    gene_set = set(expression_copy.popitem()[1].keys())
    
    while expression_copy:
        sample_expression = expression_copy.popitem()[1]
        gene_set_new = set(sample_expression.keys())
        gene_set &= gene_set_new
    
    gene_list = list(gene_set)
    gene_list.sort(key=int)
    sample_list = expression.keys()
    line_out = '\t'.join(['NAME']+sample_list)
    ofile.write(line_out+'\n')
    
    if filtering:
        for g in gene_list:
            if id_symbol_map[g][1] not in ['protein-coding','6']:
                continue
            line_out = '\t'.join([id_symbol_map[g][0]]+[expression[s][g] for s in sample_list])
            ofile.write(line_out+'\n')
    else:
        for g in gene_list:
            line_out = '\t'.join([id_symbol_map[g][0]]+[expression[s][g] for s in sample_list])
            ofile.write(line_out+'\n')
    
    ofile.close()


def writeTCGAExpressionMatrix(selected_cancer=None,source_folder = '/Users/HSH/Roche/Data/tcga_data_20140528163559/',
                              map_file = '/Users/HSH/Roche/Data/connect/ccle_tcga_mapping_v1.txt',
                              folder=None,ignore_normal=True):
    '''Write expresion matrices for all cancer types in TCGA
    Input
        -- selected_cancer: write expression matrix data for selected cancer type
        -- source_folder: folder that include corresponding expression data files
        -- map_file: file that contains information mapping aliquot_id to patientID
        -- ignore_normal: Ignore the expression data for normal samples
    '''
    import glob
    
    map_dict = {}
    f = open(map_file)
    f.readline()
    
    for line in f.readlines():
        line_split = line.split('\t')
        study = line_split[0]
        if study != 'TCGA':
            continue
        sample_type = line_split[4]
        if not sample_type:
            continue
        if ignore_normal and sample_type[0] != 'T':
            continue
        # Some of the aliquot_id contains upper case letter
        aliquot_id = line_split[17].lower()
        barcode = line_split[1]
        
        map_dict[aliquot_id] = barcode
    
    f.close()
    cancer = filter(lambda x: os.path.isdir(source_folder+x), os.listdir(source_folder))
    
    if selected_cancer:
        cancer = selected_cancer
     
    for c in cancer:
        data_folder = os.path.join(source_folder,c)
        data_folder_content = glob.glob(os.path.join(data_folder,'*.rsem.genes.results'))
        if not len(data_folder_content):
            print 'Not enough *rsem.genes.results files for ' + c
            continue
        if not folder:
            folder = os.getcwd()
        ofile_path = os.path.join(folder,c+'.gct')
        writeExpressionMatrix(data_folder,ofile_path,map_dict=map_dict)
    
    return


if __name__ == '__main__':
    pass
#    writeTCGAExpressionMatrix(['LAML'],ignore_normal=False)