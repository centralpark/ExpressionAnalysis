'''
Created on Jul 8, 2014
@author: HSH
'''

import re

def processSIF(ifile_path,ofile_path,map_file_path='/home/staff/hes14/Roche/Data/uniprot_sprot_human.dat'):
    '''Modify the sif file from Pathway Commons properly for subsequent analysis'''
    map_dict = {}
    f_map = open(map_file_path)
    content = f_map.read()
    content_sep = content.split('//\n')
    for entry in content_sep:
        record = entry.split('\n')
        list_ac = []
        gene_name = ''
        for line in record:
            if line.startswith('AC'):
                list_ac += map(lambda x:x.rstrip(';'),line.split()[1:])
            if line.startswith('GN'):
                matchObj = re.search('Name=([-\w]+);',line)
                if matchObj:
                    gene_name = matchObj.group(1)
        if (not list_ac) or (not gene_name):
            continue
        for ac in list_ac:
            map_dict[ac] = gene_name
    f_map.close()
     
    f = open(ifile_path)
    ofile = open(ofile_path,'w')
    for line in f.readlines():
        [entity_a,relate,entity_b] = line.rstrip().split('\t')
        if entity_a.startswith('http://identifiers.org/uniprot'):
            uniprot_id_a = entity_a.split('/')[-1]
        else:
            continue
        if entity_b.startswith('http://identifiers.org/uniprot'):
            uniprot_id_b = entity_b.split('/')[-1]
        else:
            continue
        try:
            line_out = '\t'.join([map_dict[uniprot_id_a],relate,map_dict[uniprot_id_b]])
        except KeyError:
            print uniprot_id_a+'\t'+uniprot_id_b
            continue
        ofile.write(line_out+'\n')
    f.close()
    ofile.close()
    

def constructNetwork(gene_list,sif_path,ofile_path,
                     map_file_path = '/home/staff/hes14/Roche/Data/uniprot_sprot_human.dat'):
    # map_dict maps Uniprot accession ID to gene symbol
    map_dict = {}
    f_map = open(map_file_path)
    content = f_map.read()
    content_sep = content.split('//\n')
    for entry in content_sep:
        record = entry.split('\n')
        list_ac = []
        gene_name = ''
        for line in record:
            if line.startswith('AC'):
                list_ac += map(lambda x:x.rstrip(';'),line.split()[1:])
            if line.startswith('GN'):
                matchObj = re.search('Name=([-\w]+);',line)
                if matchObj:
                    gene_name = matchObj.group(1)
        if (not list_ac) or (not gene_name):
            continue
        for ac in list_ac:
            map_dict[ac] = gene_name
    f_map.close()
    
    for i,val in enumerate(gene_list):
        gene_list[i] = val.upper()
    
    f = open(sif_path)
    ofile = open(ofile_path,'w')
    for line in f.readlines():
        [entity_a,relate,entity_b] = line.rstrip().split('\t')
        if entity_a.startswith('http://identifiers.org/uniprot'):
            uniprot_id_a = entity_a.split('/')[-1]
        else:
            continue
        if entity_b.startswith('http://identifiers.org/uniprot'):
            uniprot_id_b = entity_b.split('/')[-1]
        else:
            continue
        if uniprot_id_a in map_dict:
            gene_a = map_dict[uniprot_id_a]
        else:
            continue
        if uniprot_id_b in map_dict:
            gene_b = map_dict[uniprot_id_b]
        else:
            continue
        if (gene_a in gene_list) and (gene_b in gene_list):
            line_out = '\t'.join([gene_a,relate,gene_b])
            ofile.write(line_out+'\n')
        else:
            continue
    f.close()
    ofile.close()
    
    
    
if __name__ == '__main__':    
    import sys
    pathway_name = sys.argv[1]
    out_fname = sys.argv[2]    
    f = open('/home/staff/hes14/Roche/Data/ReactomePathways.gmt')
    for line in f.readlines():
        if line.startswith(pathway_name):
            break
    f.close()
    gene_list = line.rstrip().split('\t')[2:]
    constructNetwork(gene_list,'/home/staff/hes14/Roche/Data/Pathway Commons.4.All.BINARY_SIF.tsv',
                     out_fname)       