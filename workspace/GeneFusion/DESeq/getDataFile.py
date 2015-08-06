'''
This module get basic data or information and write (optionally) to file
Created on Jul 3, 2014
@author: HSH
'''

from Bio import Entrez
import imp

Entrez.email = 'seas363@gmail.com'

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
                result.append(probe_id)
                print probe_id
                continue
            title = line_split[2]
            if 'pseudogene' in title:
                result.append(probe_id)
        f.close()
        
        return result
    
    else:
        f.close()
        raise Exception('Mode not supported!')


def getAssayGeneAnnot(gene_info_path='/Users/HSH/Roche/Data/Homo_sapiens.gene_info',
                      gct_path='/Users/HSH/Roche/Data/expression/ACC.gct',
                      ofile_path='/Users/HSH/Roche/Data/expression/gene_annotation_v1.txt'):
    '''Write gene_annotation_v1.txt, probe annotation for RNASeqV2 gene expression.'''
    # Homo sapiens gene infomation
    gene_info = {}
    f_info = open(gene_info_path)
    f_info.readline()
    for line in f_info.readlines():
        line_split = line.rstrip().split('\t')
        [gene_id,symbol,gene_type] = [line_split[i] for i in [1,2,9]]
        gene_info[gene_id] = [symbol,gene_type]
    f_info.close()
    
    # list of genes used in expression assay
    f_gct = open(gct_path)
    f_gct.readline()
    probe_list = [line.rstrip().split('\t')[0] for line in f_gct.readlines()]
    f_gct.close()
    probe_list.sort(key=int)
    
    # write the file
    ofile = open(ofile_path,'w')
    header = '\t'.join(['Probe Set ID','Gene Symbol','Gene Type'])
    ofile.write(header+'\n')
    # some of the gene information are not available (obsolete)
    probe_obsolete = []
    for p in probe_list:
        if p not in gene_info:
            probe_obsolete.append(p)
            continue
        line_out = '\t'.join([p]+gene_info[p])
        ofile.write(line_out+'\n')
        
    for p in probe_obsolete:
        handle = Entrez.efetch(db='gene',id=p,retmode='xml')
        record = Entrez.read(handle)
        handle.close()
        symbol = record[0]['Entrezgene_gene']['Gene-ref']['Gene-ref_locus']
        gene_type = record[0]['Entrezgene_type']
        # sepearate from identical symbol genes
        line_out = '\t'.join([p,symbol+'|'+p,gene_type])
        ofile.write(line_out+'\n')
    
    ofile.close()
    return


def getAssayGeneAnnot_v1(ifile_path='/Users/HSH/Roche/Data/TCGA.Sept2010.09202010.gene.gaf',
                       ofile_path='/Users/HSH/Roche/Data/expression/gene_annotation_v1.txt',
                       gct_path='/Users/HSH/Roche/Data/expression/THCA.gct'):
    '''Write gene_annotation_v1.txt, probe annotation for RNASeqV2 gene expression.'''
    f = open(ifile_path)
    f.readline()
    map_dict = {}
    for line in f.readlines():
        line_split = line.rstrip().split('\t')
        feature_type = line_split[2]
        if feature_type != 'gene':
            continue
        featureID = line_split[1]
        symbol,probe_id = featureID.split('|')[0:2]
        map_dict[probe_id] = symbol
    f.close()
    
    # list of genes used in expression assay
    f = open(gct_path)
    f.readline()
    probe_list = [line.rstrip().split('\t')[0] for line in f.readlines()]
    f.close()
    probe_list.sort(key=int)
    
    ofile = open(ofile_path,'w')
    header = '\t'.join(['Probe Set ID','Gene Symbol','Gene Title','Gene Type'])
    ofile.write(header+'\n')
    for p in probe_list:
        symbol_0 = map_dict[p]
        if symbol_0 == '?':
            symbol = ''
            title = ''
            gene_type = ''
        else:
            handle = Entrez.esearch(db='gene',term=symbol_0+' AND human[ORGN]',sort='relevance')
            record = Entrez.read(handle)
            handle.close()
            if record['Count'] == '0':
                symbol = ''
                title = ''
                gene_type = ''
                print p + '\t' + symbol_0
            else:
                handle = Entrez.efetch(db='gene',id=record['IdList'][0],retmode='xml')
                record = Entrez.read(handle)
                handle.close()
                symbol = record[0]['Entrezgene_gene']['Gene-ref']['Gene-ref_locus']
                title = record[0]['Entrezgene_gene']['Gene-ref']['Gene-ref_desc']
                gene_type = record[0]['Entrezgene_type']
        line_out = '\t'.join([p,symbol,title,gene_type])
        ofile.write(line_out+'\n')
    ofile.close()
    
    notice = imp.load_source('misc','/Users/HSH/Roche/workspace/GeneFusion/misc.py')
    notice.finishNotice()
    return
        

def getTCGAAnnotation(ifile_path='/Users/HSH/Roche/Data/TCGA.Sept2010.09202010.gaf',
                       ofile_path='/Users/HSH/Roche/Data/TCGA.Sept2010.09202010.gene.gaf'):
    f = open(ifile_path)
    ofile = open(ofile_path,'w')
    ofile.write(f.readline())
    for line in f.readlines():
        line_split = line.rstrip().split('\t')
        feature_type = line_split[2]
        if feature_type != 'gene':
            continue
        ofile.write(line)
    f.close()
    ofile.close()
    
    return

        
if __name__ == '__main__':
    getAssayGeneAnnot()