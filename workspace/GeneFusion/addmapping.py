"""
Updating the ccle mapping file.

Created on Mon Aug 18 11:53:54 2014
@author: HSH
"""

import pandas as pd
import xml.etree.ElementTree as etree
import urllib2,re,os

mapping = pd.io.parsers.read_table('/Users/HSH/Roche/Data/connect/ccle_tcga_mapping.txt',
                                   low_memory=False)
codeTable = pd.io.parsers.read_csv('/Users/HSH/Roche/Data/sampleType.csv',dtype={'Code':'string'})

# The COAD has updated expression data
data_folder = '/Users/HSH/Roche/Data/tcga_data_20140811145356/COAD'
data_files = os.listdir(data_folder)

dicts = []
for df in data_files:
    aliquot_id = re.search('\w{8}-\w{4}-\w{4}-\w{4}-\w{12}',df).group()
    url = 'https://cghub.ucsc.edu/cghub/metadata/analysisFull?aliquot_id=' + aliquot_id
    f = urllib2.urlopen(url)
    tree = etree.parse(f)
    root = tree.getroot()
    barcode = root.findall('./Result/legacy_sample_id')[0].text
    disease = root.findall('./Result/disease_abbr')[0].text
    sample_type = root.findall('./Result/sample_type')[0].text
    sample_type = codeTable.ix[codeTable.Code==sample_type,2].values[0]
    analysis_id = root.findall('./Result/analysis_id')[0].text
    participant_id = root.findall('./Result/participant_id')[0].text
    sample_id = root.findall('./Result/sample_id')[0].text
    d = {'study':'TCGA','barcode':barcode,'disease':disease,'sample_type':sample_type,
         'analysis_id':analysis_id,'aliquot_id':aliquot_id,'participant_id':participant_id,
         'sample_id':sample_id}
    dicts.append(d)

mappingNew = mapping.append(dicts,ignore_index=True)
mappingNew.to_csv('/Users/HSH/Roche/Data/connect/ccle_tcga_mapping_v1.txt',index=False,
                  sep='\t')