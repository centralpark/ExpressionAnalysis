"""
Created on Tue Sep 23 11:34:27 2014
@author: HSH
"""

import re
import pandas as pd

df = pd.io.parsers.read_table('/Users/HSH/Roche/Data/all15.txt',low_memory=False)

input_str = '''Consider run gene differential analysis for fusion PPP1R14B-DNAJC4 in GBM
Consider run gene differential analysis for fusion COL1A1-MFSD10 in SARC
Consider run gene differential analysis for fusion OLFML3-COL1A1 in SARC
Consider run gene differential analysis for fusion DCAF6-MPZL1 in BRCA
Less than 3 patients harboring fusion, cluster will not be computed.
Consider run gene differential analysis for fusion NOC4L-MMP17 in SKCM
Consider run gene differential analysis for fusion BIRC6-TTC27 in STAD
Consider run gene differential analysis for fusion ABCC1-MYH11 in LUSC
Consider run gene differential analysis for fusion CHGA-RYR1 in PCPG
Consider run gene differential analysis for fusion CTSB-HDLBP in LUSC
Consider run gene differential analysis for fusion TACC3-FGFR3 in BLCA
Consider run gene differential analysis for fusion DLK1-LDLR in ACC
Consider run gene differential analysis for fusion CTRB2-CELA3A in PAAD
Consider run gene differential analysis for fusion CHI3L1-MYOG in GBM
Consider run gene differential analysis for fusion HAND2-CHGB in PCPG
Consider run gene differential analysis for fusion DLK1-NDRG4 in ACC
Consider run gene differential analysis for fusion TPM4-KLF2 in LAML
Consider run gene differential analysis for fusion THSD4-LRRC49 in BRCA
Consider run gene differential analysis for fusion CDC42BPG-CAP1 in LUSC
Consider run gene differential analysis for fusion CTRB2-CELA2A in PAAD
Consider run gene differential analysis for fusion COL1A1-MFSD10 in BRCA
Consider run gene differential analysis for fusion MRPL48-DTX4 in BRCA
Consider run gene differential analysis for fusion FURIN-PMEL in SKCM
Consider run gene differential analysis for fusion SEMA4C-ANKRD39 in SKCM
Consider run gene differential analysis for fusion CELA3A-CTRC in PAAD
Consider run gene differential analysis for fusion CTSB-AP1G1 in SKCM
Consider run gene differential analysis for fusion TOMM40L-FGA in LIHC
Consider run gene differential analysis for fusion TTYH3-MAD1L1 in SKCM
Consider run gene differential analysis for fusion SEMA4C-ANKRD39 in BLCA
Consider run gene differential analysis for fusion AP2M1-POLR2A in HNSC
Consider run gene differential analysis for fusion CD74-PSAP in LUAD
Consider run gene differential analysis for fusion KAT6B-ADK in BRCA
Consider run gene differential analysis for fusion HIPK2-GLTSCR2 in COAD
Consider run gene differential analysis for fusion PLXND1-TMCC1 in BRCA
Consider run gene differential analysis for fusion FAT1-MTNR1A in HNSC
Consider run gene differential analysis for fusion NT5DC4-CKAP2L in ESCA
Consider run gene differential analysis for fusion CAGE1-SSR1 in BRCA
Consider run gene differential analysis for fusion PPP3CB-MUC1 in OV
Consider run gene differential analysis for fusion BCAS3-CBWD6 in LUSC
Consider run gene differential analysis for fusion ITIH3-C1R in LIHC
Consider run gene differential analysis for fusion EIF4E2-EFHD1 in BRCA
Consider run gene differential analysis for fusion DNAJB6-UBE3C in BRCA
Consider run gene differential analysis for fusion UBTF-MAML3 in PCPG
Consider run gene differential analysis for fusion NOS3-APOC1 in LIHC
Consider run gene differential analysis for fusion SORL1-TBCEL in OV
Consider run gene differential analysis for fusion CAGE1-SSR1 in COAD
Consider run gene differential analysis for fusion BRAF-SND1 in THCA
Consider run gene differential analysis for fusion COL2A1-COL9A2 in BRCA
Consider run gene differential analysis for fusion CAGE1-SSR1 in ESCA
Consider run gene differential analysis for fusion CHRNG-EIF4E2 in UCS
Consider run gene differential analysis for fusion CRYBA4-CRYBB1 in DLBC
Consider run gene differential analysis for fusion ZNHIT1-ACTB in OV
Consider run gene differential analysis for fusion CDC42BPG-CAP1 in BRCA
Consider run gene differential analysis for fusion EEF2-DAPK3 in GBM
Consider run gene differential analysis for fusion HNRNPA1-MUC2 in READ
Consider run gene differential analysis for fusion SRRT-SSRP1 in SKCM
Consider run gene differential analysis for fusion AGO2-PTK2 in OV
Consider run gene differential analysis for fusion FGB-SEPP1 in LIHC
Consider run gene differential analysis for fusion SPP1-TPP1 in KIRP
Consider run gene differential analysis for fusion KLHL14-CCDC178 in OV
Consider run gene differential analysis for fusion ARGLU1-EFNB2 in OV
Consider run gene differential analysis for fusion SEPT14-LANCL2 in GBM
Consider run gene differential analysis for fusion TBC1D9-TNRC18 in BRCA
Consider run gene differential analysis for fusion CHCHD3-PRF1 in KIRP
Consider run gene differential analysis for fusion XPO1-USP34 in OV
Consider run gene differential analysis for fusion WIZ-AKAP8L in SKCM
Consider run gene differential analysis for fusion ANK3-CCDC6 in BRCA
Consider run gene differential analysis for fusion ILF3-QTRT1 in LGG
Consider run gene differential analysis for fusion SLC22A16-METTL24 in UCEC
Consider run gene differential analysis for fusion CPB1-MYH9 in BRCA
Consider run gene differential analysis for fusion FAIM2-BCDIN3D in BRCA
Consider run gene differential analysis for fusion WBP2-UNC13D in GBM
Consider run gene differential analysis for fusion CELA3A-CELA2A in PAAD
Consider run gene differential analysis for fusion PRIMA1-KRT6C in LUSC
Consider run gene differential analysis for fusion COL1A1-MFSD10 in LUAD
Consider run gene differential analysis for fusion CTSB-HDLBP in SKCM
Consider run gene differential analysis for fusion PPP1CB-PLB1 in OV
Consider run gene differential analysis for fusion CTSB-CORO1C in SKCM
Consider run gene differential analysis for fusion COL1A1-MFSD10 in LUSC
Consider run gene differential analysis for fusion CDC42BPG-CAP1 in KIRC
Consider run gene differential analysis for fusion GNAS-DALRD3 in BRCA
Consider run gene differential analysis for fusion COL1A1-MFSD10 in SKCM
Consider run gene differential analysis for fusion ABCA10-CD38 in LUSC
Consider run gene differential analysis for fusion CHCHD3-PRF1 in SKCM
Consider run gene differential analysis for fusion SORL1-TECTA in OV
Consider run gene differential analysis for fusion KAT6B-ADK in OV
Consider run gene differential analysis for fusion RGS6-HDHD1 in SKCM
Consider run gene differential analysis for fusion USP22-MYH10 in BRCA
Consider run gene differential analysis for fusion GSX2-KRT10 in BRCA
Consider run gene differential analysis for fusion GSTM2-DHRS2 in BLCA
Consider run gene differential analysis for fusion KMT2A-ELL in LAML
Consider run gene differential analysis for fusion LHX6-NDUFA8 in CESC
Consider run gene differential analysis for fusion CAGE1-SSR1 in LUSC
Consider run gene differential analysis for fusion CHI3L1-ADORA1 in GBM
Consider run gene differential analysis for fusion NCL-NDUFS2 in SKCM
Consider run gene differential analysis for fusion TTC6-MIPOL1 in BRCA
Consider run gene differential analysis for fusion AOAH-ELMO1 in LUSC
Consider run gene differential analysis for fusion CELA3A-CELA2B in PAAD
Consider run gene differential analysis for fusion PUM1-SDC3 in BRCA
Consider run gene differential analysis for fusion NPDC1-CHGA in PCPG
Consider run gene differential analysis for fusion EPS8L2-TALDO1 in OV
Consider run gene differential analysis for fusion COL1A1-ITGA3 in PAAD
Consider run gene differential analysis for fusion VCL-ADK in LUSC
Consider run gene differential analysis for fusion PTPN12-RSBN1L in STAD
Consider run gene differential analysis for fusion LHX6-NDUFA8 in BLCA
Consider run gene differential analysis for fusion VCL-ADK in STAD
Consider run gene differential analysis for fusion RNF10-PSAP in BRCA
Consider run gene differential analysis for fusion CHGB-SRRM2 in PCPG
Consider run gene differential analysis for fusion BCAS3-CBWD1 in LUSC
Consider run gene differential analysis for fusion BUB1B-EIF2AK4 in BRCA'''

lines = input_str.split('\n')
result = []
for line in lines:
    matchObj = re.match('Consider run gene differential analysis for fusion (.*)-(.*) in (.*)',line.rstrip())
    if not matchObj:
        continue
    match_row = df.ix[(df.Tumor==matchObj.group(3)) & (df.gene5p==matchObj.group(1)) & (df.gene3p==matchObj.group(2)) & (df.IsReal==1),['Tumor','gene5p','gene3p','Fusions']]    
    result.append(int(match_row.index)+1)
str(result)