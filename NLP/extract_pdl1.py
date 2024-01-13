import pandas as pd
import sys
import re
from mind_minio_client import client
import numpy as np

obj = client.get_object('cdm-data','pathology/ddp_pathology_reports.tsv')
df_path = pd.read_csv(obj,sep='\t')

def getPDL1(report):
    for s in ['PD-L1 \(','PD-L1:','PD-L1|PDL-1|PDL1']:
        if len(re.split(s,report))>1:
            return grabPercent(report,s)
    return np.nan

def grabPercent(report,s):
    pstrings = re.split(s,report)

    for pstring in pstrings[1:]:
        pstring = pstring[:150]
        if '%' in pstring or 'of 100' in pstring:
            pstring = re.split('%|\(of 100',pstring)[0]
            if 'negative' in pstring.lower():
                return 'Negative'
            pstring = pstring[-8:]
            for block in [':','(','N',';']:
                if block in pstring:
                    pstring = pstring.split(block)[1]
            return re.sub("[^0-9<>-]", "",pstring)
        elif 'negative' in pstring.lower():
            return 'Negative'
        elif 'Positive (>=1)' in pstring:
            return 'Positive'
    return ''

def parse_pdl1(s):
    if not s==s:
        return np.nan
    elif s=='Negative':
        return 0
    elif s=='Positive':
        return 1
    elif '<1' in s:
        return 0
    elif '-' in s:
        tokens = s.split('-')
        if len(tokens[0])>0:
            return (float(tokens[0])+float(tokens[1]))/2
        return float(tokens[1])
    else:
        s2 = re.sub('[^0-9]','', s)
        if len(s2)>0:
            return int(s2)
        return 0
    
df_path['report'] = df_path['path_prpt_p1'].fillna('')+df_path['path_prpt_p2'].fillna('')
df_path['E1L3N'] = df_path['report'].str.contains('E1L3N')
df_path['SP-142'] = df_path['report'].str.contains('SP-142|SP142',regex=True)
df_path['SP-263'] = df_path['report'].str.contains('SP-263|SP263',regex=True)
df_path['PD-L1 status'] = df_path['report'].apply(getPDL1)
df_path['PD-L1'] = df_path['PD-L1 status'].apply(parse_pdl1)
df_path[~df_path['PD-L1 status'].isna() & (df_path['PD-L1 status'].fillna('').apply(len)>0)\
        | df_path['E1L3N'] | df_path['SP-142'] | df_path['SP-263']][['MRN','Path Report Type','Path Procedure Date','Accession Number','PD-L1','E1L3N','SP-142','SP-263']].to_csv('pdl1.csv',index=False)

