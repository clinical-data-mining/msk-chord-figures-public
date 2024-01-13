from mind_minio_client import client
import pandas as pd
import numpy as np
from pandas import Series
import re
from io import BytesIO

obj = client.get_object('cdm-data','pathology/ddp_pathology_reports.tsv')
df_path = pd.read_csv(obj,sep='\t')

''' example:
     MLH1:     staining absent in tumor
     PMS2:     staining absent in tumor
     MSH2:     staining present in tumor
     MSH6:     staining present in tumor'''

def extractMMR(s):
    factors = ['MLH1','PMS2','MSH2','MSH6']
    for f in factors:
        if f in s:
            tokens = s.split(f)
            for i, token in enumerate(tokens[1:]):
                statement = '--'+token[:50] +'***'+ tokens[i][-50:]+'--'
                absentLoc = max(statement.lower().find('absent'),
                                statement.lower().find('loss'),
                               statement.lower().find('lost'))
                presentLoc= statement.lower().find('present')
                if absentLoc>=0 and presentLoc>=0:
                    if absentLoc<presentLoc:
                        return True
                if absentLoc>=0 and presentLoc<0:
                        return True
    return False

df_mmr = df_path[df_path['path_prpt_p1'].fillna('').str.contains('MLH1|PMS2|MSH2|MSH6',regex=True,case=False)]
df_mmr['MMR_ABSENT'] = df_mmr['path_prpt_p1'].apply(extractMMR)
df_mmr = df_mmr[~df_mmr['Accession Number'].str.contains('M')]
df_mmr[['MRN','Accession Number','Path Procedure Date','MMR_ABSENT']].to_csv('MMR.csv',index=False)

bytedata=df_mmr[['MRN','Accession Number','Path Procedure Date','MMR_ABSENT']].to_csv(index=False).encode('utf-8')
bufferdata=BytesIO(bytedata)
client.put_object('cdm-data', 'pathology/pathology_mmr_calls.csv',
                  data=bufferdata,
                  length=len(bytedata),content_type='application/csv')
