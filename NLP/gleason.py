from mind_minio_client import client
import pandas as pd
import numpy as np
from pandas import Series
import re
from io import BytesIO

obj = client.get_object('cdm-data','pathology/ddp_pathology_reports.tsv')
df_path = pd.read_csv(obj,sep='\t')

def extractGleason(s):
    splitstr = 'Gleason'
    if 'Gleason score' in s:
        splitstr = 'Gleason score'
    tokens = re.split(splitstr,s,flags=re.IGNORECASE)
    if len(tokens)>1:
        maxGleason = []
        for t in tokens[1:]:
            ppg = parsePostGleasonStr(t[:min(len(t),20)])
            if ppg==ppg:
                maxGleason+=[ppg]
        if len(maxGleason)>0:
            return max(maxGleason)
    return np.nan

def parsePostGleasonStr(s):
    if '+' in s:
        s = ''.join(s.split())
        tokens = s.split('+')
        if len(tokens)>1 and len(tokens[0])>0 and len(tokens[1])>0 and tokens[0][-1].isnumeric() and tokens[1][0].isnumeric():
            return int(tokens[0][-1])+int(tokens[1][0])
        #else:
        #    print(tokens[0])
        #    print(tokens[1])
    elif ':' in s:
        s = ''.join(s.split())
        tokens = s.split(':')
        if len(tokens)>0 and len(tokens[1])>0:
            if tokens[1][0].isnumeric():
                return int(tokens[1][0])
    return np.nan

df_path_gleason = df_path[df_path['path_prpt_p1'].fillna('').str.contains('Gleason',case=False)]
df_path_gleason['Gleason'] = df_path_gleason['path_prpt_p1'].apply(extractGleason)

bytedata=df_path_gleason[['MRN','Accession Number','Path Procedure Date','Gleason']].to_csv(index=False).encode('utf-8')
bufferdata=BytesIO(bytedata)

client.put_object('cdm-data', 'pathology/pathology_gleason_calls.csv',
                  data=bufferdata,
                  length=len(bytedata),content_type='application/csv')
