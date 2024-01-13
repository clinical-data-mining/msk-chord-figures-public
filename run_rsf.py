import pandas as pd
import numpy as np

from sksurv.datasets import load_gbsg2
from sksurv.preprocessing import OneHotEncoder
from sksurv.ensemble import RandomSurvivalForest
from sklearn.model_selection import KFold, RepeatedKFold

# harmonize treatment columns

# categorize columns
demographiccols = ['AGE','MALE','WHITE','ASIAN','BLACK','SMOKER','SEQUENCING_YEAR']

txhxcols = set(['ANY_PRIOR_TX'])

tumorregcols = ['STAGE 1','STAGE 2','STAGE 3','STAGE 4',
                'STAGE_IV_DX','STAGE_I-III_NOPROG','STAGE_I-III_PROG','entry','progressed']

metcols = ['DMETS_DX_ADRENAL',
 'DMETS_DX_BONE',
 'DMETS_DX_BRAIN',
 'DMETS_DX_LIVER',
 'DMETS_DX_LUNG',
 'DMETS_DX_LYMPH',
 'DMETS_DX_PLEURA',
 'DMETS_DX_OTHER']

pathcols = ['Gleason','HAS_Gleason','ADENOCARCINOMA','SQUAMOUS','PDL1','HAS_PDL1','HR','HER2',
           'RECTAL','ASCENDING','CECUM','NONADENOCARCINOMA','MUCINOUS',
            'MSI_OR_dMMR','HAS_MSI_OR_dMMR'] 

labcols = ['MAX_CA15-3','CA15-3','HAS_CA15-3',
          'MAX_CEA','CEA','HAS_CEA',
          'MAX_PSA','PSA','HAS_PSA',
          'MAX_CA19-9','CA19-9','HAS_CA19-9']

common_genes = ['KRAS', 'HRAS', 'RET', 'MET', 'GNAQ', 
                'PTEN', 'KIT', 'EGFR', 'FGFR3', 'PDGFRA', 
                'ERBB2', 'TP53', 'NRAS', 'NOTCH1', 'GNA11', 
                'CTNNB1', 'FGFR2', 'PIK3CA', 'IDH1', 'BRAF', 
                'FGFR1', 'ALK', 'AKT1']

cancers_of_interest = ['nsclc','brca','prostate','crc','panc']

# load 5 files into a dictionary of <cancer type>:<pandas dataframe>
cancers_of_interest = ['nsclc','brca','prostate','crc','panc']
cancer2df_master_current_tx = {}
for c in cancers_of_interest:
    cancer2df_master_current_tx[c] = pd.read_csv('data/'+c+'_dx_1st_seq_OS.csv')
    columns_noID = list(cancer2df_master_current_tx[c].columns)[1:]
    txhxcols_c = cancer2df_master_current_tx[c].columns[cancer2df_master_current_tx[c].columns.str.contains('ANY_CURRENT_|ANY_PREV_',regex=True)]
    txhxcols = txhxcols.union(set(txhxcols_c))


def runRF(df_train,df_tests,cois,endsuffix='',**kwargs):

    cois = list(set(df_train.columns).intersection(set(cois)))
    if len(cois)>0:
        X_train = df_train[cois].fillna(0).astype(int)
        y_train = df_train[['dead'+endsuffix,'stop'+endsuffix]].apply(tuple, axis=1).values.tolist()
        y_train = np.array(y_train, dtype=[('Status', '?'), ('Survival_in_days', '<f8')])

        random_state = 20
        rsf = RandomSurvivalForest(n_estimators=1000,
                                   min_samples_split=10,
                                   min_samples_leaf=15,
                                   n_jobs=-1,
                                   random_state=random_state)
        rsf.fit(X_train, y_train)
        scores = []
        for df_test in df_tests:
            X_test = df_test[cois].fillna(0).astype(int)
            y_test = df_test[['dead'+endsuffix,'stop'+endsuffix]].apply(tuple, axis=1).values.tolist()
            y_test = np.array(y_test, dtype=[('Status', '?'), ('Survival_in_days', '<f8')])
            try:
                scores += [rsf.score(X_test, y_test)]
            except:
                scores += [np.nan]
        return scores
    return [np.nan]*len(df_tests)

txhxcols = list(txhxcols)
all_cols = demographiccols+txhxcols+tumorregcols+metcols+pathcols+labcols+common_genes

# note: uncomment or comment items in variable_list and labels to include or exclude more categories of variables
variable_list = [demographiccols,
#                 txhxcols,
#                 tumorregcols,
#                 metcols,
#                 pathcols,
#                 labcols,
                 common_genes]
#                all_cols]

# for ablation studies
#variable_list = [set(all_cols)-set(v) for v in variable_list]

labels = ['Demographics',
#          'Treatment',
#          'Stage/Progression',
#          'Other Met Sites',
#          'Pathology',
#          'Tumor Markers',
          'Genomics']
#         'All']

scores = pd.DataFrame(columns=labels)

for c in ['nsclc','panc']: #testing on 2/5 but can substitute for cancers_of_interest
    print(c)
    for i, variables in enumerate(variable_list):
        print(variables)
        dftemp = cancer2df_master_current_tx[c]
        dftemp['stop_nonlt'] = dftemp['stop']-dftemp['entry']
        dftemp['dead_nonlt']=dftemp['dead']
        dftemp = dftemp[dftemp['stop_nonlt']>=0].sample(frac=1).reset_index()
        kf = KFold(n_splits=5)# Define the split - into n_splits folds
        kf.get_n_splits(dftemp) # returns the number of splitting iterations in the cross-validator 

        for f, (train_index, test_index) in enumerate(kf.split(dftemp)):
            print(f'fold: {f}')
            train = dftemp.loc[train_index]
            val = dftemp.loc[test_index]

            val_list = [val,val[val['STAGE_IV_DX'].astype(bool)],val[val['STAGE_I-III_NOPROG'].astype(bool)]]
            print([len(a) for a in val_list])
            scorelist = (runRF(train,val_list,variables,'_nonlt'))
            prefix='fold'+str(f)+'_'
            scores.loc[prefix+c,labels[i]] = scorelist[0]
            scores.loc[prefix+'stage IV '+c,labels[i]] = scorelist[1]
            scores.loc[prefix+'stage I-III '+c,labels[i]] = scorelist[2]

scores.to_csv('test_rsf_output.csv')
