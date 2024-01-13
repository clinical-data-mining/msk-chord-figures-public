import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
import numpy as np
import scipy.stats as stats
from lifelines import KaplanMeierFitter, BreslowFlemingHarringtonFitter
from lifelines.statistics import proportional_hazard_test
from datetime import datetime, timedelta, time, date
from lifelines.plotting import plot_lifetimes
from lifelines import CoxPHFitter
from lifelines.utils import median_survival_times
import matplotlib.ticker as mtick
import sys
import os

# print num at risk patients at times
def move_spines(ax, sides, dists):
    for side, dist in zip(sides, dists):
        ax.spines[side].set_position(("axes", dist))
    return ax

def remove_spines(ax, sides):
    for side in sides:
        ax.spines[side].set_visible(False)
    return ax

def n_at_risk(df,ax,label,xshift,yshift=10,rt_censor='stop',lt_trunc='start'):
    #df = df.sort_values(by='stop',ascending=False).drop_duplicates(subset='MRN')
    
    # Create another axes where we can put size ticks
    ax2 = plt.twiny(ax=ax)
    # Move the ticks below existing axes
    ax_height = (ax.get_position().y1 - ax.get_position().y0)*yshift
    ax2_ypos = -1 / ax_height

    move_spines(ax2, ["bottom"], [ax2_ypos])
    # Hide all fluff
    remove_spines(ax2, ["top", "right", "bottom", "left"])
    # Set ticks and labels on bottom
    ax2.xaxis.tick_bottom()
    # Set limit
    min_time, max_time = ax.get_xlim()
    ax2.set_xlim(min_time, max_time)
    ax2.set_xticks(ax.get_xticks())
    
    ticklabels = []
    for i in ax.get_xticks():
        ticklabels.append(sum((df[lt_trunc]/daysinmo<=i) & (df[rt_censor]/daysinmo>i)))
        
    ax2.set_xticklabels(ticklabels, ha="center")
    ax2.grid(False)
    ax2.text(xshift, ax2_ypos-0.13, label, fontsize=11)
    return ax

def calculateHR(df1,df2,fancy=True,correction_cols = []): #['Male','White','IS_PRIMARY','CVR_TMB_SCORE','TUMOR_PURITY']
    # ---calculate hazard ratio---
    columns = ['dead','start','stop','var'] + correction_cols #correct for certain variables
    df1['var'] = True
    df2['var'] = False

    df_joint = pd.concat([df1[columns],df2[columns]])
    df_joint = df_joint.replace([np.inf, -np.inf], np.nan).dropna()
    df_joint.loc[df_joint['start']<0,'start'] = 0
    df_joint.loc[df_joint['stop']<=df_joint['start'],'stop'] = df_joint[df_joint['stop']<=df_joint['start']]['start']+1

    cph = CoxPHFitter().fit(df_joint, event_col="dead", entry_col="start", duration_col="stop", show_progress=True)

    print(cph.log_likelihood_ratio_test())
    
    HR = (round(cph.hazard_ratios_['var'], 2))
    lowerCI = (round(np.exp(cph.confidence_intervals_.loc['var','95% lower-bound']),2))
    upperCI = (round(np.exp(cph.confidence_intervals_.loc['var','95% upper-bound']),2))
    if fancy:
        return 'HR = '+str(HR)+' (95% CI, '+str(lowerCI)+'-'+str(upperCI)+')'
    return (HR,lowerCI,upperCI)

daysinmo = 30.44

def generateFigure(df_master, s1,s2,label1,label2,**kwargs):
    
    xshift=-500
    figname='test.svg'
    targeted = False
    df = df_master
    ci_show = False
    null_nontargeted = True
    idcol='MRN'
    if 'idcol' in kwargs:
        idcol = kwargs['idcol']
    if 'xshift' in kwargs:
        xshift = kwargs['xshift']
    if 'figname' in kwargs:
        figname = kwargs['figname']
    if 'targeted' in kwargs:
        targeted= kwargs['targeted']
    if 'df' in kwargs:
        df = kwargs['df']
    if 'ci_show' in kwargs:
        ci_show = kwargs['ci_show']
    if 'null_nontargeted' in kwargs:
        null_nontargeted = kwargs['null_nontargeted']
    
    fig, ax = plt.subplots()
    
    xshift=xshift/daysinmo

    # ---KM curves---
    df1 = df[df[idcol].isin(s1)]
    if targeted: # select out only when treated with targeted therapy
        df1_nontargeted = df1[~df1['targeted']]
        df1 = df1[df1['targeted']]
    kmf1 = KaplanMeierFitter()
    kmf1.fit(df1["stop"]/daysinmo, event_observed=df1["dead"], entry=df1['start']/daysinmo, label=label1)
    ax = kmf1.plot(color='blue',linestyle='dashed',ci_show=ci_show,show_censors=(not ci_show),censor_styles={'ms': 3, 'marker': '|'})

    df2 = df[df[idcol].isin(s2)] 
    if targeted and null_nontargeted: #append orphan intervals pre-targeted therapy to non-targeted dataframe
        df2 = pd.concat([df1_nontargeted,df2])
    kmf2 = KaplanMeierFitter()
    kmf2.fit(df2["stop"]/daysinmo, event_observed=df2["dead"], entry=df2['start']/daysinmo, label=label2)
    kmf2.plot(ax=ax,color='black',ci_show=ci_show,show_censors=(not ci_show),censor_styles={'ms': 3, 'marker': '|'})

    ax.set_xlim([0.0, 40])
    ax.set_xlabel("Months",fontsize=12)
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
    ax.set_ylabel("Overall survival",fontsize=12)
    ax.grid(False)

    ax.text(xshift,-0.2,'Number at risk',fontsize=11)
    n_at_risk(df1,ax,label1,xshift,8)
    n_at_risk(df2,ax,label2,xshift,5)

    HR_str = calculateHR(df1,df2)

    ax.text(40/daysinmo,0,HR_str,fontsize=11)

    fig.savefig(figname, format='svg')
    print(label1+': '+str(kmf1.median_survival_time_))
    print(median_survival_times(kmf1.confidence_interval_))
    print(label2+': '+str(kmf2.median_survival_time_))
    print(median_survival_times(kmf2.confidence_interval_))
    
    return ax

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

def myRound(f):
    if (f==float('inf')) or (f==float('-inf')):
        return f
    return int(f)
        
def generateKMFSummary(df_master, s1,s2,label,targeted=False,**kwargs):
    
    df = df_master
    correctionCols = []
    idcol = 'MRN'
    if 'idcol' in kwargs:
        idcol = kwargs['idcol']
    if 'df' in kwargs:
        df = kwargs['df']
    if 'correctionCols' in kwargs:
        correctionCols = kwargs['correctionCols']
    df1 = df[df[idcol].isin(s1)]
  #  if targeted:
  #      df1_nontargeted = df1[~df1['targeted']]
  #      df1 = df1[df1['targeted']]
    try:
        kmf1 = KaplanMeierFitter()
        kmf1.fit(df1["stop"], event_observed=df1["dead"], entry=df1['start'], label='x')   
    except:
        kmf1 = BreslowFlemingHarringtonFitter()
        kmf1.fit(df1["stop"], event_observed=df1["dead"], entry=df1['start'], label='x') 
        
    df2 = df[df[idcol].isin(s2)]
  #  if targeted:
  #      df2 = pd.concat([df1_nontargeted,df2])
    kmf2 = KaplanMeierFitter()
    kmf2.fit(df2["stop"], event_observed=df2["dead"], entry=df2['start'], label='y')
    
    ci1 = median_survival_times(kmf1.confidence_interval_)/daysinmo
    str1 = (str(myRound(kmf1.median_survival_time_/daysinmo))+' ('+str(myRound(ci1.loc[0.5,'x_lower_0.95']))+'-'+str(myRound(ci1.loc[0.5,'x_upper_0.95']))+')')
    
    ci2 = median_survival_times(kmf2.confidence_interval_)/daysinmo
    str2 = (str(myRound(kmf2.median_survival_time_/daysinmo))+' ('+str(myRound(ci2.loc[0.5,'y_lower_0.95']))+'-'+str(myRound(ci2.loc[0.5,'y_upper_0.95']))+')')
    
    with HiddenPrints():
        hr = calculateHR(df1,df2,False,correctionCols)
    
    print(label+','+str(len(s1))+','+str1+','+str(len(s2))+','+str2)
    return (label,hr)


def myForest(HRs,rLabel,lLabel,name='forest.svg',figsize=(2.4,4)):
    fig, ax = plt.subplots(figsize=figsize)
    #hfont = {'fontname':'Helvetica'}

    if isinstance(HRs, pd.DataFrame):
        dfHR=HRs
    else:
        dfHR = pd.DataFrame(HRs.to_list(),columns=['HR','-95%CI','+95%CI'],index=HRs.index).iloc[::-1]
    #transform to log scale (JAMA requirement)
    dfHR['HR'] = dfHR['HR'].apply(np.log2)
    dfHR['upper'] = dfHR['+95%CI'].apply(np.log2) - dfHR['HR']
    dfHR['lower'] = dfHR['HR'] - dfHR['-95%CI'].apply(np.log2)
    ax.errorbar(dfHR['HR'],dfHR.index,xerr=dfHR[['lower','upper']].values.T,
                elinewidth=1,markeredgecolor='w',markerfacecolor='k',markersize=5,fmt='o',color='k')
    ax.grid(False)
    ax.axvline(x=0,color='k',linestyle='--')
    ax.set_xticklabels(['%.1f' % 2**i for i in ax.get_xticks()])

    ax.set_yticklabels(dfHR.index, ha = 'left')#,**hfont)
    plt.draw()
    yax = ax.get_yaxis()
    pad = max(T.label.get_window_extent().width for T in yax.majorTicks)
    yax.set_tick_params(pad=pad)

    (left,right) = ax.get_xlim()
    ax.text(left,-len(dfHR)/5,lLabel,ha='center')
    ax.text(right,-len(dfHR)/5,rLabel,ha='center')
    remove_spines(ax, ["top", "right", "bottom", "left"])

    dfHR['hr_str'] = dfHR['HR'].apply(lambda s: '%.2f' % 2**s).astype(str)
    dfHR['hr_str'] +=' ('+dfHR['-95%CI'].apply(lambda s: '%.2f' % s).astype(str)
    dfHR['hr_str'] +=','+dfHR['+95%CI'].apply(lambda s: '%.2f' % s).astype(str)+')'

    for i in range(len(dfHR)):
        ax.text(right,i,dfHR.loc[dfHR.index[i],'hr_str'])

    plt.draw()
    fig.savefig(name, format=name.split('.')[-1])
    
    
def simpleP(p):
    pstring = 'p='+str(round(p,3))
    if p>0.01:
        pstring = 'p='+str(round(p,2))
    if p>0.99:
        pstring = 'p>0.99'
    if p<0.001:
        pstring = 'p<0.001'
    return pstring

def myViolin(df,col1,col2,ax):
    df = df[[col1,col2]].dropna()
    sns.violinplot(x=col1,y=col2,data=df,ax=ax)
    ax.set(xlabel=col1)
    
    l = ax.get_xlabel()
    ax.set_xlabel(l, fontsize=12)
    
    l = ax.get_ylabel()
    ax.set_ylabel(l, fontsize=12)
    
    if df[col1].dtypes.name == 'bool':
        mask = df[col1]
        r, p = stats.mannwhitneyu(df[mask][col2],df[~mask][col2])
        ax.text(0.5,df[col2].max(),simpleP(p),horizontalalignment='center')
    ax.grid(False)
    
def myHist(df,col1,col2,ax,labels=[],nbins=10):
    df = df[[col1,col2]].dropna()
    for i in range(len(df[col1].unique())):
        c = sorted(df[col1].unique())[i]
        temp = df[df[col1]==c][col2]
        b = max(temp)-min(temp)
        if not isinstance(b,int):
            b = nbins# default
        weights = np.ones_like(temp) / len(temp)
        
        histtype='step'
        alpha=1
        rwidth=1
        if c==sorted(df[col1].unique())[0]:
            alpha=0.4
            histtype='bar'
            rwidth=0.95
        
        label=col1+' '+str(c)
        if len(labels)>0:
            label=labels[i]
        (z,x,y) = ax.hist(temp,label=label,weights=weights,
                bins=b,range=(min(temp),max(temp)),alpha=alpha,histtype=histtype,rwidth=rwidth)
        
        
    ax.set_xlabel(col2, fontsize=12)
    ax.set_ylabel('Proportion of samples', fontsize=12)
    
    if df[col1].dtypes.name == 'bool':
        mask = df[col1]
        r, p = stats.mannwhitneyu(df[mask][col2],df[~mask][col2])
        ax.text(np.median(x),max(z)/2,simpleP(p),horizontalalignment='center')
    ax.grid(False)
    ax.legend()
    
def myScatter(df,col1,col2,ax,fsize=12):
    df = df[[col1,col2]].dropna()
    sns.scatterplot(df[col1],df[col2],ax=ax,alpha=0.3)
    
    l = ax.get_xlabel()
    ax.set_xlabel(l, fontsize=12)
    
    l = ax.get_ylabel()
    ax.set_ylabel(l, fontsize=12)
    
    r, p = stats.spearmanr(df[col1],df[col2])
    dig=3
    if p>0.01:
        dig=2
    pstring = '='+str(round(p,dig))
    if p>0.99:
        pstring = '>0.99'
    if p<0.001:
        pstring = '<0.001'
    ax.text(df[col1].max(),df[col2].max(),'R='+str(round(r,2))+', p'+pstring,horizontalalignment='right',fontsize=fsize)
    ax.grid(False)
    
# generate Table 1 data
def stageint2roman(i):
    d = {1:'I',2:'II',3:'III',4:'IV'}
    if i in d:
        return d[i]
    return '-'
    
def mySummary(df,name):
    summ = []
    for c in df.columns:
        s = df[c]
        if 'HAS_'+c in df.columns:
            if sum(df['HAS_'+c].astype(int))>0:
                s = df.loc[df['HAS_'+c].astype(bool),c]
        if len(s.value_counts())<=2: #s.fillna(False).dtype=='bool' or 
            summ+=[str(sum(s.fillna(False)))+' ('+str(round(100*s.fillna(False).mean()))+'%)']
        else:
            # check if int or float
            tempsumm=''
            s = s.dropna()
            if c=='STAGE_DX_INT':
                x = s.value_counts().to_dict()
                tempsumm = [stageint2roman(int(k))+':'+str(x[k]) for k in sorted(x)]
            elif np.array_equal(s, s.astype(int)) or c=='AGE':
                tempsumm = str(int(s.median()))+\
                ' ('+str(int(s.quantile(0.05)))+\
                ','+str(int(s.quantile(0.25)))+\
                ','+str(int(s.quantile(0.75)))+\
                ','+str(int(s.quantile(0.95)))+')'
            else:
                tempsumm = str(round(s.median(),2))+\
                ' ('+str(round(s.quantile(0.05),2))+\
                ','+str(round(s.quantile(0.25),2))+\
                ','+str(round(s.quantile(0.75),2))+\
                ','+str(round(s.quantile(0.95),2))+')'
            summ+=[tempsumm]
    return pd.Series(summ,index=df.columns,name=name+' N='+str(len(df)))