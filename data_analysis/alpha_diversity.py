#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 10:36:26 2024

#Relative abundance calculation and alpha diversity

@author: ekateria
"""

import pandas as pd
import numpy as np
import pyreadr as pyr
import math
import seaborn as sb
import scipy.stats as stats
from itertools import combinations as comb
from statsmodels.formula.api import ols



from statsmodels.stats.multitest import multipletests
#import numpy as np

wdir='PATH_TO_MANUS_FOLDER' #set working directory
prefix='MAGs' #prefix to be used for the output files
#prefix='vOTUs' #prefix to be used for the output files
#prefix='pOTUs' #prefix to be used for the output files

#-----------------------------------------------------------------------------------------------------------
#Load metadata        

####NB that beforeBL variable codes for a/b use 4 months prior to baseline, not if it was at ANY time before baseline  
               
meta=pd.read_csv('/'.join([wdir,'participant_data/screening_data.tsv']), sep='\t')
meta['beforeBL'] = pd.Categorical(meta['beforeBL'], categories=['Yes','No'], ordered=True)
meta['kjonn'] = pd.Categorical(meta['kjonn'], categories=['Male','Female'], ordered=True)
meta['age_cat'] = pd.Categorical(meta['age_cat'], categories=['50-59','60-69','>=70'], ordered=True)
meta['final_result']=meta['final_result'].apply(lambda row: row.split('. ')[1])
meta['final_result']=pd.Categorical(meta['final_result'], categories=['Negative', 'Non neoplastic findings',
                                                                                 '>= 3 Non-advanced adenomas','Advanced adenoma','Advanced serrated',
                                                                                 'Cancer','Other lesions'], ordered=True)

samples=pyr.read_r('/'.join([wdir, 'participant_data/sample_meta.Rds']))[None]

#load metagenome count data

def load_data(prefix):
    
    def data_plasmids(cov,fold):
   
        cov=cov.rename(columns={'Unnamed: 0': 'sample_id'})
        cov=cov.set_index('sample_id')
        fold=fold.rename(columns={'Unnamed: 0': 'sample_id'})
        fold=fold.set_index('sample_id')
        cov=cov.applymap(lambda el: 0 if el < 75 else el)
        fold=fold.where(cov != 0,0)
        
        return fold


    if prefix=='MAGs':
        data=pd.read_csv('/'.join([wdir.replace('papers/Baseline_descriptive','datasets'),'metagenome/MAGs/filtered/counts/median_coverage_genomes.tsv']), sep='\t')


    elif prefix=='vOTUs':
        data=pd.read_csv('/'.join([wdir, 'datasets/Phages_coverage_min75.tsv']), sep='\t')
        data=data.set_index('ID')
        data=data.T
        data=data.reset_index(names='sample_id')
    
        
    elif prefix=='pOTUs':
        data=pd.read_csv('/'.join([wdir, 'datasets/Plasmids_coverage_min75.tsv']), sep='\t')
        
    return data

data=load_data(prefix)

#rename the samples
data['sample_id']=data['sample_id'].apply(lambda row: row.replace('-','_'))

#Keep only samples that are in the selected dataset
data=data.loc[data['sample_id'].isin(samples['sample_id'].to_list())]
data=data.set_index('sample_id')

#------------------------------------------------------------------------------------------------------------
#RELATIVE ABUNDANCE TABLE
#Remove taxa that are not present in any of the samples
taxa=data.sum(axis=0)
taxa=taxa[taxa!=0].index.to_list()
data=data[taxa]

#Scale the data

def scale_data(row):
    total=row.sum()
    row=row/total*100
    return row

relab=data.apply(scale_data, axis=1)

#Save relative abundance table
relab.to_csv('/'.join([wdir, 'datasets',prefix+'_relab.tsv']), sep='\t', index=True)

relab=pd.read_csv('/'.join([wdir, 'datasets',prefix+'_relab.tsv']), sep='\t').set_index('sample_id')

for p in ['MAGs','vOTUs','pOTUs']:
    relab=pd.read_csv('/'.join([wdir, 'datasets',p+'_relab.tsv']), sep='\t').set_index('sample_id')
    if 'allrelab' not in globals():
        allrelab=relab

    else:
        allrelab=allrelab.merge(relab,left_index=True,right_index=True)
        
#------------------------------------------------------------------------------------------------------------
#ALPHA DIVERSITY

def obs_sp(row):
    row=row[row!=0]
    obs_ind=len(row)
    return obs_ind

def inv_simpson(row):
    row=row/100
    s=row**2
    isim=1/s.sum()
    return isim

def shannon(row):
    row=row/100
    row=row[row!=0]
    shan_ind=-sum(row.apply(lambda x: x*math.log(x)))
    return shan_ind


obs=allrelab.apply(obs_sp,axis=1)
print('Observed species: '+ '{:.2f}'.format(obs.mean()) + '±'+'{:.2f}'.format(obs.std()))
shan_ind=allrelab.apply(shannon,axis=1)
print('Shannon index: '+ '{:.2f}'.format(shan_ind.mean()) + '±'+'{:.2f}'.format(shan_ind.std()))
inv_sim=allrelab.apply(inv_simpson,axis=1)
print('Inversed Simpsons: '+ '{:.2f}'.format(inv_sim.mean()) + '±'+'{:.2f}'.format(inv_sim.std()))


#Save data
prefix='combined'
alpha_div=pd.DataFrame({'ObsSp':obs, 'InvSimpson':inv_sim, 'Shannon': shan_ind})
alpha_div=alpha_div.reset_index().rename(columns={'index':'sample_id'})
alpha_div=alpha_div.merge(samples[['sample_id','deltaker_id','Prøvetype','FIT_value']], on='sample_id',how='left')
alpha_div=alpha_div.merge(meta, on='deltaker_id',how='left')

alpha_div.to_csv('/'.join([wdir, 'datasets/'+prefix+'_AlphaDiv.tsv']), sep='\t',index=False)

alpha_div=pd.read_csv('/'.join([wdir, 'datasets/'+prefix+'_AlphaDiv.tsv']), sep='\t')
alpha_div['final_result']=pd.Categorical(alpha_div['final_result'], categories=['Negative', 'Non neoplastic findings',
                                                                               '>= 3 Non-advanced adenomas','Advanced adenoma','Advanced serrated',
                                                                                'Cancer','Other lesions'], ordered=True)


## Statistics for alpha diversity
##Adjust for sequencing depth

def diff_adjusted(samples, y, col, adj):
    
    model = ols(f'{y} ~ C({col}) + {adj}', data=samples).fit()
    
    print(model.summary())

    results=pd.DataFrame({'Y':[y], 'Group':[col], 'Adj':[adj], 'Rsq':[model.rsquared], 'RsqAdj':[model.rsquared_adj],
                          'Fstat':[model.fvalue],'Prob_F':[model.f_pvalue],'PvalGroupVar':[model.pvalues[1]],'PvalAgjVar':[model.pvalues[2]],
                          'PvalIntercept':[model.pvalues[0]]})
    
    return results

divs=['ObsSp','Shannon','InvSimpson']
catg=['beforeBL','age_cat','kjonn','senter']
samples=pd.read_csv('/'.join([wdir,'datasets/data_by_samples_crispr.csv']),sep = '\t')

All_OLS=pd.DataFrame()
for p in ['MAGs','vOTUs','pOTUs','combined']:
    alpha_div=pd.read_csv('/'.join([wdir, 'datasets/'+p+'_AlphaDiv.tsv']), sep='\t')
    alpha_div=alpha_div.merge(samples[['deltaker_id','Total_Bases_QC_ATLAS']], on='deltaker_id', how='left')

    OLS_differ=pd.DataFrame()
    for y in divs:
        for gr in catg:
            olsres=diff_adjusted(alpha_div, y, gr, 'Total_Bases_QC_ATLAS')
            OLS_differ=pd.concat([OLS_differ, olsres])
    OLS_differ['Domain']=p
    All_OLS=pd.concat([All_OLS, OLS_differ])


All_OLS.to_csv('/'.join([wdir, 'results/OLS_alphadiv_adjusted.tsv']), sep='\t',index=False)

sb.catplot(All_OLS, col='Domain',x='RsqAdj',hue='Y', y='Group', kind='bar',sharex=False, palette=['#52a1b5','#9791c6','#acc0c3'])
