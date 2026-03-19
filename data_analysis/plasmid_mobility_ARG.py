#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Make a summary for plasmid mobility and ARG load 

"""

import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
import pyreadr as pyr
import scipy.stats as sst
from statsmodels.formula.api import ols
import statsmodels.api as sm
from itertools import combinations

wdir='PATH_TO_MANUS_FOLDER'

args=pd.read_csv('/'.join([wdir, 'datasets/plasmids_dereplicated_0.9_ARGs.txt']), sep='\t')
mob=pd.read_csv('/'.join([wdir, 'datasets/plasmids_dereplicated_0.9_MOBtyper.txt']), sep='\t')
antib=pd.read_csv('/'.join([wdir,'antibiotics', 'Antibiotics_per_id_categorical_timeadded.tsv']), sep='\t')

meta=pyr.read_r('/'.join([wdir, 'datasets/metadata/data_by_sample.rds']))[None]
relab=pd.read_csv('/'.join([wdir,'datasets/PTUs_relab.tsv']),sep='\t')

MobSum={'PTU':[], 'Mobility':[], 'NumARG_total':[],'NumARG_unique':[],
                     'Num_efflux_pump':[],'Num_drug_class':[]}

for _,p in mob.iterrows():
    arg=args.loc[args['Contig']==p.sample_id]
    MobSum['PTU'].append(p.sample_id)
    MobSum['Mobility'].append(p.predicted_mobility)

    if len(arg)>0:
        MobSum['NumARG_total'].append(len(arg))
        MobSum['NumARG_unique'].append(len(arg['Best_Hit_ARO'].unique().tolist()))
        MobSum['Num_efflux_pump'].append(len(arg.loc[arg['Resistance Mechanism'].str.contains('efflux')]))
        MobSum['Num_drug_class'].append(len(arg['Drug Class'].unique().tolist()))
    else:
        MobSum['NumARG_total'].append(0)
        MobSum['NumARG_unique'].append(0)
        MobSum['Num_efflux_pump'].append(0)
        MobSum['Num_drug_class'].append(0)

MobSum=pd.DataFrame(MobSum)

#Keep plasmids from the manuscript dataset
PTUlist=pd.read_csv('/'.join([wdir, 'datasets/PTU_rename_key.csv']), sep='\t')
MobSum=MobSum.loc[MobSum['PTU'].isin(PTUlist['PTU'].tolist())]

#Count percentage of plasmids in each mobility class

perc=(MobSum.loc[MobSum['NumARG_total']>0,'Mobility'].value_counts()/MobSum['Mobility'].value_counts()*100).reset_index()

fig,ax=plt.subplots(4,1)
sb.barplot(perc,x='count', y='Mobility', color='#ceba7d',ax=ax[0],order=['non-mobilizable','mobilizable','conjugative'])
sb.boxplot(MobSum.query('NumARG_total>0'),y='Mobility', x='NumARG_total', color='#ceba7d', ax=ax[1],
           order=['non-mobilizable','mobilizable','conjugative'])
sb.boxplot(MobSum.query('NumARG_total>0'),y='Mobility', x='Num_efflux_pump', color='#ceba7d', ax=ax[2],
           order=['non-mobilizable','mobilizable','conjugative'])
sb.boxplot(MobSum.query('NumARG_total>0'),y='Mobility', x='Num_drug_class', color='#ceba7d', ax=ax[3],
           order=['non-mobilizable','mobilizable','conjugative'])
ax[0].set(ylabel='',xlabel='Fraction of PTUs that contain ARGs')
ax[1].set(ylabel='',xlabel='Number of ARGs per PTU')
ax[2].set(ylabel='',xlabel='Number of efflux pumps per PTU')
ax[3].set(ylabel='',xlabel='Number of drug classes per PTU')

plt.savefig('/'.join([wdir, 'results/PTU_mobility_ARG_efflux_drugs.pdf']))
#------------------------------------------------------------------------------------
#Statistical tests
#Binomial test for fractions of plasmids; Null - non-mobilizable

perc=perc.rename(columns={'count':'Perc_plasmids'})
perc=perc.merge(MobSum.loc[MobSum['NumARG_total']>0,'Mobility'].value_counts().reset_index(),on='Mobility', how='left')
perc=perc.rename(columns={'count':'Num_plasmids_ARG'})
perc=perc.merge(MobSum['Mobility'].value_counts().reset_index(), on='Mobility', how='left')
perc=perc.rename(columns={'count':'Num_plasmids_total'})

nullval=perc.loc[perc['Mobility']=='non-mobilizable','Perc_plasmids'][0]/100

def bintest(row, nullval):
    res=sst.binomtest(row.Num_plasmids_ARG, row.Num_plasmids_total,nullval, alternative='two-sided')
    return res.pvalue

perc['BinomialP']=perc.apply(bintest,axis=1, nullval=nullval)
perc=perc.set_index('Mobility')

perc['Num_plasmids_nonARG']=perc['Num_plasmids_total']-perc['Num_plasmids_ARG']

chi2, p, dof, exp =sst.chi2_contingency(perc[['Mobility','Num_plasmids_ARG','Num_plasmids_nonARG']].set_index('Mobility')) 

print('\n')
print('Chi-squared test')
print(f'Chi2={chi2:.2f}, pval={p:.4f}')
print('\n')
print('Observed values:')
print(perc[['Mobility','Num_plasmids_ARG','Num_plasmids_nonARG']].set_index('Mobility'))
print('\n')
print('Expected values:')
print(pd.DataFrame(exp, index=perc['Mobility'].tolist(), columns=['Num_plasmids_ARG','Num_plasmids_nonARG']))

#Number of ARGs adjusted for sample sequencing depth

MobSum['sample_id']=MobSum['PTU'].apply(lambda row: row.split('_')[0])
MobSum['sample_id']=MobSum['sample_id'].apply(lambda row: row.replace('-','_'))
MobSum=MobSum.merge(meta[['sample_id','Total_Bases_QC_ATLAS']], on='sample_id', how='left')

####LINEAR REGRESSION, Number of resistance genes is a response (outcome), plasmid type - independent variable
def diff_adjusted(data, y, col, adj):
    
    model = ols(f'{y} ~ C({col}) + {adj}', data=data).fit()
    
    print(model.summary())

    results=pd.DataFrame({'Y':[y], 'Group':['_vs_'.join(data[col].unique().tolist())], 'Adj':[adj], 'Rsq':[model.rsquared], 'RsqAdj':[model.rsquared_adj],
                          'Fstat':[model.fvalue],'Prob_F':[model.f_pvalue],'PvalGroupVar':[model.pvalues[1]],'PvalAgjVar':[model.pvalues[2]],
                          'PvalIntercept':[model.pvalues[0]]})
    
    return results

OLS_differ=pd.DataFrame()

mobil=['non-mobilizable','mobilizable','conjugative']

detected=MobSum.query('NumARG_total>0')
for y in ['NumARG_total', 'NumARG_unique', 'Num_efflux_pump','Num_drug_class']:
    for gr in list(combinations(mobil,2)):
        olsres=diff_adjusted(detected.loc[detected['Mobility'].isin(gr)], y, 'Mobility', 'Total_Bases_QC_ATLAS')
        OLS_differ=pd.concat([OLS_differ, olsres])
        
OLS_differ.to_csv('/'.join([wdir, 'results/OLS_PTUs_ARGs_vs_mobility_onlyARGcontain.tsv']),sep='\t',index=False)

MobSum.to_csv('/'.join([wdir, 'results/PTUs_mobility_ARG_summary.tsv']),sep='\t',index=False)

#Get the list of all antibiotics that are detected ARGs against
all_antib=list(set(args['Antibiotic'].str.split('; ').explode().tolist()))
all_antib=[a for a in all_antib if type(a).__name__!='float']

#------------------------------------------------------------------------------------------------------
##Make a linear regression between days since last a/b use and a) number of ARG-containing plasmids; b) number of ARGs in them

def modify_days(val):
    if pd.isna(val):
        return 6205 #if no record of antibiotic use in the dataset, then turn to over 17 years
    else:
        return int(val.replace(' days',''))

antib['TimeFromLast']=antib['TimeFromLast'].apply(modify_days)

#Overview over which plasmids are detected in which samples
args=args.rename(columns={'Contig':'PTU'})
relab=relab.set_index('sample_id')

ARG_by_sample=pd.DataFrame({'NumPlas':[],'NumPlasARG':[],'NumARGtotal':[], 'NumARGunique':[]})
for ix, s in relab.iterrows():
    plas=s[s>0].index.tolist()
    plasarg=args.loc[args['PTU'].isin(plas)]
    ARG_by_sample.loc[ix,'NumPlas']=len(plas)
    ARG_by_sample.loc[ix,'NumPlasARG']=len(plasarg)
    ARG_by_sample.loc[ix,'NumARGtotal']=len(plasarg['Best_Hit_ARO'].tolist())
    ARG_by_sample.loc[ix,'NumARGunique']=len(plasarg['Best_Hit_ARO'].unique().tolist())

ARG_by_sample['PercPlasARG']=ARG_by_sample['NumPlasARG']/ARG_by_sample['NumPlas']*100
ARG_by_sample=ARG_by_sample.reset_index() 
ARG_by_sample=ARG_by_sample.rename(columns={'index':'sample_id'})

ARG_by_sample=ARG_by_sample.merge(meta[['sample_id','deltaker_id','Total_Bases_QC_ATLAS']], on='sample_id', how='left')   
ARG_by_sample=ARG_by_sample.merge(antib[['deltaker_id','TimeFromLast']], on='deltaker_id', how='left')   

#ARG_by_sample.to_csv('/'.join([wdir,'results/PlasmidARGs_by_sample.csv']), sep='\t', index=False)
ARG_by_sample=pd.read_csv('/'.join([wdir,'results/PlasmidARGs_by_sample.csv']), sep='\t')

ARG_by_sample=ARG_by_sample.query('TimeFromLast!=6205') #remove those that did not take antibiotics at all


models=dict()
model_stats={'Y':[],'AdjVar':[],'Rsq':[],'RsqAdj':[],'CoefDaysLastUse':[],'Fstat':[],'Prob-F':[],'PvalDaysLastUse':[],'PvalAdjVar':[],'PvalConst':[]}
adj_var='Total_Bases_QC_ATLAS'

for yvar in ['NumPlas','NumPlasARG','NumARGtotal','NumARGunique','PercPlasARG']:
    x=ARG_by_sample[['deltaker_id','TimeFromLast',adj_var]].set_index('deltaker_id')
    x=sm.add_constant(x)
    y=ARG_by_sample[['deltaker_id',yvar]].set_index('deltaker_id')

    #Fit ordinary least squares regression model
    model=sm.OLS(y,x).fit()
                  
    models[f'{y}']=model
                
    #Get model stats
    model_stats['Y'].append(yvar)
    model_stats['AdjVar'].append(adj_var)
    model_stats['Rsq'].append(models[f'{y}'].rsquared)
    model_stats['RsqAdj'].append(models[f'{y}'].rsquared_adj)
    model_stats['CoefDaysLastUse'].append(models[f'{y}'].params['TimeFromLast'])
    model_stats['Fstat'].append(models[f'{y}'].fvalue)
    model_stats['Prob-F'].append(models[f'{y}'].f_pvalue)
    model_stats['PvalDaysLastUse'].append(models[f'{y}'].pvalues['TimeFromLast'])
    model_stats['PvalAdjVar'].append(models[f'{y}'].pvalues[adj_var])
    model_stats['PvalConst'].append(models[f'{y}'].pvalues['const'])
                
model_stats=pd.DataFrame(model_stats)    
model_stats.to_csv('/'.join([wdir, f'results/Antibiotics_vs_NumPlasARGs_AdjustedLinearModel_{adj_var}.csv']),index=False)

fig=sb.regplot(ARG_by_sample, x='TimeFromLast',y='PercPlasARG', color='#CEBA7D')
fig.set(xlabel='Time from last antibiotics use, days', ylabel='Plasmids with ARG load, %')
plt.savefig('/'.join([wdir,'results/PercPlasARG_TimeFromLast_LinReg.pdf']))

##get the statistics on number of plasmids detected per sample (by SCAPP)

numplas=pd.read_csv('PATH_TO/scapp_snakemake_wf/data/dereplication/cluster_map.txt',sep='\t')
#Keep only PTUs that are in the manus
numplas=numplas.loc[numplas['PTU'].isin(relab.columns.tolist())]
numplas['sample']=numplas['genome'].apply(lambda row: row.split('_')[0])
persam=numplas['sample'].value_counts().reset_index()
persam['sample']=persam['sample'].apply(lambda row: row.replace('-','_'))
persam=persam.loc[persam['sample'].isin(relab.index.tolist())]
persam['count'].describe()
persam.to_csv('/'.join([wdir,'results/NumPlasmids_assembled_persample.csv']),index=False)

