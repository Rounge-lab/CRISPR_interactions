#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Diversity comparisons between MAGs and Phages

"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
from itertools import combinations as comb
import pyreadr as pyr
import statsmodels.api as sm
import os

wdir='PATH_TO_MANUS_FOLDER' #set working directory

prefix='MAGs_vOTUs_PTUs' #prefix to be used for the output files
os.makedirs('/'.join([wdir, prefix]), exist_ok=True)
#------------------------------------------------------------------------------------------------------------
#Load alpha diversity
mag_alpha=pd.read_csv('/'.join([wdir, 'datasets/MAGs_AlphaDiv.tsv']), sep='\t')
votu_alpha=pd.read_csv('/'.join([wdir, 'datasets/vOTUs_AlphaDiv.tsv']), sep='\t')
PTU_alpha=pd.read_csv('/'.join([wdir, 'datasets/PTUs_AlphaDiv.tsv']), sep='\t')

#Load metadata

meta=pd.read_csv('/'.join([wdir,'participant_data/screening_data.tsv']), sep='\t')
meta['beforeBL'] = pd.Categorical(meta['beforeBL'], categories=['Yes','No'], ordered=True)
meta['kjonn'] = pd.Categorical(meta['kjonn'], categories=['Male','Female'], ordered=True)
meta['age_cat'] = pd.Categorical(meta['age_cat'], categories=['50-59','60-69','>=70'], ordered=True)
meta['final_result']=meta['final_result'].apply(lambda row: row.split('. ')[1])

#------------------------------------------------------------------------------------------------------------
#Check Linear correlation between mags, phages and plasmid alpha diversity
#Linear model adjusted for the sequencing depth

samples=pyr.read_r('/'.join([wdir, 'participant_data/sample_meta.Rds']))[None]
mag_alpha=mag_alpha.merge(samples[['sample_id','Total_Bases_QC_ATLAS','N50']], on='sample_id',how='left')
votu_alpha=votu_alpha.merge(samples[['sample_id','Total_Bases_QC_ATLAS','N50']], on='sample_id',how='left')
PTU_alpha=PTU_alpha.merge(samples[['sample_id','Total_Bases_QC_ATLAS','N50']], on='sample_id',how='left')

alpha={'MAG':mag_alpha, 'vOTU': votu_alpha, 'PTU':PTU_alpha}
sids=samples['sample_id'].tolist()

def adjust_model(key1,key2, div, adj_var, alpha=alpha, sids=sids):
    a1=alpha[key1]
    a2=alpha[key2]
    
    #make sure the order is the same
    a1=a1.set_index('sample_id').reindex(sids)
    a2=a2.set_index('sample_id').reindex(sids)
    
    x=a1[[div,adj_var]]
    x=sm.add_constant(x)
    y=a2[div]
    
    #Fit ordinary least squares regression model
    model=sm.OLS(y,x).fit()
    return model

def plot_models(div,adj_var,alpha=alpha):
    models=dict()
    model_stats={'Comparison':[],'DivInd':[],'AdjVar':[],'Rsq':[],'RsqAdj':[],'Fstat':[],'Prob-F':[],'PvalDaysLastUse':[],'PvalAdjVar':[],'PvalConst':[]}

    pairs=list(comb(list(alpha.keys()),2))

    fig,axs=plt.subplots(3,1,figsize=(5,10))
    
    ax=0 #set the counter for the axes

    for pair in pairs:
        model=adjust_model(pair[0],pair[1],div,adj_var)
        models[f'{pair[0]}_{pair[1]}']=model
        
        #Get model stats
        model_stats['Comparison'].append(f'{pair[0]}_{pair[1]}')
        model_stats['DivInd'].append(div)
        model_stats['AdjVar'].append(adj_var)
        model_stats['Rsq'].append(models[f'{pair[0]}_{pair[1]}'].rsquared)
        model_stats['RsqAdj'].append(models[f'{pair[0]}_{pair[1]}'].rsquared_adj)
        model_stats['Fstat'].append(models[f'{pair[0]}_{pair[1]}'].fvalue)
        model_stats['Prob-F'].append(models[f'{pair[0]}_{pair[1]}'].f_pvalue)
        model_stats['PvalDaysLastUse'].append(models[f'{pair[0]}_{pair[1]}'].pvalues[div])
        model_stats['PvalAdjVar'].append(models[f'{pair[0]}_{pair[1]}'].pvalues[adj_var])
        model_stats['PvalConst'].append(models[f'{pair[0]}_{pair[1]}'].pvalues['const'])
        
        #Plot the model
        x=alpha[pair[0]].set_index('sample_id').reindex(sids)
        x=x[[div,adj_var]]
        x=sm.add_constant(x)
        preds = models[f'{pair[0]}_{pair[1]}'].get_prediction(x)
        pred_summary = preds.summary_frame(alpha=0.05)  # 95% confidence intervals
        pred_summary=pred_summary.reindex(sids)
        x=x.join(pred_summary, how='inner')

        # Plotting Model
        # Scatter plot of the original data points
        y=alpha[pair[1]].set_index('sample_id').reindex(sids)
        toscatter=pd.DataFrame({pair[0]:x[div].tolist(),pair[1]:y[div].tolist()})
        sb.scatterplot(data=toscatter, x=pair[0], y=pair[1], color='#5fc0bf', ax=axs[ax],alpha=0.8,size=2,legend=False)
        sb.regplot(data=x, x=div, y='mean', color='#5fc0bf', ax=axs[ax], scatter=False, ci=95)
        axs[ax].set(xlabel=pair[0],ylabel=pair[1])
        ax+=1

    model_stats=pd.DataFrame(model_stats)
    return model_stats

for div in ['InvSimpson','ObsSp','Shannon']:
    for adj_var in ['Total_Bases_QC_ATLAS','N50']:
        model_stats=plot_models(div,adj_var)

        model_stats.to_csv('/'.join([wdir, prefix, f'{div}_AdjustedLinearModel_{adj_var}.csv']),index=False)
        plt.tight_layout()
        plt.savefig('/'.join([wdir, prefix, f'{div}_AdjustedLinearModel_{adj_var}.pdf']))

#Make a graph for alpha diversity vs antibiotics use

antib=pd.read_csv('/'.join([wdir,'antibiotics', 'Antibiotics_per_id_categorical_timeadded.tsv']), sep='\t')

mag_alpha['Domain']='MAGs'
votu_alpha['Domain']='vOTUs'
PTU_alpha['Domain']='PTUs'
alpha=pd.concat([mag_alpha,votu_alpha,PTU_alpha])

alpha=alpha.merge(antib[['deltaker_id','TimeFromLast']], on='deltaker_id')
alpha.to_csv('/'.join([wdir, 'datasets/MAGsvOTUsPTUs_AlphaDiv.tsv']), sep='\t', index=False)

def modify_days(val):
    if pd.isna(val):
        return 6205 #if no record of antibiotic use in the dataset, then turn to over 17 years
    else:
        return int(val.replace(' days',''))

alpha['DaysLastUse']=alpha['TimeFromLast'].apply(modify_days)
alpha=alpha.query('DaysLastUse<6205') #exclude unknown values

alpha=alpha.merge(samples[['sample_id','Total_Bases_QC_ATLAS']], on='sample_id', how='left')

#Fit linear regression models (days since last a/b use to alpha diversity), adjusted for seq depth

colors=['#4385BF','#5B7571','#CEBA7D']

models=dict()
model_stats={'Domain':[],'DivInd':[],'AdjVar':[],'Rsq':[],'RsqAdj':[],'CoefDaysLastUse':[],'Fstat':[],'Prob-F':[],'PvalDivInd':[],'PvalAdjVar':[],'PvalConst':[]}

col=0 #initialize color count

adj_var='Total_Bases_QC_ATLAS'
for p in ['MAGs','vOTUs', 'PTUs']:
    fig,axs=plt.subplots(3,1,figsize=(5,10))
    ax=0

    for div in ['InvSimpson','ObsSp','Shannon']:
        x=alpha.loc[alpha['Domain']==p, ['deltaker_id','DaysLastUse',adj_var]]
        x=x.sort_values(by='deltaker_id')
        ids=x['deltaker_id'].tolist()
        x=x.set_index('deltaker_id').reindex(ids)
        x=sm.add_constant(x)
        y=alpha.loc[alpha['Domain']==p,['deltaker_id',div]]
        y=y.set_index('deltaker_id').reindex(ids)

        #Fit ordinary least squares regression model
        model=sm.OLS(y,x).fit()
                    
        models[f'{p}_{div}']=model
                
        #Get model stats
        model_stats['Domain'].append(p)
        model_stats['DivInd'].append(div)
        model_stats['AdjVar'].append(adj_var)
        model_stats['Rsq'].append(models[f'{p}_{div}'].rsquared)
        model_stats['RsqAdj'].append(models[f'{p}_{div}'].rsquared_adj)
        model_stats['CoefDaysLastUse'].append(models[f'{p}_{div}'].params['DaysLastUse'])
        model_stats['Fstat'].append(models[f'{p}_{div}'].fvalue)
        model_stats['Prob-F'].append(models[f'{p}_{div}'].f_pvalue)
        model_stats['PvalDivInd'].append(models[f'{p}_{div}'].pvalues['DaysLastUse'])
        model_stats['PvalAdjVar'].append(models[f'{p}_{div}'].pvalues[adj_var])
        model_stats['PvalConst'].append(models[f'{p}_{div}'].pvalues['const'])
                
        # Plotting Model
        # Scatter plot of the original data points
        toscatter=alpha.loc[alpha['Domain']==p,['DaysLastUse',div]]
        sb.scatterplot(data=toscatter, x='DaysLastUse', y=div, color=colors[col], ax=axs[ax],alpha=0.8,size=2,legend=False)
        sb.regplot(data=toscatter, x='DaysLastUse', y=div, color=colors[col], ax=axs[ax], scatter=False, ci=95)
        axs[ax].set(xlabel='Days since last a/b use',ylabel=f'{div} index')
        ax+=1
    
    plt.tight_layout()
    #plt.savefig('/'.join([wdir,prefix,f'Antibiotics_LinReg_adjSeqDepth_{p}.pdf']))
    col+=1
model_stats=pd.DataFrame(model_stats)    
model_stats.to_csv('/'.join([wdir, prefix, f'Antibiotics_vs_Alpha_AdjustedLinearModel_{adj_var}.csv']),index=False)

##Plot diversity indices at different groups

alpha=alpha[['Domain','sample_id','InvSimpson','age_cat','kjonn','senter','beforeBL']]
alpha=pd.melt(alpha, id_vars=['Domain','sample_id','InvSimpson'],var_name='Cat',value_name='Value')

y_tick_values = {'age_cat': ['50-59','60-69','>=70'], 'kjonn': ['Male','Female'],'senter': ['Bærum','Moss'],
    'beforeBL': ['Yes','No']}

for c in y_tick_values:
    fig=sb.catplot(alpha.loc[alpha['Cat']==c],col='Domain',hue='Domain',y='Value',x='InvSimpson',sharey=False,
                   sharex=False,kind='violin', col_order=['MAGs','vOTUs','PTUs'], palette=colors, order=y_tick_values[c])
    fig.set(ylabel='',xlabel='Inverse Simpson',title='')
    plt.savefig('/'.join([wdir, f'results/InvSimpson_Violin_{c}.pdf']))

#Plot only antibiotics and region

fig,ax=plt.subplots(1,3)

dom=['MAGs','PTUs','vOTUs']
toplot=alpha.loc[alpha['Cat'].isin(['beforeBL','senter'])]
for i in range(len(dom)):
    sb.violinplot(toplot.loc[toplot['Domain']==dom[i]], hue='Value',y='Cat',x='InvSimpson',
               ax=ax[i], palette=[colors[i]]*2,split=True, gap=.1, inner='quart', 
               hue_order=['Yes','No','Moss','Bærum'],legend=False,order=['beforeBL','senter'])
    ax[i].set(ylabel='',xlabel='Inverse Simpson')
    if i==0:
        ax[i].set_yticklabels(['a/b use','Region'])
    else:
        ax[i].set_yticklabels('')

plt.savefig('/'.join([wdir, 'results/InvSimpson_Violin_split_antib_senter.pdf']))

#-------------------------------------------------------------------------------------------------
#Compare BC indexes across domains

#Load beta diversity
mag_bray=pd.read_csv('/'.join([wdir, 'datasets/MAGs_BrayCurtis.tsv']), sep='\t')
votu_bray=pd.read_csv('/'.join([wdir, 'datasets/vOTUs_BrayCurtis.tsv']), sep='\t')
PTU_bray=pd.read_csv('/'.join([wdir, 'datasets/PTUs_BrayCurtis.tsv']), sep='\t')


def melt_bray(bc,domain):
    bc_m=bc.melt(id_vars='Unnamed: 0',var_name='Sample2', value_name='BC')
    bc_m['Domain']=domain
    bc_m=bc_m.rename(columns={'Unnamed: 0':'Sample1'})
    
    return bc_m

mag_bray=melt_bray(mag_bray,'MAGs')
votu_bray=melt_bray(votu_bray,'vOTUs')
PTU_bray=melt_bray(PTU_bray,'PTUs')

all_bray=pd.concat([mag_bray,votu_bray,PTU_bray])
all_bray=all_bray.query('BC>0') #remove compare to self

all_bray['Comb']=all_bray.apply(lambda row: [row.Sample1, row.Sample2, row.Domain], axis=1)
all_bray['Comb']=all_bray['Comb'].apply(lambda row: '_'.join(list(sorted(row))))

all_bray=all_bray.drop_duplicates(subset='Comb', keep='first')

sb.boxplot(all_bray,x='BC',y='Domain',hue='Domain',hue_order=['MAGs', 'vOTUs', 'PTUs'], palette=colors)
plt.tight_layout()
plt.savefig('/'.join([wdir, prefix, 'BrayCurtis_by_domains.pdf']))


def kruskal_group(data,cat_col,y):
    query=data.copy()
    query=query.dropna(subset=[cat_col,y])
    categ=query[cat_col].unique().tolist()
    for_stats={}
    for cat in categ:
        for_stats[cat]=query.loc[query[cat_col]==cat,y].tolist()

    h_st,pval=stats.kruskal(*list(for_stats.values())) # *indicates to treat each element in list as a group, the list will contain as many elements as there are dict entries, each entry is a separate list within a list

    return h_st,pval

def kruskal_pairwise(domains,all_bray=all_bray):
    Kruskaldb=pd.DataFrame()
    dompairs=list(comb(domains,2))
    for d in dompairs:
        test=all_bray.loc[all_bray['Domain'].isin(d)]
        h_st,pval=kruskal_group(test,'Domain','BC')
        k=pd.DataFrame.from_dict({'Domains':[f'{d[0]}_vs_{d[1]}'],'Hstat': [h_st], 'Pval': [pval]})
        Kruskaldb=pd.concat([Kruskaldb,k])
    _,padj,_,_=multipletests(Kruskaldb['Pval'], method='fdr_bh')
    Kruskaldb['FDRp']=padj 
    
    return Kruskaldb

KruskalBeta=kruskal_pairwise(['MAGs','PTUs','vOTUs'])
KruskalBeta.to_csv('/'.join([wdir, prefix, 'Kruskal_BrayCurtis.tsv']),index=False)

#Compare those that had a/b within 4 months prior to baseline vs those that didn't have 

samples=samples.merge(meta[['deltaker_id','beforeBL']], on='deltaker_id', how='left')
all_bray=all_bray.merge(samples[['sample_id', 'beforeBL']], left_on='Sample1', right_on='sample_id', how='left')
all_bray=all_bray.merge(samples[['sample_id', 'beforeBL']], left_on='Sample2', right_on='sample_id', how='left')

all_bray['beforeBL_comb']=all_bray.apply(lambda row: '_'.join(list(set([samples.loc[samples['sample_id']==row.Sample1,'beforeBL'].values[0], 
                                                                        samples.loc[samples['sample_id']==row.Sample2,'beforeBL'].values[0]]))), axis=1)

sb.boxplot(all_bray,x='BC',y='Domain',hue='beforeBL_comb',hue_order=['Yes', 'No', 'No_Yes'], palette=['#ff5757','#ace1a5','#ffbd59'])
plt.tight_layout()
plt.savefig('/'.join([wdir, prefix, 'BrayCurtis_by_domains_ABuse_4moprior.pdf']))

Kruskal_ab=pd.DataFrame()
for d in ['MAGs', 'vOTUs', 'PTUs']:
    qbray=all_bray.loc[all_bray['Domain']==d]
    h_st,pval=kruskal_group(qbray,'beforeBL_comb','BC')
    k=pd.DataFrame.from_dict({'Domain':[f'{d}'],'Hstat': [h_st], 'Pval': [pval]})
    Kruskal_ab=pd.concat([Kruskal_ab, k])

_,padj,_,_=multipletests(Kruskal_ab['Pval'], method='fdr_bh')
Kruskal_ab['FDRp']=padj 

Kruskal_ab.to_csv('/'.join([wdir, prefix, 'Kruskal_BrayCurtis_domains_abuse4mo.tsv']),index=False)

Kruskal_ab=pd.DataFrame()
for d in ['MAGs', 'vOTUs', 'PTUs']:
    qbray=all_bray.loc[all_bray['Domain']==d]
    pairs=list(comb(['Yes','No','No_Yes'],2))
    for p in pairs:
        h_st,pval=kruskal_group(qbray.loc[qbray['beforeBL_comb'].isin(p)],'beforeBL_comb','BC')
        k=pd.DataFrame.from_dict({'Domain':[f'{d}'],'Combination':[f'{p[0]} vs {p[1]}'],'Hstat': [h_st], 'Pval': [pval]})
        Kruskal_ab=pd.concat([Kruskal_ab, k])

_,padj,_,_=multipletests(Kruskal_ab['Pval'], method='fdr_bh')
Kruskal_ab['FDRp']=padj 
