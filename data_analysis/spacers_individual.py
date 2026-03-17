#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Find how many spacer clusters are shared or new in an individual 

"""

import pandas as pd
import pyreadr as pyr
import seaborn as sb
from scipy.stats import pearsonr
import scipy.stats as stats
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
import random


#Load all needed data
wdir="PATH_TO_MANUS_FOLDER"

spacers=pd.read_csv('/'.join([wdir,'datasets/spacers_manus_table.csv']),sep='\t')

#----------------------------------------------------------------------------
## Get number of spacers with/without targets detected in individual
#----------------------------------------------------------------------------

PTUtab=pd.read_csv('/'.join([wdir, 'datasets/PTUs_relab.tsv']), sep='\t').set_index('sample_id')
votutab=pd.read_csv('/'.join([wdir, 'datasets/vOTUs_relab.tsv']), sep='\t').set_index('sample_id')
PTUlen=pd.read_csv('/'.join([wdir, 'datasets/PTUs_lengths.csv']), sep=',')
votulen=pd.read_csv('/'.join([wdir, 'datasets/vOTUs_lengths.csv']), sep=',')

reltabs={'PTUs':PTUtab,'vOTUs':votutab}
lengths={'PTUs':PTUlen,'vOTUs':votulen}

def individuality_target(tartype, reltabs=reltabs, lengths=lengths, spacers=spacers):
    tarsp=spacers.dropna(subset=tartype)
    tarsp=tarsp[['Cluster','Spacers','Sample',tartype]]
    reltab=reltabs[tartype]
    tarmanus=reltab.columns.tolist()
    
    qlen=lengths[tartype]
    totallen=qlen['Length, bp'].sum()
    
    for ix,t in tarsp.iterrows():
        print(f'Checking targets for spacer nr {ix}')
        qtar=t[tartype].replace('[','')
        qtar=qtar.replace(']','')
        qtar=qtar.replace("'","")
        qtar=qtar.split(', ')
        qtar=[q for q in qtar if q in tarmanus]
        tarsp.loc[ix,'TotalNumTargets']=len(qtar)
        
        #Find which spacer targets are present in the individual
        relab=reltab.loc[t.Sample,qtar]
        #targets present in individual
        ptar=relab[relab>0].index.tolist()
        tarsp.loc[ix,'IndivNumTargets']=len(ptar)
        
        relab_ind=reltab.loc[t.Sample]
        inddb=relab_ind[relab_ind>0].index.tolist()
        indlen=qlen.loc[qlen[tartype[:-1]].isin(inddb),'Length, bp'].sum()
        tarsp.loc[ix,'IndivDatabase']=indlen

        #general targets - keep the target if it is NOT detected in individual, ie remove all vOTUs that are detected in individual
        gtar=[g for g in qtar if g not in inddb]
        tarsp.loc[ix,'PopulationNumTargets']=len(gtar)
        tarsp.loc[ix,'PopulationDatabase']=totallen-indlen
        
    return tarsp

virspacers=individuality_target('vOTUs')
virspacers=virspacers.query('TotalNumTargets>0') #drop those spacers where targets are only found in a total population, and not in manus samples
virspacers['IndividualityScore']=1-virspacers['PopulationNumTargets']/virspacers['TotalNumTargets']
virspacers['IndividualityNorm']=virspacers['IndivNumTargets']/virspacers['IndivDatabase']*10**6 #number of targets per 1MBp database
virspacers['PopulationNorm']=virspacers['PopulationNumTargets']/virspacers['PopulationDatabase']*10**6 #number of targets per 1MBp database
virspacers.to_csv('/'.join([wdir,'results/vOTU_spacers_targets_individuality.csv']),sep='\t',index=False)

fig=sb.scatterplot(virspacers,x='PopulationNorm',y='IndividualityNorm', size=.05, legend=False)

plasspacers=individuality_target('PTUs')
plasspacers=plasspacers.query('TotalNumTargets>0') #drop those spacers where targets are only found in a total population, and not in manus samples
plasspacers['IndividualityScore']=1-plasspacers['PopulationNumTargets']/plasspacers['TotalNumTargets']
plasspacers['IndividualityNorm']=plasspacers['IndivNumTargets']/plasspacers['IndivDatabase']*10**6 #number of targets per 1MBp database
plasspacers['PopulationNorm']=plasspacers['PopulationNumTargets']/plasspacers['PopulationDatabase']*10**6 #number of targets per 1MBp database
plasspacers.to_csv('/'.join([wdir,'results/PTU_spacers_targets_individuality.csv']),sep='\t',index=False)

fig=sb.scatterplot(plasspacers,x='PopulationNorm',y='IndividualityNorm', size=.05)

fig,ax=plt.subplots(1,2)
sb.scatterplot(virspacers,x='PopulationNorm',y='IndividualityNorm',size=0.2,marker='o',color='#6b519d',legend=False,ax=ax[0])
sb.scatterplot(plasspacers,x='PopulationNorm',y='IndividualityNorm',size=0.2,marker='o',color='#4d4f50',legend=False,ax=ax[1])
plt.savefig('/'.join([wdir,'results/Spacers_Individuality_Norm.pdf']))

virspacers['Type']='vOTU'
plasspacers['Type']='PTU'

alltargets=pd.concat([virspacers,plasspacers], ignore_index=True)
toplot=alltargets[['Cluster','Spacers','Sample','Type','IndividualityNorm','PopulationNorm']].melt(id_vars=['Cluster','Spacers','Sample','Type'],var_name='Level',value_name='NumTargets')

fig=sb.violinplot(toplot,y='Type',x='NumTargets',hue='Level',split=True, gap=0.1,inner='quart',cut=1, palette=['#5fc0bf','#acc0c3'],legend=False, inner_kws=dict(color='#be2928'))
fig.set(xlabel='Number of targets per Mb viral genomes',ylabel='')
plt.savefig('/'.join([wdir,'results/Spacers_Individuality_NumberOfHits_perMB.pdf']))

h_st,pval=kruskal_group(toplot.loc[toplot['Type']=='vOTU'],'Level','NumTargets') 
print(f'Kruskal-Wallis for vOTUs: H={h_st:.2f},pval={pval:.3f}')

h_st,pval=kruskal_group(toplot.loc[toplot['Type']=='PTU'],'Level','NumTargets') 
print(f'Kruskal-Wallis for PTUs: H={h_st:.2f},pval={pval:.3f}')

## Make a contingency table for the individuality of spacers

def or_chi2(tartype, alltargets=alltargets, spacers=spacers):
    query=spacers[['Cluster','Spacers']]
    query=query.merge(alltargets.loc[alltargets['Type']==tartype,['Spacers','IndividualityNorm','PopulationNorm']], on='Spacers',how='left')
    query['IndTar']=query['IndividualityNorm'].apply(lambda row: 'Yes' if row>0 else 'No')
    query['PopTar']=query['PopulationNorm'].apply(lambda row: 'Yes' if row>0 else 'No')
    crtab=pd.crosstab(query['PopTar'],query['IndTar'])
    
    OR=crtab.loc['Yes','Yes']*crtab.loc['No','No']/crtab.loc['Yes','No']*crtab.loc['No','Yes']
    chi2, p_chi2, dof, expected = stats.chi2_contingency(crtab, correction=False)
    
    expected=pd.DataFrame(expected,index=crtab.index.tolist(),columns=crtab.columns.tolist())
    return OR,chi2,p_chi2,crtab,expected

tartype='PTU'
OR,chi2,p_chi2,crtab,expected=or_chi2(tartype)

print('----------------------------------')
print(f'Contingency table for {tartype}')
print('----------------------------------\n')
print(f'Odds ratio: {OR.astype(int)}')
print(f'Chi-square test:chi2={chi2:.2f}; pvalue={p_chi2:.2f}')
print('\nObserved values:')
print(crtab)
print('\nExpected values:')
print(expected)
print('----------------------------------')

##Check the individuality of targets vs other random individuals

seed=42
random.seed=seed

def individ_vs_individ(tartype, ranind, reltabs=reltabs, spacers=spacers, samples=samples):
    tarsp=spacers.dropna(subset=tartype)
    tarsp=tarsp[['Cluster','Spacers','Sample',tartype]]
    reltab=reltabs[tartype]
    tarmanus=reltab.columns.tolist()

    for ix,t in tarsp.iterrows():
        print(f'Checking targets for spacer nr {ix}')
        qtar=t[tartype].replace('[','')
        qtar=qtar.replace(']','')
        qtar=qtar.replace(' ','')
        qtar=qtar.split(',')
        qtar=[q for q in qtar if q in tarmanus] # spacer targets
        tarsp.loc[ix,'TotalNumTargets']=len(qtar)

        # Find which spacer targets are present in the individual
        relab=reltab.loc[t.Sample,qtar]
        # spacer targets present in individual
        ptar=relab[relab>0].index.tolist()
        tarsp.loc[ix,'IndivNumTargets']=len(ptar)
        tarsp.loc[ix,'IndivNumTargets_perGBseq']=len(ptar)/samples.loc[samples['Sample']==t.Sample,'Total_Bases_QC_ATLAS'].item()*10**9

        # take {ranind} random individuals from the population
        indiv=reltab.index.tolist()
        indiv=[i for i in indiv if i!=t.Sample]
        indlist=random.sample(indiv,ranind)
        itarlens=[]
        itarnonorm=[]
        for i in indlist:
            irelab=reltab.loc[i,qtar]
            itar=irelab[irelab>0].index.tolist()
            itarnonorm.append(len(itar))
            itarlens.append(len(itar)/samples.loc[samples['Sample']==i,'Total_Bases_QC_ATLAS'].item()*10**9)

        tarsp.loc[ix,'RandomIndNumTargets_perGBseq']=itarlens
        tarsp.loc[ix,'RandomIndNumTargets_Median']=np.median(itarnonorm)
        tarsp.loc[ix,'RandomIndNumTargets_Median_perGBseq']=np.median(itarlens)

    return tarsp


ranind=10
virindsp=individ_vs_individ('vOTUs',ranind)
plasindsp=individ_vs_individ('PTUs',ranind)

virindsp['Delta'] = virindsp['IndivNumTargets_perGBseq'] - virindsp['RandomIndNumTargets_Median_perGBseq']
plasindsp['Delta'] = plasindsp['IndivNumTargets_perGBseq'] - plasindsp['RandomIndNumTargets_Median_perGBseq']

virindsp.to_csv('/'.join([wdir, f'results/vOTU_spacers_targets_individ_vs_individ_{ranind}_random_ind.csv']), sep='\t', index=False)
plasindsp.to_csv('/'.join([wdir, f'results/PTU_spacers_targets_individ_vs_individ_{ranind}_random_ind.csv']), sep='\t', index=False)

virindsp['Sum'] = virindsp['IndivNumTargets'] + virindsp['RandomIndNumTargets_Median']  # exclude those that don't have any targets
plasindsp['Sum'] = plasindsp['IndivNumTargets'] + plasindsp['RandomIndNumTargets_Median']  # exclude those that don't have any targets


def make_plot(toplot):
    toplot = toplot.sort_values(by='Delta', ascending=False).reset_index()
    plt.scatter(toplot.index.tolist(), toplot['Delta'], color='#5fc0bf', s=5, zorder=1)
    plt.vlines(toplot.index.tolist(), ymin=0, ymax=toplot['Delta'], color='#5fc0bf', linewidth=1)
    plt.axhline(0, color='black', linewidth=1, linestyle='--')

    # plt.yscale('log')


make_plot(virindsp.query('Sum>0'))
plt.ylim([-2, 8])
plt.savefig('/'.join([wdir, f'results/vOTU_spacers_targets_individ_vs_individ_{ranind}_random_ind.pdf']))

make_plot(plasindsp.query('Sum>0'))
plt.savefig('/'.join([wdir, f'results/PTU_spacers_targets_individ_vs_individ_{ranind}_random_ind.pdf']))


mges = pd.concat([plasindsp, virindsp])
make_plot(mges.query('Sum>0'))
plt.savefig('/'.join([wdir, f'results/MGE_spacers_targets_individ_vs_individ_{ranind}_random_ind.pdf']))

# Find how many clusters have target detected either in individual or random individual
len(mges.query('Sum>0')['Cluster'].unique().tolist())


# Find how many have targets detected in individuals
def find_stats(indsp, tartype):
    print(f"Total number of spacer clusters with {tartype} targets: {len(indsp['Cluster'].unique().tolist())}")
    print(f"Number of spacer clusters with {tartype} target in individual: {len(indsp.query('IndivNumTargets>0')['Cluster'].unique().tolist())}")


find_stats(virindsp, 'vOTUs')
find_stats(plasindsp, 'PTUs')
find_stats(pd.concat([plasindsp, virindsp]), 'MGE')


#Check if there is significant difference between random and individual
def kruskal_group(data,cat_col,y):
    categ=data[cat_col].unique().tolist()
    for_stats={}
    for cat in categ:
        for_stats[cat]=data.loc[data[cat_col]==cat,y].tolist()

    h_st,pval=stats.kruskal(*list(for_stats.values())) # *indicates to treat each element in list as a group, the list will contain as many elements as there are dict entries, each entry is a separate list within a list

    return h_st,pval


virindsp = pd.melt(
    virindsp[['Cluster','Sample','vOTUs','IndivNumTargets_perGBseq','RandomIndNumTargets_Median_perGBseq']],
    id_vars=['Cluster','Sample','vOTUs'],
    var_name='Group',
    value_name='NumTargets'
)

h_st, pval = kruskal_group(virindsp, 'Group', 'NumTargets')
print(f'H={h_st:.2f}; p={pval:.2f}')


plasindsp = pd.melt(
    plasindsp[['Cluster','Sample','PTUs','IndivNumTargets_perGBseq','RandomIndNumTargets_Median_perGBseq']],
    id_vars=['Cluster','Sample','PTUs'],
    var_name='Group',
    value_name='NumTargets'
)

h_st, pval = kruskal_group(plasindsp, 'Group', 'NumTargets')
print(f'H={h_st:.2f}; p={pval:.2f}')