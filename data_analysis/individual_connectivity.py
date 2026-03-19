#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Individual connectivity based on the MAG-target pairs

"""

import pandas as pd
import numpy as np
from scipy.spatial.distance import jaccard
import pyreadr as pyr
import matplotlib.pyplot as plt
from itertools import combinations as comb

from skbio import DistanceMatrix
from skbio.stats.distance import permanova
from statsmodels.formula.api import ols
import seaborn as sb

from statannotations.Annotator import Annotator as annot


wdir='PATH_TO_MANUS_FOLDER'

meta=pd.read_csv('/'.join([wdir,'participant_data/screening_data.tsv']), sep='\t')
meta['beforeBL'] = pd.Categorical(meta['beforeBL'], categories=['Yes','No'], ordered=True)
meta['sex'] = pd.Categorical(meta['kjonn'], categories=['Male','Female'], ordered=True)
meta['age_cat'] = pd.Categorical(meta['age_cat'], categories=['50-59','60-69','>=70'], ordered=True)

samples=pyr.read_r('/'.join([wdir, 'participant_data/sample_meta.Rds']))[None]
samples=samples.merge(meta[['deltaker_id','beforeBL','age_cat','sex','senter']],on='deltaker_id',how='left')

spacers=pd.read_csv('/'.join([wdir,'datasets/spacers_manus_table_new.csv']), sep='\t')

edges=spacers[['Cluster','Sample','MAG','vOTUs','PTUs']]
edges=pd.melt(edges, id_vars=['Cluster','MAG','Sample'], value_vars=['vOTUs','PTUs'], var_name='domain', value_name='Taxon')
edges=edges.dropna(subset='Taxon')
edges=edges.dropna(subset='MAG')

edges['Taxon']=edges['Taxon'].apply(lambda row: row.replace("'",''))
edges['Taxon']=edges['Taxon'].apply(lambda row: row.replace("[",''))
edges['Taxon']=edges['Taxon'].apply(lambda row: row.replace("]",''))
edges['Taxon']=edges['Taxon'].str.split(', ')
edges = edges.explode('Taxon')

#Keep those that are in the manus
mags=pd.read_csv('/'.join([wdir, 'datasets/MAGs_relab.tsv']),sep='\t').set_index('sample_id')
votus=pd.read_csv('/'.join([wdir, 'datasets/vOTUs_relab.tsv']),sep='\t').set_index('sample_id')
PTUs=pd.read_csv('/'.join([wdir, 'datasets/PTUs_relab.tsv']),sep='\t').set_index('sample_id')

mges=votus.columns.tolist()
mges.extend(PTUs.columns.tolist())

edges=edges.loc[edges['Taxon'].isin(mges)]
edges=edges.loc[edges['MAG'].isin(mags.columns.tolist())]

#Remove MGEs that have crisprs
plcrispr=pd.read_csv('/'.join([wdir,'results/CCTyper_plasmids_rerun/crisprs_all.tab']),sep='\t')
vircrispr=pd.read_csv('/'.join([wdir,'results/CCTyper_viruses/crisprs_all.tab']),sep='\t')
edges=edges.loc[~edges['Taxon'].isin(plcrispr['Contig'].unique().tolist())]
edges=edges.loc[~edges['Taxon'].isin(vircrispr['Contig'].unique().tolist())]

edges['pair']=edges['MAG']+'_'+edges['Taxon']
numpairs=edges['pair'].value_counts()

edges.to_csv('/'.join([wdir,'results/Edges_with_samples_network_crisprfree.csv']),sep='\t', index=False)

prestab=pd.crosstab(edges['Sample'],edges['pair'])
prestab[prestab>0]=1
prestab.to_csv('/'.join([wdir,'datasets/Individual_interactions_PresAbs_CrisprF.csv']),sep='\t',index=True)

#Get number of individuals where each interaction is detected
numind=prestab.sum(axis=0).reset_index()
numind.columns=['pair','NumInd']

numpairs=numpairs.reset_index()
numpairs=numpairs.rename(columns={'count':'NumInteractions'})
numpairs=numpairs.merge(numind,on='pair', how='left')

numpairs['delta']=numpairs['NumInteractions']-numpairs['NumInd']
numpairs.sort_values(by='delta',ascending=True)

numpairs[['pair','NumInteractions','NumInd']].to_csv('/'.join([wdir,'results/Interactions_counts_by_spacer_and_indiv.csv']), index=False)

#Number of interactions per sample
numint=prestab.sum(axis=1)
numint=numint.rename(columns={'Sample':'sample_id', 0:'NumInt'})
numint.to_csv('/'.join([wdir,'datasets/Individual_interactions_NumPerInd_CrisprF.csv']),sep='\t',index=True)

#Number of individuals per interaction
prestab=pd.read_csv('/'.join([wdir,'datasets/Individual_interactions_PresAbs_CrisprF.csv']),sep='\t').set_index('Sample')
numind=prestab.sum(axis=0).reset_index()
numind.columns=['pair','NumIndividuals']
numind=numind.sort_values(by='NumIndividuals', ascending=False)
numind.to_csv('/'.join([wdir,'results/Individuals_per_interaction_crisprfree.csv']), index=False)
fig=sb.histplot(numind, bins=15, color='#ebf3f1', legend=False)
fig.set(xlabel='Number of individuals with same host-target pairs', ylabel='Number of host-target pairs')
fig.set(yscale='log')

#---------------------------------------------------------------------------
#check if number of interactions per individual is significantly different between groups
#---------------------------------------------------------------------------

numint=pd.read_csv('/'.join([wdir,'datasets/Individual_interactions_NumPerInd_CrisprF.csv']),sep='\t')
numint=numint.merge(samples[['sample_id','age_cat','beforeBL','sex','center','Total_Bases_QC_ATLAS']],
                    on='sample_id', how='left')

#add info on DMMs to the samples
dmm=pd.read_csv('/'.join([wdir.replace('Baseline_descriptive','CRCbiome_main'), 'analyses/data/dmm/crcbiome_dmm.tsv']), sep='\t')
numint=numint.merge(dmm, on='sample_id',how='left')
numint=numint.rename(columns={'gr':'dmm'})

def diff_adjusted(data, y, col, adj):
    
    model = ols(f'{y} ~ C({col}) + {adj}', data=data).fit()
    
    print(model.summary())

    results=pd.DataFrame({'Y':[y], 'Group':[col], 'Adj':[adj], 'Rsq':[model.rsquared], 'RsqAdj':[model.rsquared_adj],
                          'Fstat':[model.fvalue],'Prob_F':[model.f_pvalue],'PvalGroupVar':[model.pvalues[1]],'PvalAgjVar':[model.pvalues[2]],
                          'PvalIntercept':[model.pvalues[0]]})
    
    return results

OLS_differ=pd.DataFrame()

for gr in ['beforeBL','age_cat','sex','center','dmm']:
    olsres=diff_adjusted(numint, 'NumInt', gr, 'Total_Bases_QC_ATLAS')
    OLS_differ=pd.concat([OLS_differ, olsres])
        
OLS_differ.to_csv('/'.join([wdir, 'results/OLS_NumberInteractions_SeqDepthAdj_CrisprF.tsv']),sep='\t',index=False)

#Plot variation between number of interactions at different age categories
ax=sb.swarmplot(numint,y='NumInt',x='age_cat',color='#5fc0bf',zorder=1)
sb.boxplot(numint,y='NumInt',x='age_cat',fill=True, width=0.1,fliersize=0,showcaps=True, 
           boxprops={'facecolor':'none', 'edgecolor':'red', 'zorder':2}, showfliers=False, 
           whiskerprops={'zorder':2,'color':'red'},capprops={'color':'red'}, medianprops={'zorder':2,'color':'red'})

pairs=list(comb(numint['age_cat'].unique().tolist(),2))
ann=annot(ax=ax,data=numint,x='age_cat',y='NumInt',pairs=pairs)
ann.configure(test='Kruskal',text_format='star',loc='inside')
ann.apply_and_annotate()

ax.set(xlabel='Age category',ylabel='Number of interactions')

plt.savefig('/'.join([wdir,'results/Num_Interactions_AgeCat.pdf']))

#Plot variation between number of interactions vs dmms

ax=sb.swarmplot(numint,y='NumInt',x='dmm',color='#5fc0bf',zorder=1, order=['V1','V2','V3','V4'])
sb.boxplot(numint,y='NumInt',x='dmm',fill=True, width=0.1,fliersize=0,showcaps=True, order=['V1','V2','V3','V4'],
           boxprops={'facecolor':'none', 'edgecolor':'red', 'zorder':2}, showfliers=False,
           whiskerprops={'zorder':2,'color':'red'},capprops={'color':'red'}, medianprops={'zorder':2,'color':'red'})

pairs=list(comb(numint['dmm'].unique().tolist(),2))
ann=annot(ax=ax,data=numint,x='dmm',y='NumInt',pairs=pairs)
ann.configure(test='Kruskal',text_format='star',loc='inside')
ann.apply_and_annotate()

ax.set(xlabel='DMM',ylabel='Number of interactions')

plt.savefig('/'.join([wdir,'results/Num_Interactions_DMM.pdf']))

#---------------------------------------------------------------------------
#Jaccard distance between individuals
#---------------------------------------------------------------------------

def jaccard_d(data):
    dist=pd.DataFrame(index=data.index.to_list(),columns=data.index.to_list())
    for ix,rowa in data.iterrows():
        for iy, rowb in data.iterrows():
            j = jaccard(rowa, rowb)
            dist.at[ix,iy]=j
    
    return dist

jacdist=jaccard_d(prestab)
jacdist.to_csv('/'.join([wdir, 'datasets','Individual_interactions_Jaccard_CrisprF.tsv']), sep='\t',index=True)

##PERMANOVA ANALYSIS

def calc_permanova(skdm,groups,cols):
    
    perm_res=pd.DataFrame()
    
    for col in cols:
        perm=permanova(skdm,groups[col].tolist())
        perm['Grouping']=col
        
        perm_res=pd.concat([perm_res,perm.to_frame().T],ignore_index=True)

    return perm_res
    

cols=['sex','age_cat','beforeBL','center','dmm']

jacdist=pd.read_csv('/'.join([wdir, 'datasets','Individual_interactions_Jaccard_CrisprF.tsv']), sep='\t').set_index('Unnamed: 0')

jacdist=pd.read_csv('/'.join([wdir, 'datasets','dist_by_shared_spacer_targets_jaccard.tsv']), sep='\t')
jacdist['sample_id']=jacdist['sample_id'].apply(lambda row: row.replace('-','_'))
jacdist=jacdist.set_index('sample_id')
jacdist.columns=[c.replace('-','_') for c in jacdist.columns.tolist()]


skdis=DistanceMatrix(np.ascontiguousarray(jacdist.values), ids=jacdist.index.to_list())
groups=pd.DataFrame({'sample_id':jacdist.index.tolist()})
groups=groups.merge(samples[['sample_id','deltaker_id','beforeBL','age_cat','sex','center']], on='sample_id', how='left')
groups=groups.merge(dmm,on='sample_id',how='left')
groups=groups.rename(columns={'gr':'dmm'})

dis_perm=calc_permanova(skdis,groups,cols)

dis_perm.to_csv('/'.join([wdir, 'results/Individual_targets_Jaccard_PERMANOVA_CrisprF.tsv']), sep='\t')

##Plot clustermap grouping the data

col='sex'
labels=groups[col].unique().tolist()
colors=['#a5c0de','#b71c80']
colmap=dict(zip(labels,colors))
groups=groups.loc[groups['sample_id'].isin(jacdist.index.tolist())]
groups['color']=groups[col].apply(lambda row: colmap[row])
groups=groups.sort_values(by=col)
plotmap=groups[['sample_id','color']].set_index('sample_id')['color']

jacdist_p=jacdist.reindex(groups['sample_id'].tolist())
jacdist_p=jacdist_p[groups['sample_id'].tolist()]

sb.clustermap(jacdist_p,row_colors=plotmap,col_colors=plotmap, row_cluster=True,col_cluster=True)

##Jaccard distance within-between the DMM groups

JacDMM=pd.DataFrame()
for d in ['V1','V2','V3','V4']:
    sam=dmm.loc[dmm['gr']==d,'sample_id'].tolist()
    sam=[s for s in sam if s in jacdist.index.tolist()]
   
    for s in sam:
        #compare to samples within the dmm
        j=jacdist.loc[s,[x for x in sam if x!=s]]
        inavg=j.mean()
        #compare to samples outside the dmm
        out=[s for s in jacdist.index.tolist() if s not in sam]
        j=jacdist.loc[s,out]
        utavg=j.mean()
        jacdmm=pd.DataFrame({'sample_id':[s],'dmm':[d],'SelfDMMJac':[inavg],'OutDMMJac':[utavg]})
        JacDMM=pd.concat([JacDMM,jacdmm])
        
JacDMM=pd.melt(JacDMM,id_vars=['sample_id','dmm'],value_name='AvgJac',var_name='DMMgroup')

ax=sb.violinplot(JacDMM,x='dmm',y='AvgJac',hue='DMMgroup',order=['V1','V2','V3','V4'],
                 legend=False, palette=['#5fc0bf','#b9b9b9'],split=True,gap=.2,cut=0,
                 inner='quart')      

ji=['SelfDMMJac','OutDMMJac']   
pairs=[((d,ji[0]),(d,ji[1])) for d in ['V1','V2','V3','V4']]
ann=annot(ax=ax,data=JacDMM,x='dmm',y='AvgJac',pairs=pairs,hue='DMMgroup')
ann.configure(test='Kruskal',text_format='star',loc='outside',verbose=True)
ann.apply_and_annotate()

ax.set(ylabel='Average Jaccard index',xlabel='')
plt.savefig('/'.join([wdir,'results/AvgJaccard_DMM_SelfvsOut.pdf']))
