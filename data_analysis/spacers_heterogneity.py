#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Spacers heterogeneity - in how many different mOTUs are they detected

"""
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator
from itertools import combinations as comb
from statsmodels.stats.multitest import multipletests
import scipy.stats as stats
from scipy.stats import binomtest

#Load all needed data
wdir="PATH_TO_MANUS_FOLDER"

spacers=pd.read_csv('/'.join([wdir,'datasets/spacers_manus_table_new.csv']),sep='\t')

samples=pd.read_csv('/'.join([wdir, 'datasets/data_by_samples_crispr.csv']),sep = '\t')
samples=samples.rename(columns={'sample_id':'Sample'})

#Find the size of the clusters

clustmap=spacers[['Cluster','Spacers']].drop_duplicates(subset='Spacers')
clustmap=clustmap['Cluster'].value_counts().reset_index()
clustmap=clustmap.rename(columns={'count':'NumSpacers'})

clustind=spacers[['Cluster','Sample']].drop_duplicates()
clustind=clustind['Cluster'].value_counts().reset_index()
clustind=clustind.rename(columns={'count':'NumIndividuals'})

clustmap=clustmap.merge(clustind, on='Cluster',how='left')

sb.histplot(clustmap, x='NumSpacers')
sb.lineplot(clustmap, x='NumIndividuals', y='NumSpacers')

clustmap.to_csv('/'.join([wdir, 'results/SpacerCluster_Sizes.csv']),sep='\t',index=False)

#Find how many are singleton/non-singleton

singleton=clustmap.query('NumSpacers==1')
non_singleton=clustmap.query('NumSpacers>1')
non_singleton=non_singleton.merge(clustind, on='Cluster', how='left')
len(non_singleton.query('NumIndividuals>1'))

#Check if those clusters that are not singletons, are detected in same/different MAGs

shared=clustmap.query('NumSpacers>1')
shared=shared.merge(spacers[['Cluster','Spacers','MAG']], on='Cluster',how='left')
shared=shared.drop_duplicates().dropna(subset='MAG')
sharedsum=shared.groupby(['Cluster','MAG'])['MAG'].value_counts().reset_index()
sharedsum.to_csv('/'.join([wdir,'results/SpacersHeterogeneity_MAGs.csv']),sep='\t', index=False)

#For each cluster, find heterogeneity of its spacers
clusthet=sharedsum['Cluster'].value_counts().reset_index().rename(columns={'count':'NumMAGs'})
clustsum=sharedsum.groupby(['Cluster'])['count'].sum().reset_index().rename(columns={'count':'NumSpacersInMAGs'})
clusthet=clusthet.merge(clustsum, on='Cluster',how='left')
clusthet=clusthet.merge(clustmap[['Cluster','NumSpacers']], on='Cluster', how='left')
clusthet=clusthet.rename(columns={'NumSpacers':'NumSpacersTotal'})
clusthet.to_csv('/'.join([wdir,'results/SpacersHeterogeneity_MAGs_summary.csv']),sep='\t', index=False)
clusthet['NumMAGsScaled']=clusthet['NumMAGs']*10

###for those clusters that are detected in >1 mOTU, check LCA of these motus
magtax=pd.read_csv('/'.join([wdir,'datasets/MAG_taxonomy_full.tsv']),sep=',')

def check_lca(reps,magtax=magtax):
    qtax=magtax.loc[magtax['MAG'].isin(reps)]
    cols=qtax.columns.tolist()
    cols=cols[:0:-1]
    stat='no common'
    for c in cols:
        if len(qtax[c].unique().tolist())==1:
            stat=qtax[c].unique().tolist()[0]
            break
    return c,stat

def check_family(reps,magtax=magtax):
    qtax=magtax.loc[magtax['MAG'].isin(reps)]
    if len(qtax['family'].unique().tolist())==1:
        stat='same family'
    return stat

for ix,r in clusthet.iterrows():
    if r.NumMAGs>1:
        c,stat=check_lca(shared.loc[shared['Cluster']==r.Cluster,'MAG'].unique().tolist())
        clusthet.at[ix,'MAG_LCA']=stat
        clusthet.at[ix,'MAG_LCAlevel']=c
        
        c=check_family(shared.loc[shared['Cluster']==r.Cluster,'MAG'].unique().tolist())
        clusthet.at[ix,'Same_family']=stat
      
clusthet.query('NumMAGs>1')['MAG_LCAlevel'].value_counts()            
clusthet.query('NumMAGs>1')['Same_family'].value_counts()       

ax=sb.scatterplot(clusthet.query('NumMAGs<3'), x='NumSpacersTotal', y='NumSpacersInMAGs',size='NumMAGs',sizes=(10,50), color='#6f748880', legend=False)
sb.scatterplot(clusthet.query('NumMAGs>=3'), x='NumSpacersTotal', y='NumSpacersInMAGs',size='NumMAGs',sizes=(100,300), hue='Same_family', palette=['#5fc0bf80','#e2b5b580'])
sb.lineplot(x=[0, 260], y=[0, 260], linestyle="--", color='#bbbbbb')
ax.set(xlabel='Number of spacers in a cluster', ylabel='Number of spacers in MAGs')
plt.savefig('/'.join([wdir,'results/SpacersHeterogeneity.pdf']))

########
#get clusters that are detected in over 10% of the population (n=103 individuals)

prevclus=clustmap.query('NumIndividuals>=103')
prevclus=prevclus.merge(spacers[['Cluster','MAG']], on='Cluster', how='left').drop_duplicates()
knownmag=prevclus[['Cluster','MAG']].dropna(subset='MAG').drop_duplicates()
magtax=pd.read_csv('/'.join([wdir,'datasets/MAG_taxonomy_full.tsv']),sep=',')
knownmag=knownmag.merge(magtax,on='MAG',how='left')
prevclusbymag=knownmag.groupby(['Cluster','MAG'])['species'].value_counts().reset_index()

#####
#Assess spacer heterogeneity in mOTUs
#####
mags=spacers[['MAG','Cluster','Sample']].drop_duplicates()
nummags=mags[['MAG','Sample']].drop_duplicates()
nummags=nummags['MAG'].value_counts().reset_index()
nummags=nummags.rename(columns={'count':'NumMAGsWithSpacers'})

motus=mags[['MAG','Cluster']].groupby(['MAG','Cluster']).value_counts().sort_values(ascending=False).reset_index()

motus=motus.rename(columns={'count':'DetectedInNumMAGs'})

motus=nummags.merge(motus, on='MAG', how='left')
motus=motus.sort_values(by='DetectedInNumMAGs', ascending=False)
motus['SpacerPrev']=motus['DetectedInNumMAGs']/motus['NumMAGsWithSpacers']*100
motus=motus.merge(magtax[['MAG','family']], on='MAG', how='left')

#Plot the prevalence of spacers in mOTUs with >10 MAGs with spacers detected 
toplot=motus.query('NumMAGsWithSpacers>=10')

#Order plotting
famorder=toplot[['MAG','family']].drop_duplicates()
famorder=famorder['family'].value_counts().sort_values(ascending=False).reset_index()

famorder=famorder['family'].tolist()[0:10]
famorder.append('others')

toplot['family']=toplot['family'].apply(lambda row: 'others' if row not in famorder else row)

magorder=[]
for f in famorder:
    #find all mags within the group
    fmag=toplot.loc[motus['family']==f]
    #order them based on their median prevalence
    fmag=fmag.groupby(['MAG'])['SpacerPrev'].median().sort_values(ascending=False)
    magorder.extend(fmag.index.tolist())

colors=['#b71c80','#edf4b9','#007588','#fd953f','#8c9ac7','#fed52b','#6eb872','#9399ff','#646104','#eac0bf','#bebebe']
colors=dict(zip(famorder,colors))

fig=sb.boxplot(toplot, x='MAG',y='SpacerPrev',order=magorder,hue='family',palette=colors,
               flierprops=dict(marker='o', markersize=1, markerfacecolor='black', markeredgecolor='black', alpha=0.7))
plt.xticks(rotation=90)
labels=[t.get_text() for t in fig.get_xticklabels()]
labels=[t.replace('MAG','mOTU_') for t in labels]
fig.set_xticklabels(labels)
plt.savefig('/'.join([wdir,'results/mOTU_spacers_heterogenity.pdf']))

#Perform Kruskal_Wallis test with FDR correction for each pair of the MAGs

def kruskal_group(data,cat_col,y):
    query=data.copy()
    query=query.dropna(subset=[cat_col,y])
    categ=query[cat_col].unique().tolist()
    for_stats={}
    for cat in categ:
        for_stats[cat]=query.loc[query[cat_col]==cat,y].tolist()

    h_st,pval=stats.kruskal(*list(for_stats.values())) # *indicates to treat each element in list as a group, the list will contain as many elements as there are dict entries, each entry is a separate list within a list

    return h_st,pval

pairs=list(comb(magorder,2))

KruskalMAGs=pd.DataFrame()
for p in pairs:
    h_st,pval=kruskal_group(toplot.loc[toplot['MAG'].isin(list(p))],'MAG','SpacerPrev')
    KruskalMAGs=pd.concat([KruskalMAGs, pd.DataFrame({'MAG1':p[0],'MAG2':p[1],'Hst':[h_st],'Pval':[pval]})])

_,padj,_,_=multipletests(KruskalMAGs['Pval'],method='fdr_bh')
KruskalMAGs['FDRp']=padj

KruskalMAGs=KruskalMAGs.sort_values(by='FDRp',ascending=True)

KruskalMAGs.to_csv('/'.join([wdir,'results/mOTU_Spacers_heterogenity_Kruskal.csv']),index=False)

fig=sb.boxplot(toplot, x='family',y='SpacerPrev',order=famorder,color='#5fc0bf',
               flierprops=dict(marker='o', markersize=1, markerfacecolor='black', markeredgecolor='black', alpha=0.7))
plt.xticks(rotation=90)
pairs=list(comb(famorder,2))

annotator = Annotator(fig, pairs, data=toplot, x='family', y='SpacerPrev', order=famorder)
annotator.configure(test='Kruskal', text_format='star', comparisons_correction='fdr_bh',
                    correction_format='replace',loc='inside', hide_non_significant=True,
                    line_height=0.01,line_offset_to_group=0.01,line_offset=0.0, line_width=0.5)
annotator.apply_and_annotate()

##Find B.animalis spacer-samples
banim=spacers.loc[spacers['MAG']=='mOTU1812']
banim=banim['participant_id'].drop_duplicates().reset_index()
banim=banim.rename(columns={'participant_id':'deltaker_id'})

diet=pd.read_csv('PATH_TO_FFQ_DATA/foodgroups_n1616.csv', sep=';')

banim=banim.merge(diet[['deltaker_id','CULTB']], on='id',how='left')
banim.to_csv('/'.join([wdir,'results/Banimalis_ids_Biola.csv']),sep='\t',index=False)

col='CULTB'
yoghurt=diet[['id',col]]
deltid=deltid.loc[deltid['deltaker_id'].isin(spacers['participant_id'].unique().tolist())]
yoghurt['banim']=yoghurt['id'].apply(lambda row: 'Yes' if row in banim['id'].tolist() else 'No')
yoghurt=yoghurt.loc[yoghurt['id'].isin(deltid['id'].tolist())]

#remove those for who we don't have metagenomic data

sb.boxplot(yoghurt, x='CULTB',y='banim', color='#5fc0bf')

h_st,pval=kruskal_group(yoghurt,'banim',col)

yoghurt['bin']=yoghurt[col].apply(lambda row:'ja' if row>0 else 'nei')

yoghurt.to_csv('/'.join([wdir,'results/CultBiola_all_ids.csv']),sep='\t',index=False)

cr=pd.crosstab(yoghurt['bin'],yoghurt['banim'])

binomtest(cr.loc['ja','Yes'],33,p=cr.loc['ja','No']/cr.sum(axis=0)[0],alternative='two-sided')
