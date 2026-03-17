#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 14:14:26 2025

Connectivity between MAGs

@author: ekateria
"""

import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
from itertools import combinations
from scipy.stats import chi2_contingency

wdir='PATH_TO_MANUS_FOLDER'

nodes=pd.read_csv('/'.join([wdir,'results/Nodes_network_crisprfree.csv']))
edges=pd.read_csv('/'.join([wdir,'results/Edges_network_crisprfree.csv']))

##Make a MAG-MAG clustermap

def count_distance(edges, nodes=nodes):
    maglist=edges['MAG'].unique().tolist()
    dist=pd.DataFrame(data=0,index=maglist, columns=maglist)
    for link in edges['Taxon'].unique().tolist():
        query=edges.loc[edges['Taxon']==link]
        if len(query)>1:
            mcons=list(combinations(query['MAG'].unique().tolist(),2))
            for mc in mcons:
                dist.loc[mc[0],mc[1]]+=1
                dist.loc[mc[1],mc[0]]+=1 #make square matrix
                
    magnodes=nodes.loc[nodes['Taxon'].isin(maglist),['Taxon','Label']]
    
    return dist,magnodes

alldist,magnodes=count_distance(edges)

#Makecolmap only once for all the combinations (both votus and potus)
labels=magnodes['Label'].unique().tolist()
colors=['#fd8d3c','#bdbdbd','#8c96c6', '#ae017e','#edf8b1', '#02818a','#666600','#9292ff',
        '#ff92ff', '#f781bf','#b4639f', '#74c476','#d7301f','#b2182b','#ffff99','#92ffff']
colmap=dict(zip(labels,colors))

def plot_clusters(plotdist, clustbool, magnodes=magnodes, colmap=colmap):
    plotnodes=magnodes.loc[magnodes['Taxon'].isin(plotdist.index.tolist())]
    plotnodes['family']=plotnodes['Label'].apply(lambda row: colmap[row])
    plotcols=plotnodes[['Taxon','family']].set_index('Taxon')['family']
    
    fig=sb.clustermap(data=plotdist, row_colors=plotcols, col_colors=plotcols, 
                      cmap=sb.color_palette('light:#5A9'),row_cluster=clustbool, col_cluster=clustbool)
    fig.ax_heatmap.set_xticklabels('')
    fig.ax_heatmap.set_yticklabels('')
    fig.ax_heatmap.tick_params(right=False, bottom=False)
    
    if not clustbool:
        fig.cax.set_position([0.1, 0.2, 0.02, 0.3])
    
#Filter based on the number of connections
tokeep=alldist.sum(axis=0)
tokeep=tokeep[tokeep>=200]
plotdist=alldist.loc[tokeep.index.tolist(), tokeep.index.tolist()]
plot_clusters(plotdist,True)
plt.savefig('/'.join([wdir,'results/MAG-to-MAG_connectivity/Connections_atleast200_targetpair.pdf']))

#Make a presence/absence of connection
plotdist=alldist.copy()
plotdist[plotdist>0]=1
plot_clusters(plotdist,True)
plt.savefig('/'.join([wdir,'results/MAG-to-MAG_connectivity/Presence_absence_targetpair.pdf']))

#vOTUs
virdist,_=count_distance(edges.loc[edges['Taxon'].str.contains('vOTU')])
#pOTUs
plasdist,_=count_distance(edges.loc[edges['Taxon'].str.contains('plasmid')])

type='pOTU'
qdist=plasdist.copy()
#number of connections
tokeep=qdist.sum(axis=0)
tokeep=tokeep[tokeep>=200]
plot_clusters(qdist.loc[tokeep.index.tolist(),tokeep.index.tolist()],True)
plt.savefig('/'.join([wdir,f'results/MAG-to-MAG_connectivity/Connections_{type}_only_at_least200_targetpairs.pdf']))

#presence/absence of connection
plotdist=qdist.copy()
plotdist[plotdist>0]=1
tokeep=qdist.sum(axis=0)
tokeep=tokeep[tokeep>=1]
plot_clusters(plotdist.loc[tokeep.index.tolist(),tokeep.index.tolist()], True)
plt.savefig('/'.join([wdir,f'results/MAG-to-MAG_connectivity/Presence_absence_{type}_only.pdf']))

#Arrange by taxonomy and make a heatmap
#by number of connections
tokeep=qdist.sum(axis=0)
tokeep=tokeep[tokeep>=100] #200 for all and vOTU, 100 for plasmids
plotdist=qdist.loc[tokeep.index.tolist(), tokeep.index.tolist()]
plotnodes=magnodes.loc[magnodes['Taxon'].isin(plotdist.index.tolist())]
plotnodes=plotnodes.sort_values(by='Label', ascending=True)
plotdist=plotdist.loc[plotnodes['Taxon'].tolist(),plotnodes['Taxon'].tolist()]
plot_clusters(plotdist,False)
plt.savefig('/'.join([wdir,'results/MAG-to-MAG_connectivity/Heatmap_Connections_all_at_least200_targetpairs.pdf']))


#by presence/absence of connections
sb.color_palette("light:#5A9", as_cmap=True)
plotdist=qdist.copy()
plotdist[plotdist>0]=1
tokeep=plotdist.sum(axis=0)
tokeep=tokeep[tokeep>=1]
plotdist=plotdist.loc[tokeep.index.tolist(),tokeep.index.tolist()]
plotnodes=magnodes.loc[magnodes['Taxon'].isin(plotdist.index.tolist())]
plotnodes=plotnodes.sort_values(by='Label', ascending=True)
plot_clusters(plotdist.loc[plotnodes['Taxon'].tolist(),plotnodes['Taxon'].tolist()], False)
plt.savefig('/'.join([wdir,'results/MAG-to-MAG_connectivity/Heatmap_Presence_absence_pOTU_only.pdf']))


##Make a boxplot for MAG connectivity to same family/other families

magnodes=magnodes.rename(columns={'Taxon':'MAG'})
edges=edges.merge(magnodes[['MAG','Label']], on='MAG', how='left')

def connectivity(edges):
    connect=edges[['MAG','Label']].drop_duplicates()
    connect['Same_family']=0
    connect['Single_hit']=0
    connect['Other_family']=0
    
    for ix, m in connect.iterrows():
        print(f'Checking connections for {m.MAG}...')
        q=edges.loc[edges['MAG']==m.MAG]
        for link in q['Taxon'].tolist():
            qt=edges.loc[edges['Taxon']==link]
            qt=qt.loc[qt['MAG']!=m.MAG]
            if len(qt)>0:
                for f in qt['Label'].tolist():
                    if f!=m.Label:
                        connect.loc[ix,'Other_family']+=1
                    else:
                        connect.loc[ix,'Same_family']+=1
            else:
                connect.loc[ix,'Single_hit']+=1
    
    connect_melt=pd.melt(connect,id_vars=['MAG','Label'],var_name='ConnectionType',value_name='NumConnections')
    connect['SameFraction']=connect['Same_family']/(connect['Same_family']+connect['Other_family'])*100

    return connect, connect_melt

virconnect,virconnect_melt=connectivity(edges.loc[edges['Taxon'].str.contains('vOTU')])  
virconnect['Target']='vOTU'
virconnect_melt['Target']='vOTU'
plasconnect, plasconnect_melt=connectivity(edges.loc[edges['Taxon'].str.contains('plasmid')])      
plasconnect['Target']='pOTU'
plasconnect_melt['Target']='pOTU'


allconnect=pd.concat([virconnect, plasconnect])
allconnect.to_csv('/'.join([wdir, 'results/MAG_connectivity_crisprfree.csv']), index=False, sep='\t')

#Order the data for plotting
order=allconnect.loc[allconnect['Target']=='vOTU',['Label','SameFraction']].groupby('Label').median().sort_values(by='SameFraction',ascending=False)
cats=order.index.tolist()
cats=[c for c in cats if c!='Other']
cats.append('Other')
allconnect['Label']=pd.Categorical(allconnect['Label'], categories=cats, ordered=True)
fig=sb.catplot(allconnect, kind='box', hue=None, x='SameFraction',y='Label',col='Target', col_order=['vOTU','pOTU'], color='#5fc0bf')
fig.set(xlabel='Same family MAGs among connected, %',ylabel='')
for ax in fig.axes.flat:
    for label in ax.get_yticklabels():
        label.set_fontstyle('italic')
        
plt.savefig('/'.join([wdir,'results/Connectivity_MAGs_boxplots.pdf']))

#Find how many MAGs have single connections 
single=allconnect.loc[allconnect['Single_hit']!=0]
single['all']=single['Same_family']+single['Single_hit']+single['Other_family']
single=single.loc[single['all']==single['Single_hit']]
sd=single['Target'].value_counts() #130 MAGs target individual vOTUs and 74 - pOTUs
num_shared=len(allconnect['MAG'].unique().tolist())-len(single['MAG'].unique().tolist())
print(f'Number of MAGs that only has individual vOTU targets: {sd.vOTU}')
print(f'Number of MAGs that only has individual pOTU targets: {sd.pOTU}')
print(f'Number of MAGs that has shared targets: {num_shared}')


## Do conjugative plasmid tend to connect many families vs other plasmid types

plasmob=pd.read_csv('/'.join([wdir, 'datasets/plasmids_dereplicated_0.9_MOBtyper.txt']), sep='\t')
plasmob=plasmob.rename(columns={'sample_id':'Taxon'})

plasmob=plasmob.loc[plasmob['Taxon'].isin(edges['Taxon'].unique().tolist()),['Taxon','predicted_mobility']]

nodes=nodes.rename(columns={'Taxon':'MAG'})
edges=edges.merge(nodes[['MAG','Taxonomy']], on='MAG', how='left') #Add taxonomy to avoid that if there are several targets to 'Other' family, that they are listed as one

for ix, p in plasmob.iterrows():
    query=edges.loc[edges['Taxon']==p.Taxon]
    plasmob.at[ix,'NumConnectedFamilies']=len(query['Taxonomy'].unique().tolist())

plasmob['NumFamCat']=plasmob['NumConnectedFamilies'].apply(lambda row: 'intrafamily' if row==1 else 'interfamily')

#####Calculate fraction of multiple families vs one family for each mobility class; binomial p-value afterwards
plasfam=plasmob.groupby('predicted_mobility')['NumFamCat'].value_counts().reset_index()
plasfam=plasfam.pivot(index='predicted_mobility',columns='NumFamCat',values='count')
plasfam['Fraction']=plasfam['interfamily']/(plasfam['intrafamily']+plasfam['interfamily'])*100

chi2,p,dof,expected=chi2_contingency(plasfam[['interfamily','intrafamily']])
chi2,p,dof,expected=chi2_contingency(plasfam.loc[['mobilizable','conjugative'],['interfamily','intrafamily']])

##Do provirus/virus tend to connect many/single families
votus=pd.read_csv('/'.join([wdir, 'datasets/viral_contigs_list_organized_onlymanus_votus.csv']), sep=',')

votu_connect=votus[['new_id']].drop_duplicates()

for ix, v in votu_connect.iterrows():
    query=edges.loc[edges['Taxon']==v.new_id]
    checkV=votus.loc[votus['new_id']==v.new_id,'checkV'].unique().tolist()
    if 'Yes' in checkV:
        votu_connect.loc[ix,'checkV']='Yes'
    else:
        votu_connect.loc[ix,'checkV']='No'
    votu_connect.at[ix,'NumConnectedFamilies']=len(query['Taxonomy'].unique().tolist())
    

votu_connect['NumFamCat']=votu_connect['NumConnectedFamilies'].apply(lambda row: 'intrafamily' if row==1 else 'interfamily')

#####Calculate fraction of multiple families vs one family for each mobility class; binomial p-value afterwards
votufam=votu_connect.groupby('checkV')['NumFamCat'].value_counts().reset_index()
votufam=votufam.pivot(index='checkV',columns='NumFamCat',values='count')
votufam['Fraction']=votufam['interfamily']/(votufam['intrafamily']+votufam['interfamily'])*100

chi2,p,dof,expected=chi2_contingency(votufam[['interfamily','intrafamily']])

# Summarize fraction of MAGs within mOTUs with the same target

mges=pd.read_csv('/'.join([wdir,'datasets/targeting_genomes_per_target.tsv']), sep='\t')
motus_single_target=mges.loc[mges['n_cassette_genomes']>=10]

# Make a list of species that are potential species-inferred hosts

magtax=pd.read_csv('/'.join([wdir, '/datasets/MAG_taxonomy_full.tsv']), sep=',')

def add_magsp(level,magtax=magtax,edges=edges,nodes=nodes):
    mge=nodes.loc[nodes['Domain']==level,'Taxon'].to_frame()
    mge=mge.drop_duplicates()
    mge.columns=[level]
    for ix, m in mge.iterrows():
        query=edges.loc[edges['Taxon']==m[level],'MAG'].tolist()
        query=list(set(query)) #keep only unique entries
        species=magtax.loc[magtax['MAG'].isin(query),'species'].fillna('unknown').tolist()
        species=list(set(species))
        mge.at[ix,'TargetedBySpecies']=';'.join(species)
        
    return mge

virtarget=add_magsp('vOTU')
plastarget=add_magsp('pOTU')

virtarget.to_csv('/'.join([wdir, '/results/vOTUs_targeted_by_species.csv']),sep='\t', index=False)
plastarget.to_csv('/'.join([wdir, '/results/PTUs_targeted_by_species.csv']),sep='\t', index=False)

    
    
        

