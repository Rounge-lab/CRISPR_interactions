#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 14:59:23 2024

Beta-diversity calculations 

@author: ekateria
"""

import pandas as pd
import pyreadr as pyr
import numpy as np
from scipy.spatial.distance import braycurtis, jaccard
import matplotlib.pyplot as plt


from skbio import DistanceMatrix
from skbio.stats.distance import permanova
from skbio.stats.ordination import pcoa


wdir='PATH_TO_MANUS_FOLDER' #set working directory
prefix='MAGs' #prefix to be used for the output files
#prefix='vOTUs' #prefix to be used for the output files
#prefix='pOTUs'
#Load metadata--------------------------------------------------------------

####NB that beforeBL variable codes for a/b use 4 months prior to baseline, not if it was at ANY time before baseline  

meta=pd.read_csv('/'.join([wdir,'participant_data/screening_data.tsv']), sep='\t')
meta['beforeBL'] = pd.Categorical(meta['beforeBL'], categories=['Yes','No'], ordered=True)
meta['kjonn'] = pd.Categorical(meta['kjonn'], categories=['Male','Female'], ordered=True)
meta['age_cat'] = pd.Categorical(meta['age_cat'], categories=['50-59','60-69','>=70'], ordered=True)

samples=pyr.read_r('/'.join([wdir, 'participant_data/sample_meta.Rds']))[None]
samples=samples.merge(meta[['deltaker_id','beforeBL','age_cat','kjonn','senter']],on='deltaker_id',how='left')

#------------------------------------------------------------------------------------------------------------
#BETA DIVERSITY
#Calculate beta diversity (BC)

def bray_curtis(data):
    dist=pd.DataFrame(index=data.index.to_list(),columns=data.index.to_list())
    for ix,rowa in data.iterrows():
        for iy, rowb in data.iterrows():
            bc = braycurtis(rowa, rowb)
            dist.at[ix,iy]=bc
    
    return dist

def calc_prev(relab, domain):
    prev=relab.copy()
    if domain!='MAGs':
        prev[prev>0]=1 #keep all results over 0 since each virus/plasmid are covered on at least 75% of its length
    else:
        prev[prev>=0.001]=1 #set the threshold to 0.001 % (it makes 10000 bp out of a minimum threshold of 1Gb data per sample)
        prev[prev<0.001]=0
    
    return prev


def jaccard_d(data):
    dist=pd.DataFrame(index=data.index.to_list(),columns=data.index.to_list())
    for ix,rowa in data.iterrows():
        for iy, rowb in data.iterrows():
            j = jaccard(rowa, rowb)
            dist.at[ix,iy]=j
    
    return dist


#------------------------ 
#Load data and calculate beta-div
for p in ['MAGs','vOTUs','pOTUs', 'Combined']:
    
    if p=='Combined':
        mags=pd.read_csv('/'.join([wdir, 'datasets','MAGs_relab.tsv']), sep='\t').set_index('sample_id')
        votu=pd.read_csv('/'.join([wdir, 'datasets','vOTUs_relab.tsv']), sep='\t').set_index('sample_id')
        potu=pd.read_csv('/'.join([wdir, 'datasets','pOTUs_relab.tsv']), sep='\t').set_index('sample_id')

        relab=mags.merge(votu, left_index=True, right_index=True)
        relab=relab.merge(potu, left_index=True, right_index=True)
        
        prev=calc_prev(mags,'MAGs')
        vprev=calc_prev(votu,'vOTUs')
        pprev=calc_prev(potu,'pOTUs')
        
        prev=prev.merge(vprev,left_index=True,right_index=True)
        prev=prev.merge(pprev,left_index=True,right_index=True)

    else:
        relab=pd.read_csv('/'.join([wdir, 'datasets', p+'_relab.tsv']), sep='\t').set_index('sample_id')
        prev=calc_prev(relab,p)
        

    #Calculate BC and Jaccard
    print(f'Calculating BrayCurtis for {p}...')
    bray=bray_curtis(relab)
    print(f'Calculating Jaccard for {p}...')
    jac=jaccard_d(prev)
    
    #Save data
    bray.to_csv('/'.join([wdir, 'datasets', p+'_BrayCurtis.tsv']), sep='\t')
    jac.to_csv('/'.join([wdir, 'datasets', p+'_Jaccard.tsv']), sep='\t')

#Metagenome_bray (mags, votus and potus combined)

##PERMANOVA ANALYSIS

def calc_permanova(skdm,groups,cols):
    
    perm_res=pd.DataFrame()
    
    for col in cols:
        perm=permanova(skdm,groups[col].tolist())
        perm['Grouping']=col
        
        perm_res=pd.concat([perm_res,perm.to_frame().T],ignore_index=True)
    return perm_res
    
distances=['BrayCurtis','Jaccard']
cols=['kjonn','age_cat','beforeBL','senter']

for p in ['MAGs','vOTUs','pOTUs', 'Combined']:
    for d in distances:
        print(f'Performing PERMANOVA for {p}_{d}')
        dis=pd.read_csv('/'.join([wdir, 'datasets',f'{p}_{d}.tsv']), sep='\t').set_index('Unnamed: 0')
        skdis=DistanceMatrix(np.ascontiguousarray(dis.values), ids=dis.index.to_list())
    
        groups=pd.DataFrame({'sample_id':dis.index.tolist()})
        groups=groups.merge(samples[['sample_id','deltaker_id','beforeBL','age_cat','kjonn','senter']], on='sample_id', how='left')
    
        dis_perm=calc_permanova(skdis,groups,cols)
        dis_perm.to_csv('/'.join([wdir, f'results/{p}_PERMANOVA_{d}.tsv']), sep='\t')

# Perform PCoA and make the plot
def pcoa_analysis(distancematrix, groups,col):
    
    pcoa_res = pcoa(distancematrix)

    # Visualize PCoA results
    sc = pcoa_res.samples
    sc = sc[['PC1','PC2','PC3']].merge(groups,left_index=True,right_index=True)
    
    var_exp_pc1 = pcoa_res.proportion_explained[0]*100
    var_exp_pc2 = pcoa_res.proportion_explained[1]*100

    plt.figure()

    gr=sc[col].unique().tolist()
    paints=['#9791c6', '#1c6462','#77bdc2','#5fc0bf','#eac0bf','#ace1a5','#b2b3cd']
       
    # Plot PC1, PC2 values for each sample, calculate and plot centroids
    for g in range(len(gr)):
        centroid = sc.loc[sc[col]==gr[g],['PC1','PC2']].mean(axis=0)
        plt.scatter(sc.loc[sc[col]==gr[g],'PC1'], sc.loc[sc[col]==gr[g],'PC2'], label=gr[g],color=paints[g])
        plt.scatter(centroid['PC1'], centroid['PC2'], color=paints[g], marker='^',s=200, edgecolor='k')
        
    plt.legend(title=col)
    plt.xlabel(f'PC1, {var_exp_pc1:.1f}% variance explained')
    plt.ylabel(f'PC2, {var_exp_pc2:.1f}% variance explained')
    plt.show()
    
    return pcoa_res,sc

p='Combined'
d='BrayCurtis'
dis=pd.read_csv('/'.join([wdir, 'datasets',f'{p}_{d}.tsv']), sep='\t').set_index('Unnamed: 0')
skdis=DistanceMatrix(np.ascontiguousarray(dis.values), ids=dis.index.to_list())
pcoa_bray,sc_bray=pcoa_analysis(skdis,groups,'beforeBL')
plt.savefig('/'.join([wdir, 'PCoA_Combined_BrayBeforeBL.pdf']))


#Post hoc permanova 

def calc_permanova_posthoc(distmat,groups,col,vsgroup):
    
    gr=groups[col].unique().tolist()
    gr=[g for g in gr if g!=vsgroup]
    perm_res=pd.DataFrame()
 
    for g in gr:
        
        ids=groups.loc[groups[col]==g]
        ids=pd.concat([ids, groups.loc[groups[col]==vsgroup]],ignore_index=False)
        dm=distmat.loc[ids.index.tolist(),ids.index.tolist()]
        dm=DistanceMatrix(np.ascontiguousarray(dm.values), ids=dm.index.to_list())
        perm=permanova(dm,ids[col].tolist())
        perm['Target']=g
        perm['Contrast']=vsgroup
        
        perm_res=pd.concat([perm_res,perm.to_frame().T],ignore_index=True)
    return perm_res

phperm=calc_permanova_posthoc(bray,groups,'age_cat','60-69')
phperm.to_csv('/'.join([wdir, 'PostHoc_PERMANOVA_Combined_Bray_Agecat.tsv']), sep='\t')

