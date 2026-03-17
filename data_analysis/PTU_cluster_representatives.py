#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Find cluster sizes for each PTU

"""

import pandas as pd
import pyreadr as pyr
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt

wdir='PATH_TO_MANUS_FOLDER' #set working directory

clustermap=pd.read_csv('PATH_TO/scapp_snakemake_wf/data/dereplication/cluster_map.txt', 
                       sep='\t',header=None)

clustermap.columns=['genome','pOTU','Identity','Length','Plasmid length','Evalue']

clustermap.to_csv('PATH_TO/scapp_snakemake_wf/data/dereplication/cluster_map.txt', 
                       sep='\t',index=False)

num_plasmids=clustermap.groupby(['pOTU'])['genome'].count().reset_index()
num_plasmids.to_csv('/'.join([wdir,'datasets/Plasmid_NumGenomes.csv']),sep='\t',index=False)
clustermap.to_csv('/'.join([wdir,'datasets/pOTU_Clustermap.csv']),sep='\t',index=False)

##get the statistics on number of plasmids detected per sample (by SCAPP)

numplas=pd.read_csv('PATH_TO/scapp_snakemake_wf/data/dereplication/cluster_map.txt',sep='\t')
relab=pd.read_csv('/'.join([wdir, 'datasets/pOTUs_relab.tsv']),sep='\t')
relab=relab.set_index('sample_id')

#Keep only pOTUs that are in the manus
numplas=numplas.loc[numplas['pOTU'].isin(relab.columns.tolist())]
numplas['sample']=numplas['genome'].apply(lambda row: row.split('_')[0])
persam=numplas['sample'].value_counts().reset_index()
persam['sample']=persam['sample'].apply(lambda row: row.replace('-','_'))
persam=persam.loc[persam['sample'].isin(relab.index.tolist())]
persam['count'].describe()
persam.to_csv('/'.join([wdir,'results/NumPlasmids_assembled_persample.csv']),index=False)

#Check correlations between number of plasmids and sequencing depth
samples=pyr.read_r('/'.join([wdir, 'participant_data/sample_meta.Rds']))[None]
persam=persam.rename(columns={'sample':'sample_id','count':'Number assembled plasmids'})
persam=persam.merge(samples[['sample_id','Total_Bases_QC_ATLAS','N50']], on='sample_id',how='left')

mag_alpha=pd.read_csv('/'.join([wdir, 'datasets/MAGs_AlphaDiv.tsv']), sep='\t')
persam=persam.merge(mag_alpha[['sample_id','ObsSp','InvSimpson','Shannon']], on='sample_id')

corrtab=persam.set_index('sample_id').corr(method='spearman')

# Plot the correlation tab
mask = np.tril(np.ones_like(corrtab, dtype=bool))  # Lower triangle mask (True = hide)
fig=sb.heatmap(corrtab,mask=mask,cmap="coolwarm", annot=True, square=True,  cbar_kws={"shrink": 0.8})
fig.set_xticklabels(fig.get_xticklabels(),rotation=90)
plt.savefig('/'.join([wdir,'results/PlasmidAssembly_DepthCorrelation.pdf']))