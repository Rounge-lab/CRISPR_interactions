#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

PTU taxonomy delineation based on Spearman correlations

"""

import pandas as pd
# import scipy.stats as stats
# from statsmodels.stats.multitest import multipletests
import seaborn as sb
import matplotlib.pyplot as plt

wdir='PATH_TO_MANUS_FOLDER' #set working directory

potu_taxon=pd.read_csv('/'.join([wdir, 'datasets/pOTU_taxonomy_IMGPR_0.9_modified.txt']),sep=',')
mag_taxon=pd.read_csv('/'.join([wdir, 'datasets/MAG_taxonomy_full.tsv']), sep=',')
potu_corr=pd.read_csv('/'.join([wdir, 'results/Relab_Spearman_correlation_MAGs_pOTUs_square.csv']),sep=',')
potus=pd.read_csv('/'.join([wdir,'datasets/pOTUs_db_summary.csv']))

potu_corr=potu_corr.rename(columns={'Unnamed: 0':'MAG'})
potu_corr=potu_corr.set_index('MAG')

#Keep only potus with known family delineation and those that are in baseline samples
potu_taxon=potu_taxon.loc[potu_taxon['pOTU'].isin(potus['pOTU'].tolist())]
known=potu_taxon.loc[~potu_taxon['Hit_family'].isin(['Unknown','Unclassified'])]
cols=known['pOTU'].tolist()
#filter the list to keep only those that are >20% prevalent
cols=[c for c in cols if c in potu_corr.columns.tolist()]
known_corr=potu_corr[cols]

#Find strongest positively correlated MAG for each pOTU without any threshold
def find_corr_mags(corr_df):
    Corr_df=pd.DataFrame()
    for p in corr_df:
        mag=corr_df.loc[corr_df[p]==corr_df[p].max()].index.tolist()
        cor=pd.DataFrame({'pOTU':[p],'MAG':mag,'Corr':[corr_df[p].max()]})
        Corr_df=pd.concat([Corr_df,cor])
                          
    Corr_df=Corr_df.merge(potu_taxon[['pOTU','Hit_family']], on='pOTU',how='left')
    Corr_df=Corr_df.merge(mag_taxon[['MAG','family']], on='MAG',how='left')

    return Corr_df

known_df=find_corr_mags(known_corr)
known_df['correspondence']=known_df.apply(lambda row: 'correspond' if row.Hit_family==row.family else 'do not correspond', axis=1)
known_df.to_csv('/'.join([wdir, 'results/pOTUs_MAG_taxonomy_correlation.csv']), sep='\t', index=False)

step=0.1
thr_it=0.1

CorrRates=pd.DataFrame()
while thr_it<0.99:
    filt=known_df.query(f'Corr>={thr_it}')
    rates=pd.DataFrame(filt['correspondence'].value_counts()).T
    rates.index=[thr_it]
    CorrRates=pd.concat([CorrRates,rates], ignore_index=False)
    thr_it+=step
CorrRates=CorrRates.reset_index()
CorrRates=CorrRates.rename(columns={'index':'CorrThr'})

CorrRates['TotalN']=CorrRates['correspond']+CorrRates['do not correspond']
CorrRates['Correct%']=CorrRates['correspond']/CorrRates['TotalN']*100
CorrRates['Label']=CorrRates.apply(lambda row: f'{row.CorrThr:.2f}\n (n={row.TotalN})',axis=1)
fig=sb.lineplot(CorrRates,x='CorrThr',y='Correct%',color='#5fc0bf')
fig=sb.scatterplot(CorrRates,x='CorrThr',y='Correct%',color='#5fc0bf')
labels=['0']
labels.extend(CorrRates['Label'].tolist())
fig.set_xticklabels(labels)
fig.set(xlabel='Spearman correlation threshold',ylabel='Corresponding family delineations,%')
plt.savefig('/'.join([wdir,'results/PlasmidTaxonomyDelineation_MaxSpearmanCorrThreshold.pdf']))

##Find best matching MAGs  for pOTUs with unknown or unclassified status
unknown=potu_taxon.loc[potu_taxon['Hit_family'].isin(['Unknown','Unclassified'])]
cols=unknown['pOTU'].tolist()
#filter the list to keep only those that are >20% prevalent
cols=[c for c in cols if c in potu_corr.columns.tolist()]
unknown_corr=potu_corr[cols]

thr=0.3
Un_df=find_corr_mags(unknown_corr)
Un_df=Un_df.query(f'Corr>={thr}')
Un_df=Un_df.rename(columns={'family':'MAG_inferred_family'})

#Add cross-correlating pOTUs to the list

def add_crosspotu(corr_df):
    
    crosspotu=pd.read_csv('/'.join([wdir, 'results/Cross_correlating_pOTUs.csv']))
    #check if any of the no-familypredicted pOTUs are among cross-correlating
    crosspotu=crosspotu.loc[crosspotu['Tax1'].isin(corr_df['pOTU'].tolist())]
    
    crosspotu=crosspotu[['Tax1','Tax2']]
    crosspotu=crosspotu.merge(corr_df,left_on='Tax1',right_on='pOTU',how='left')
    crosspotu=crosspotu.drop(columns=['Tax1','pOTU'])
    crosspotu=crosspotu.rename(columns={'Tax2':'pOTU'})
    corr_df=pd.concat([corr_df,crosspotu], ignore_index=True)
    
    return corr_df
    
Un_df=add_crosspotu(Un_df)
Un_df['MAG_inferred_family'].value_counts()
Un_df.to_csv('/'.join([wdir, 'results/pOTUs_MAG_taxonomy_inferance_unknown.csv']), sep='\t', index=False)

#Compare class predictions to unknown pOTUs
unknown_df=Un_df.loc[Un_df['Hit_family']=='Unknown']
#find class prediction for the unknown hits
unknown_df['pred_class']=unknown_df['host_taxonomy'].apply(lambda row: row.split('c__')[1] if 'c__' in row else 'no class')
unknown_df['pred_class']=unknown_df['pred_class'].apply(lambda row: row.split(';')[0] if ';' in row else row)

unknown_df['pred_order']=unknown_df['host_taxonomy'].apply(lambda row: row.split('o__')[1] if 'o__' in row else 'no order')
unknown_df['pred_order']=unknown_df['pred_order'].apply(lambda row: row.split(';')[0] if ';' in row else row)
unknown_df['class_correspond']=unknown_df.apply(lambda row: 'Yes' if row['pred_class']==row['class'] else 'No', axis=1)


#Combine both known and unknown taxonomies
known_df=pd.read_csv('/'.join([wdir, 'results/pOTUs_MAG_taxonomy_correlation.csv']), sep='\t')
known_df=known_df.rename(columns={'family':'MAG_inferred_family'})
known_df=add_crosspotu(known_df)

#Add pOTUs with no reference in IMG/PR
potu_corr=pd.read_csv('/'.join([wdir, 'results/Relab_Spearman_correlation_MAGs_pOTUs_square.csv']),sep=',')
potu_corr=potu_corr.rename(columns={'Unnamed: 0':'MAG'})
potu_corr=potu_corr.set_index('MAG')

cols=[c for c in potu_corr.columns.tolist() if c not in potu_taxon['pOTU'].unique().tolist()]
noref_df=find_corr_mags(potu_corr[cols]) #7754 pOTUs
noref_df=noref_df.query(f'Corr>={thr}')
noref_df=add_crosspotu(noref_df)
noref_df=noref_df.rename(columns={'family':'MAG_inferred_family'})

dfcols=['pOTU','Hit_family','MAG_inferred_family','Corr']
inftax=pd.concat([known_df[dfcols],Un_df[dfcols],noref_df[dfcols]],ignore_index=True)
inftax.to_csv('/'.join([wdir, 'results/pOTUs_MAG_taxonomy_inferance_all.csv']), sep='\t', index=False)


