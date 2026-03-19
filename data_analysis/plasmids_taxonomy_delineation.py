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

PTU_taxon=pd.read_csv('/'.join([wdir, 'datasets/PTU_taxonomy_IMGPR_0.9_modified.txt']),sep=',')
mag_taxon=pd.read_csv('/'.join([wdir, 'datasets/MAG_taxonomy_full.tsv']), sep=',')
PTU_corr=pd.read_csv('/'.join([wdir, 'results/Relab_Spearman_correlation_MAGs_PTUs_square.csv']),sep=',')
PTUs=pd.read_csv('/'.join([wdir,'datasets/PTUs_db_summary.csv']))

PTU_corr=PTU_corr.rename(columns={'Unnamed: 0':'MAG'})
PTU_corr=PTU_corr.set_index('MAG')

#Keep only PTUs with known family delineation and those that are in baseline samples
PTU_taxon=PTU_taxon.loc[PTU_taxon['PTU'].isin(PTUs['PTU'].tolist())]
known=PTU_taxon.loc[~PTU_taxon['Hit_family'].isin(['Unknown','Unclassified'])]
cols=known['PTU'].tolist()
#filter the list to keep only those that are >20% prevalent
cols=[c for c in cols if c in PTU_corr.columns.tolist()]
known_corr=PTU_corr[cols]

#Find strongest positively correlated MAG for each PTU without any threshold
def find_corr_mags(corr_df):
    Corr_df=pd.DataFrame()
    for p in corr_df:
        mag=corr_df.loc[corr_df[p]==corr_df[p].max()].index.tolist()
        cor=pd.DataFrame({'PTU':[p],'MAG':mag,'Corr':[corr_df[p].max()]})
        Corr_df=pd.concat([Corr_df,cor])
                          
    Corr_df=Corr_df.merge(PTU_taxon[['PTU','Hit_family']], on='PTU',how='left')
    Corr_df=Corr_df.merge(mag_taxon[['MAG','family']], on='MAG',how='left')

    return Corr_df

known_df=find_corr_mags(known_corr)
known_df['correspondence']=known_df.apply(lambda row: 'correspond' if row.Hit_family==row.family else 'do not correspond', axis=1)
known_df.to_csv('/'.join([wdir, 'results/PTUs_MAG_taxonomy_correlation.csv']), sep='\t', index=False)

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

##Find best matching MAGs  for PTUs with unknown or unclassified status
unknown=PTU_taxon.loc[PTU_taxon['Hit_family'].isin(['Unknown','Unclassified'])]
cols=unknown['PTU'].tolist()
#filter the list to keep only those that are >20% prevalent
cols=[c for c in cols if c in PTU_corr.columns.tolist()]
unknown_corr=PTU_corr[cols]

thr=0.3
Un_df=find_corr_mags(unknown_corr)
Un_df=Un_df.query(f'Corr>={thr}')
Un_df=Un_df.rename(columns={'family':'MAG_inferred_family'})

#Add cross-correlating PTUs to the list

def add_crossPTU(corr_df):
    
    crossPTU=pd.read_csv('/'.join([wdir, 'results/Cross_correlating_PTUs.csv']))
    #check if any of the no-familypredicted PTUs are among cross-correlating
    crossPTU=crossPTU.loc[crossPTU['Tax1'].isin(corr_df['PTU'].tolist())]
    
    crossPTU=crossPTU[['Tax1','Tax2']]
    crossPTU=crossPTU.merge(corr_df,left_on='Tax1',right_on='PTU',how='left')
    crossPTU=crossPTU.drop(columns=['Tax1','PTU'])
    crossPTU=crossPTU.rename(columns={'Tax2':'PTU'})
    corr_df=pd.concat([corr_df,crossPTU], ignore_index=True)
    
    return corr_df
    
Un_df=add_crossPTU(Un_df)
Un_df['MAG_inferred_family'].value_counts()
Un_df.to_csv('/'.join([wdir, 'results/PTUs_MAG_taxonomy_inferance_unknown.csv']), sep='\t', index=False)

#Compare class predictions to unknown PTUs
unknown_df=Un_df.loc[Un_df['Hit_family']=='Unknown']
#find class prediction for the unknown hits
unknown_df['pred_class']=unknown_df['host_taxonomy'].apply(lambda row: row.split('c__')[1] if 'c__' in row else 'no class')
unknown_df['pred_class']=unknown_df['pred_class'].apply(lambda row: row.split(';')[0] if ';' in row else row)

unknown_df['pred_order']=unknown_df['host_taxonomy'].apply(lambda row: row.split('o__')[1] if 'o__' in row else 'no order')
unknown_df['pred_order']=unknown_df['pred_order'].apply(lambda row: row.split(';')[0] if ';' in row else row)
unknown_df['class_correspond']=unknown_df.apply(lambda row: 'Yes' if row['pred_class']==row['class'] else 'No', axis=1)


#Combine both known and unknown taxonomies
known_df=pd.read_csv('/'.join([wdir, 'results/PTUs_MAG_taxonomy_correlation.csv']), sep='\t')
known_df=known_df.rename(columns={'family':'MAG_inferred_family'})
known_df=add_crossPTU(known_df)

#Add PTUs with no reference in IMG/PR
PTU_corr=pd.read_csv('/'.join([wdir, 'results/Relab_Spearman_correlation_MAGs_PTUs_square.csv']),sep=',')
PTU_corr=PTU_corr.rename(columns={'Unnamed: 0':'MAG'})
PTU_corr=PTU_corr.set_index('MAG')

cols=[c for c in PTU_corr.columns.tolist() if c not in PTU_taxon['PTU'].unique().tolist()]
noref_df=find_corr_mags(PTU_corr[cols]) #7754 PTUs
noref_df=noref_df.query(f'Corr>={thr}')
noref_df=add_crossPTU(noref_df)
noref_df=noref_df.rename(columns={'family':'MAG_inferred_family'})

dfcols=['PTU','Hit_family','MAG_inferred_family','Corr']
inftax=pd.concat([known_df[dfcols],Un_df[dfcols],noref_df[dfcols]],ignore_index=True)
inftax.to_csv('/'.join([wdir, 'results/PTUs_MAG_taxonomy_inferance_all.csv']), sep='\t', index=False)


