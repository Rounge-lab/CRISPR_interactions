#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  4 11:23:42 2025

Make a summary for plasmids detected in CRCbiome

@author: ekateria
"""

import pandas as pd
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
import seaborn as sb
import matplotlib.pyplot as plt

wdir='PATH_TO_MANUS_FOLDER' #set working directory

potus=pd.read_csv('/'.join([wdir,'datasets/pOTUs_lengths.csv']))

potu_taxon=pd.read_csv('/'.join([wdir, 'datasets/pOTU_taxonomy_IMGPR_0.9.txt']),sep=',')

potu_taxon['host_taxonomy']=potu_taxon['host_taxonomy'].fillna('Unclassified')
potu_taxon['Hit_family']=potu_taxon['host_taxonomy'].apply(lambda row: row.split('f__')[1] if 'f__' in row else (row if row=='Unclassified' else 'Unknown'))
potu_taxon['Hit_family']=potu_taxon['Hit_family'].apply(lambda row: row.split(';')[0] if ';' in row else row)
potu_taxon=potu_taxon.rename(columns={'Query':'pOTU'})
potu_taxon=potu_taxon.query('LongestMatchPairwiseId>=90') #keep only hits with>90% identity

wholedb=pd.read_csv('PATH_TO_DATABASES/IMGPR/IMGPR_plasmid_data.tsv', sep='\t')

#add mag taxonomy
mag_taxon=pd.read_csv('/'.join([wdir, 'datasets/MAG_taxonomy_full.tsv']), sep=',')

#Rename pOTUs (add it to pOTU_taxon and save renaming)
for ix, p in potus.iterrows():
    potus.at[ix,'pOTU_new']=f'CRCbiome-pOTU_{ix+1:05}'

potu_key=potus[['pOTU','pOTU_new']]
potu_key.to_csv('/'.join([wdir,'datasets/pOTU_rename_key.csv']),sep='\t', index=False)

def modify_conjugation(df):
    
    df['origin_of_transfer']=df['origin_of_transfer'].fillna('Not found')
    df['origin_of_transfer']=df['origin_of_transfer'].apply(lambda row: row if 'Not found' in row else 'Yes')
    
    conj_cols=['mob_genes','t4cp_genes','t4ss_atpase_genes','other_conjugation_genes']
    for c in conj_cols:
        df[c]=df[c].fillna('Not found')
        df[c]=df[c].apply(lambda row: row if 'Not found' in row else 'Yes')
    
    df['conjugation']=df.apply(lambda row: 'Yes' if 'Yes' in [row.mob_genes, row.t4cp_genes, row.t4ss_atpase_genes, row.other_conjugation_genes]
                               else 'Not found', axis=1)
    
    return df
    
potu_taxon=modify_conjugation(potu_taxon)
potu_taxon.to_csv('/'.join([wdir, 'datasets/pOTU_taxonomy_IMGPR_0.9_modified.txt']), index=False)
wholedb=modify_conjugation(wholedb)
    
potus=potus.merge(potu_taxon[['pOTU', 'source_type','ecosystem','Hit_topology','putatively_complete','putative_phage_plasmid',
                              'conjugation','origin_of_transfer']],
                  on='pOTU',how='left')

def split_ecos(df, col):
    assoc=df[col].str.split(';', expand=True)
    assoc.columns=[f'L{str(i+1)}' for i in range(assoc.shape[1])]
    df=pd.concat([df, assoc], axis=1)
    df=df.drop(columns=col)
    
    return df

potus=split_ecos(potus,'ecosystem')
wholedb=split_ecos(wholedb,'ecosystem')

potus.to_csv('/'.join([wdir,'datasets/pOTUs_db_summary.csv']), index=False)
potus=pd.read_csv('/'.join([wdir,'datasets/pOTUs_db_summary.csv']))

#Exclude those that are not in the manus dataset

relab=pd.read_csv('/'.join([wdir, 'datasets/pOTUs_relab.tsv']),sep='\t')
relab=relab.set_index('sample_id')

sums=relab.sum(axis=0).reset_index()
sums.columns=['pOTU','sum']
torem=sums.query('sum==0') #all plasmids are mapping to the samples

#Make a descriptive summary of our databases

def make_summary(cols, potus=potus, wholedb=wholedb):
    
    Summary=pd.DataFrame()
    for c in cols:
        csum=potus[c].value_counts()
        csum=pd.DataFrame(csum).reset_index().rename(columns={'index':'Value',c:'Count'})
        
        if 'Hit_' in c:
            wsum=wholedb[c.replace('Hit_','')].value_counts()
        else:
            wsum=wholedb[c].value_counts()
        wsum=pd.DataFrame(wsum)
        wsum.columns=['Count']

        csum['Column']=c
        csum['CRCbiomedb_prop']=None
        csum['Wholedb_prop']=None
        csum['Binomial_p_greater']=None
        csum['Binomial_p_twosided']=None
        
        for ix, row in csum.iterrows():
            prop=wsum.at[row.Value,'Count']/wsum['Count'].sum()
            csum.at[ix,'Wholedb_prop']=prop
            
            csum.at[ix,'CRCbiomedb_prop']=row.Count/csum['Count'].sum()

            gp=stats.binomtest(row.Count, csum['Count'].sum(), p=prop, alternative='greater')
            csum.at[ix, 'Binomial_p_greater']=gp.pvalue
            
            tsp=stats.binomtest(row.Count, csum['Count'].sum(), p=prop, alternative='two-sided')
            csum.at[ix, 'Binomial_p_twosided']=tsp.pvalue
        
        Summary=pd.concat([Summary, csum])
        
    return Summary

cols=potus.columns.tolist()[2:]  
Summary=make_summary(cols)

Summary=Summary[['Column','Value','Count','CRCbiomedb_prop','Wholedb_prop', 'Binomial_p_greater','Binomial_p_twosided']]
Summary.to_csv('/'.join([wdir,'datasets/pOTUs_db_summary_Counts.csv']), index=False)

#Plot the values

Summary=pd.read_csv('/'.join([wdir,'datasets/pOTUs_db_summary_Counts.csv']))

#Prepare for plotting

forplot=pd.DataFrame()
for c in ['L1','L2','L3']:
    df=Summary.loc[Summary['Column']==c]
    if c!='L1':
        if c=='L2':
            lim=100
        else:
            lim=50
        
        df['Value']=df.apply(lambda row: row.Value if row.Count>lim else 'Other', axis=1)
        cols=df.columns.tolist()[2:]
        df = df.groupby(['Column', 'Value'], as_index=False)[cols].sum()
        df['isother']=df['Value'].apply(lambda row: True if row=='Other' else False)
        df=df.sort_values(by='isother', ascending=True).drop(columns=['isother'])
    forplot=pd.concat([forplot,df])

forplot = forplot[['Column','Value','CRCbiomedb_prop','Wholedb_prop']].melt(id_vars=['Column','Value'], value_vars=['CRCbiomedb_prop','Wholedb_prop'], 
                    var_name='Database', value_name='Proportion')

forplot['Proportion']=forplot['Proportion']*100

##Plot only our database and a star for the IMG/PR database
fp=forplot.loc[forplot['Column']=='L2']
sb.barplot(data=fp.loc[fp['Database']=='CRCbiomedb_prop'], y='Value', x='Proportion', 
           color='#CEBA7D', width=0.7, native_scale=True, dodge=True, legend=False)
plt.ylabel('Source')
plt.xlabel('Fraction of database records,%')

sb.scatterplot(data=fp.loc[fp['Database']=='Wholedb_prop'], y='Value', x='Proportion', 
           color='k', edgecolor='k', marker='*', s=200,legend=False)

plt.savefig('/'.join([wdir,'pOTUs_withRef_source_proportions_ONLY_L2.pdf']))  

##Find if there is a correspondence between Spearman correlation predictions and taxonomy delineation

potu_taxon=pd.read_csv('/'.join([wdir, 'datasets/pOTU_taxonomy_IMGPR_0.9_modified.txt']),sep=',')
potu_corr=pd.read_csv('/'.join([wdir, 'results/Relab_Spearman_correlation_MAGs_pOTUs_square.csv']),sep=',')
potu_corr=potu_corr.rename(columns={'Unnamed: 0':'MAG'})
potu_corr=potu_corr.set_index('MAG')
#Keep only potus with known family delineation and those that are in baseline samples
potu_taxon=potu_taxon.loc[potu_taxon['pOTU'].isin(potus['pOTU'].tolist())]
potu_taxon=potu_taxon.loc[~potu_taxon['Hit_family'].isin(['Unknown','Unclassified'])]
cols=potu_taxon['pOTU'].tolist()
#filter the list to keep only those that are >20% prevalent
cols=[c for c in cols if c in potu_corr.columns.tolist()]
potu_corr=potu_corr[cols]

#Find strongest positively correlated MAG for each pOTU
Corr_df=pd.DataFrame()
for p in potu_corr:
    mag=potu_corr.loc[potu_corr[p]==potu_corr[p].max()].index.tolist()
    cor=pd.DataFrame({'pOTU':[p],'MAG':mag,'Corr':[potu_corr[p].max()]})
    Corr_df=pd.concat([Corr_df,cor])
                      
Corr_df=Corr_df.merge(potu_taxon[['pOTU','Hit_family']], on='pOTU',how='left')
Corr_df=Corr_df.merge(mag_taxon[['MAG','family']], on='MAG',how='left')

Corr_df['correspondence']=Corr_df.apply(lambda row: 'correspond' if row.Hit_family==row.family else 'do not correspond', axis=1)

#Make labels
num_repr=pd.DataFrame(Corr_df['Hit_family'].value_counts()).reset_index()
num_repr=num_repr.rename(columns={'index':'Hit_family','Hit_family':'NumRepr'})
num_repr['Label']=num_repr.apply(lambda row: row['Hit_family'] if row.NumRepr>=10 else 'Other', axis=1)

#Make labels for the plot
other=num_repr.loc[num_repr['Label']=='Other','NumRepr'].sum()

num_repr['NumRepr']=num_repr['NumRepr'].apply(lambda row: row if row>=10 else other)

order=num_repr.loc[num_repr['Label']!='Other']
order['Label']=order.apply(lambda row: f'{row.Label} (n={row.NumRepr})', axis=1)
order=order['Label'].tolist()

order.append(f'Other (n={other})')

#Add number of representatives to the summary
Corr_df=Corr_df.merge(num_repr[['Hit_family','Label','NumRepr']], on='Hit_family', how='left')
Corr_df['Label']=Corr_df.apply(lambda row: f'{row.Label} (n={row.NumRepr})', axis=1)
Corr_df.to_csv('/'.join([wdir, 'results/Family_delineation_correspondence_pOTUs_MAGs_SpCor.csv']),index=False)

fig=sb.boxplot(Corr_df, x='Label', y='Corr',hue='correspondence', hue_order=['correspond', 'do not correspond'],
              palette=['#ace1a5', '#eac0bf'], order=order)
plt.xticks(rotation=90, fontstyle='italic')
fig.set(xlabel='',ylabel='Highest Spearman correlation', ylim=[0, 1])
plt.savefig('/'.join([wdir, 'results/Family_delineation_correspondence_pOTUs_MAGs_SpCor.pdf']))

#Make a crosstab between family and correspondence
crosstab=pd.crosstab(Corr_df['Label'],Corr_df['correspondence'])
crosstab['prop_yes']=crosstab['correspond']/(crosstab['correspond']+crosstab['do not correspond'])*100
crosstab['prop_no']=crosstab['do not correspond']/(crosstab['correspond']+crosstab['do not correspond'])*100

crosstab=crosstab.sort_values(by='prop_yes',ascending=True)
ax=crosstab[['prop_yes','prop_no']].plot(kind='barh', stacked=True, color=['#5895a8','#a4a8a8'], legend=False)
plt.yticks(fontstyle= 'italic')
plt.xlabel('Proportion pOTUs, %')
plt.ylabel('')
plt.savefig('/'.join([wdir,'results/Family_delineation_correspondence_pOTUs_MAGs.pdf']))

#Check if there is higher correlation between those that correspond than those that do not correspond
def kruskal_group(data,cat_col,y):
    query=data.copy()
    query=query.dropna(subset=[cat_col,y])
    categ=query[cat_col].unique().tolist()
    for_stats={}
    for cat in categ:
        for_stats[cat]=query.loc[query[cat_col]==cat,y].tolist()

    h_st,pval=stats.kruskal(*list(for_stats.values())) # *indicates to treat each element in list as a group, the list will contain as many elements as there are dict entries, each entry is a separate list within a list

    return h_st,pval

h_st, pval=kruskal_group(Corr_df,'correspondence','Corr') 

KruskalCorr=pd.DataFrame()
for f in order:
    data=Corr_df.loc[Corr_df['Label']==f]
    if len(data['correspondence'].unique().tolist())>1:
        h_st,p=kruskal_group(data,'correspondence','Corr')
        kr=pd.DataFrame({'Label':[f],'Hst':[h_st],'Pval':[p]})
        KruskalCorr=pd.concat([KruskalCorr,kr])
_, padj, _, _ = multipletests(KruskalCorr['Pval'], method='fdr_bh')
KruskalCorr['FDRp']=padj
KruskalCorr.to_csv('/'.join([wdir, 'results/Family_delineation_correspondence_pOTUs_MAGs_SpCor_KruskalWallis.csv']),index=False)

##------------------------------------------
#for each host family (most prevalent pOTU in each family) find closely correlated MAGs and their delineation and plot

potu_summary=pd.read_csv('/'.join([wdir,'results/MAGs_vOTUs_pOTUs/pOTU_prevalence_query_hit_lengthratio.csv']), sep='\t')

#Get correlation matrix
potu_corr=pd.read_csv('/'.join([wdir, 'results/Relab_Spearman_correlation_MAGs_pOTUs_square.csv']),sep=',')
potu_corr=potu_corr.rename(columns={'Unnamed: 0':'MAG'})
potu_corr=potu_corr.set_index('MAG')

#Find most prevalent pOTUs with each host delineation

taxfam=potu_summary['Hit_family'].unique().tolist()
taxfam=[t for t in taxfam if t not in ['No reference', 'Unclassified','Unknown','Other']]
taxfam=[t for t in taxfam if 'CAG' not in t]
taxfam=[t for t in taxfam if 'U' not in t]


Most_prev=pd.DataFrame()
for fam in taxfam:
    qfam=potu_summary.loc[potu_summary['Hit_family']==fam]
    mostprev=qfam.loc[qfam['Prevalence, %']==max(qfam['Prevalence, %'])]
    Most_prev=pd.concat([Most_prev,mostprev])
    
Most_prev=Most_prev.drop_duplicates(subset='Hit_family')
    
Most_prev=Most_prev.loc[Most_prev['Prevalence, %']>=20]

plist=Most_prev['pOTU'].tolist()
corr_toplot=potu_corr[plist]

mags_toplot=[]
for p in plist:
    #find strongly correlating mags
    mag=corr_toplot.loc[corr_toplot[p]==max(corr_toplot[p])].index.tolist()
    mags_toplot.extend(mag)
    
corr_toplot=corr_toplot.loc[mags_toplot]

fig=sb.heatmap(corr_toplot.T,cmap='vlag', center=0)
fig.set_xticklabels(fig.get_xticklabels(),fontstyle='italic')
fig.set_yticklabels(fig.get_yticklabels(),fontstyle='italic')
fig.set_xlabel('')
fig.set_ylabel('')

plt.savefig('/'.join([wdir,'results/MostPrevalent_pOTUs_to_MAGs_SpearmanCorr.pdf']))


#Find strongest correlation to unknown hosts
potu_taxon=pd.read_csv('/'.join([wdir, 'datasets/pOTU_taxonomy_IMGPR_0.9_modified.txt']),sep=',')
potu_corr=pd.read_csv('/'.join([wdir, 'results/Relab_Spearman_correlation_MAGs_pOTUs_square.csv']),sep=',')
potu_corr=potu_corr.rename(columns={'Unnamed: 0':'MAG'})
potu_corr=potu_corr.set_index('MAG')
#Keep only potus with known family delineation and those that are in baseline samples
potu_taxon=potu_taxon.loc[potu_taxon['pOTU'].isin(potus['pOTU'].tolist())]
potu_taxon=potu_taxon.loc[potu_taxon['Hit_family'].isin(['Unknown','Unclassified'])]
cols=potu_taxon['pOTU'].tolist()
#filter the list to keep only those that are >20% prevalent
cols=[c for c in cols if c in potu_corr.columns.tolist()]
potu_corr=potu_corr[cols]

#Find strongest positively correlated MAG for each pOTU
Corr_df=pd.DataFrame()
for p in potu_corr:
    mag=potu_corr.loc[potu_corr[p]==potu_corr[p].max()].index.tolist()
    cor=pd.DataFrame({'pOTU':[p],'MAG':mag,'Corr':[potu_corr[p].max()]})
    Corr_df=pd.concat([Corr_df,cor])
                      
Corr_df=Corr_df.merge(potu_taxon[['pOTU','Hit_family']], on='pOTU',how='left')
Corr_df=Corr_df.merge(mag_taxon[['MAG','family']], on='MAG',how='left')

sb.boxplot(Corr_df, y='family',x='Corr')