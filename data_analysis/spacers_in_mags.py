#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 25 12:10:33 2025

Make a summary of which MAG families contain CRISPRs

@author: ekateria
"""

import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt

wdir='PATH_TO_MANUS_FOLDER'
nfcassettes=pd.read_csv('/'.join([wdir, 'datasets/Combined_crisprs_all_NOT_filtered.tab']), sep = '\t')

mags_contigs=pd.read_csv('/'.join([wdir, 'datasets/MAGs_contig_list_OneByOne.csv']), sep=',')
#keep only mags detected in the dataset
magsmanus=pd.read_csv('/'.join([wdir, 'datasets/MAGs_relab.tsv']), sep='\t').set_index('sample_id')
magsmanus=magsmanus.columns.tolist()
mags_contigs=mags_contigs.loc[mags_contigs['MAG'].isin(magsmanus)]

nfcassettes=nfcassettes.merge(mags_contigs, on='Contig', how='left')

cas_mags=nfcassettes.dropna(subset='MAG')['MAG'].unique().tolist()

taxonomy=pd.read_csv('/'.join([wdir,'datasets/MAG_taxonomy_full.tsv']))
mags=pd.read_csv('/'.join([wdir,'datasets/MAGs_relab.tsv']),sep='\t')
maglist=mags.columns.tolist()[1:]
taxonomy=taxonomy.loc[taxonomy['MAG'].isin(maglist)]

taxonomy['CRISPRDet']=taxonomy['MAG'].apply(lambda row: 'Yes' if row in cas_mags else 'No')

taxlev='family'
sum_tax=taxonomy[[taxlev,'CRISPRDet']].pivot_table(index=taxlev, columns='CRISPRDet', aggfunc='size',fill_value=0)
sum_tax['Total']=sum_tax['Yes']+sum_tax['No']
sum_tax=sum_tax.sort_values(by='Yes', ascending=False)

sum_tax['PropYes']=sum_tax['Yes']/sum_tax['Total']*100
sum_tax['PropNo']=sum_tax['No']/sum_tax['Total']*100

sum_tax[['Total','PropYes']].to_csv('/'.join([wdir, 'results/MAG_families_CRISPRprop.csv']), sep='\t', index=True)

sum_tax['label']=sum_tax['Total'].apply(lambda row: f'(n={str(row)})')
sum_tax=sum_tax.reset_index()
sum_tax[taxlev]=sum_tax.apply(lambda row: f'{row[taxlev]} {row.label}', axis=1)

nocrispr=sum_tax.query('PropYes==0')
sum_tax=sum_tax.set_index(taxlev)
sum_tax=sum_tax.query('Total>=10') #Plot families that have more than 10 representatives
sum_tax=sum_tax.sort_values(by='PropYes',ascending=True)
ax=sum_tax[['PropYes','PropNo']].plot(kind='barh', stacked=True, color=['#5fc0bf','#aab4c2'], legend=False)
plt.yticks(fontstyle= 'italic')
plt.savefig('/'.join([wdir,f'results/CRISPR_in_MAG{taxlev}_prop.pdf']))

sum_tax.to_csv('/'.join([wdir,f'results/CRISPR_in_MAG{taxlev}_prop.csv']), index=True)

##Check if there is any difference in repeats length with regards to genome lengths and taxonomy
mag_length=pd.read_csv('/'.join([wdir, 'datasets/MAGs_lengths.csv']))
mag_cassettes=nfcassettes.dropna(subset='MAG')
mag_cassettes=mag_cassettes.merge(mag_length, on='MAG', how='left')
mag_cassettes=mag_cassettes.merge(taxonomy[['MAG','class','order','family', 'genus', 'species']], on='MAG', how='left')

sb.regplot(mag_cassettes, x='Length, bp', y='N_repeats')
pr, p = stats.pearsonr(mag_cassettes['Length, bp'],mag_cassettes['N_repeats'])

#Plot number of repeats with regards to taxonomic family

#make an order by avg number of repeats
order = mag_cassettes.groupby('family')['N_repeats'].mean().sort_values(ascending=False).index

#make a label
labels=pd.DataFrame()
for fam in order:
    m=mag_cassettes.loc[mag_cassettes['family']==fam]
    nmags=len(m['MAG'].unique().tolist())
    nsams=len(m['sample_id'].unique().tolist())
    
    if nmags>2 and nsams>=5:
        toplot='yes'
    else:
        toplot=None
    
    lab=pd.DataFrame({'family':[fam],'label':[f'{fam} (n={nmags},s={nsams})'], 'plot':[toplot]})
    labels=pd.concat([labels, lab])
    
labels=labels.dropna(subset='plot')
    
fig=sb.boxplot(mag_cassettes, x='N_repeats', y='family', color='#4385BF', order=labels['family'].tolist())
fig.set_yticklabels(labels['label'].tolist(),fontsize=8, fontstyle='italic')
fig.set_xlabel('Number of repeats')
fig.set_ylabel('')
plt.savefig('/'.join([wdir,'results/NumRepeats_inMAGs_by_family.pdf']))

#Get the list of families with most and least number of repeats
repcount=mag_cassettes.groupby('family')['N_repeats'].std().sort_values(ascending=False).reset_index(name='count')
repcount=repcount.loc[repcount['family'].isin(labels['family'].tolist())]

#Number of cassettes in MAGs and their avg length
mag_cassettes['sample_id']=mag_cassettes['Contig'].apply(lambda row: row.split('_')[0])
counts=mag_cassettes.groupby(['sample_id','MAG']).size().reset_index(name='count')
counts=counts.merge(taxonomy[['MAG','class','order','family','genus','species']], on='MAG', how='left')
order = counts.groupby('family')['count'].mean().sort_values().index
fig=sb.boxplot(counts, x='count', y='family', color='#4385BF', order=order)

unique_mags=taxonomy.loc[taxonomy['MAG'].isin(counts['MAG'].unique().tolist())]

