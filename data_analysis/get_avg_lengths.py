#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 24 15:41:08 2025

Calculate average lengths for MAGs, vOTUs
Make summary plots for taxonomy, genome lengths and completeness

@author: ekateria
"""

import pandas as pd
from Bio import SeqIO
import seaborn as sb
import matplotlib.pyplot as plt
import matplotlib.patches as patches

import pyreadr as pyr
import numpy as np
from scipy.stats import pearsonr


wdir='PATH_TO_MANUS_FOLDER' #set working directory

prefix='MAGs_vOTUs_pOTUs'

## Get the lengths of all contigs for general dataset stats
contigs=pd.read_csv('/'.join([wdir, 'datasets/All_contig_lengths.csv']),sep='\t')
contigs[['sample_id','num']]=contigs['Contig'].str.split('_',expand=True)

samples=pyr.read_r('/'.join([wdir, 'participant_data/sample_meta.Rds']))[None]
baseline=samples['sample_id'].tolist()

baseline=[b.replace('_','-') for b in baseline]
contigs=contigs.loc[contigs['sample_id'].isin(baseline)]
contigstats=contigs['Length'].describe()

##Make a summary for MAGs and vOTUs
mags=pd.read_csv('/'.join([wdir, 'datasets/MAGs_relab.tsv']),sep='\t')
maglist=mags.columns.tolist()[1:]

magids=pd.read_csv('PATH_TO_MAGS/genomes/clustering/old2newID.tsv', sep='\t')
magids=magids.loc[magids['MAG'].isin(maglist)]

#Read fasta file and get the sum of all sequence lengths
magsloc='PATH_TO_MAGS/genomes/Dereplication/dereplicated_genomes'

def sum_length(fasta):
    tot_len=0
    for seq in SeqIO.parse(fasta,'fasta'):
        tot_len+=len(seq.seq)
    return tot_len

mag_lengths={'MAG':[],'Length, bp':[]}

for ix,mag in magids.iterrows():
    maglen=sum_length(f'{magsloc}/{mag.BinID}.fasta')
    mag_lengths['MAG'].append(mag.MAG)
    mag_lengths['Length, bp'].append(maglen)

mag_lengths=pd.DataFrame(mag_lengths)

mag_lengths.to_csv('/'.join([wdir,'datasets/MAGs_lengths.csv']), index=False)

#Extract vOTU lengths
votus=pd.read_csv('/'.join([wdir, 'datasets/vOTUs_relab.tsv']),sep='\t')
votulist=votus.columns.tolist()[1:]

votu_fasta='PATH_TO_VIRUSES/dereplication/repr_viral_seqs.fasta'
votu_lengths={'vOTU':[],'Length, bp':[]}
for seq in SeqIO.parse(votu_fasta, 'fasta'):
    votu_lengths['vOTU'].append(seq.id)
    votu_lengths['Length, bp'].append(len(seq.seq))
    
votu_lengths=pd.DataFrame(votu_lengths)
votu_lengths=votu_lengths.loc[votu_lengths['vOTU'].isin(votulist)]
votu_lengths.to_csv('/'.join([wdir,'datasets/vOTUs_lengths.csv']), index=False)

votu_lengths['Length, bp'].describe()

#Extract pOTU lengths
potus=pd.read_csv('/'.join([wdir, 'datasets/pOTUs_relab.tsv']),sep='\t')
potulist=potus.columns.tolist()[1:]
potu_fasta='PATH_TO/scapp_snakemake_wf/data/dereplication/plasmids_dereplicated_0.9.fasta'
potu_lengths={'pOTU':[],'Length, bp':[]}
for seq in SeqIO.parse(potu_fasta, 'fasta'):
    potu_lengths['pOTU'].append(seq.id)
    potu_lengths['Length, bp'].append(len(seq.seq))
    
potu_lengths=pd.DataFrame(potu_lengths)
potu_lengths=potu_lengths.loc[potu_lengths['pOTU'].isin(potulist)]
potu_lengths.to_csv('/'.join([wdir,'datasets/pOTUs_lengths.csv']), index=False)
potu_lengths['Length, bp'].describe()

##Find how many genomes are included into each MAG, vOTU and pOTU
mag_numgen=pd.read_csv('/'.join([wdir,'datasets/MAGs_NumGenomes.csv']), sep='\t')
votu_numgen=pd.read_csv('/'.join([wdir,'datasets/vOTUs_NumGenomes.csv']), sep='\t')
potu_numgen=pd.read_csv('/'.join([wdir,'datasets/Plasmid_NumGenomes.csv']),sep='\t')

##------------------------------------------------------------------------
##Make plots for prevalence, lengths and completeness
mag_lengths=pd.read_csv('/'.join([wdir,'datasets/MAGs_lengths.csv']))
votu_lengths=pd.read_csv('/'.join([wdir,'datasets/vOTUs_lengths.csv']))

mag_taxon=pd.read_csv('/'.join([wdir,'datasets/MAG_taxonomy_full.tsv']), sep=',')
vir_taxon=pd.read_csv('/'.join([wdir,'datasets/vOTUs_taxonomy.csv']),sep=',')
vir_taxon=vir_taxon.drop(columns='Unnamed: 0')
vir_taxon=vir_taxon.rename(columns={'Scaffold':'vOTU'})

mag_comp=pd.read_csv('PATH_TO_MAGS/genomes/checkm/completeness.tsv', sep='\t')
mag_comp=mag_comp[['Bin Id','Completeness','# genomes']].rename(columns={'Bin Id':'MAG', '# genomes':'genome'})

vir_comp=pd.read_csv('PATH_TO_VIRUSES/dereplication/checkV_summary.tsv', sep='\t')
vir_rename=pd.read_csv('PATH_TO_VIRUSES/dereplication/old_to_new_ids.tsv', sep='\t')

vir_comp=vir_comp.merge(vir_rename, left_on='virus_id',right_on='old_id', how='left')
vir_comp=vir_comp[['new_id','completeness']].rename(columns={'new_id':'vOTU','completeness':'Completeness'})

mags=mags.set_index('sample_id')
votus=votus.set_index('sample_id')

colors={'MAG':'#4385BF', 'vOTU':'#5B7571', 'pOTU':'#CEBA7D'}

def plot_stats(relab,length,comp,taxon,numgen,groupname,taxlev,repr_thr,genome_lims,colors=colors):
    prev=relab.copy()
    if groupname=='vOTU':
        prev[prev>0]=1 #keep all results over 0 since each virus is covered on at least 75% of its length
    else:
        prev[prev>=0.001]=1 #set the threshold to 0.001 % (it makes 10000 bp out of a minimum threshold of 1Gb data per sample)
        prev[prev<0.001]=0
    summary=pd.DataFrame(prev.sum(axis=0)/prev.shape[0]*100, columns=['Prevalence, %'])
    summary=summary.reset_index().rename(columns={'index':groupname})
    
    summary=summary.merge(length, on=groupname, how='left')
    summary=summary.merge(comp, on=groupname, how='left')
    if groupname!='MAG':
        summary=summary.merge(numgen,on=groupname,how='left')
    summary=summary.merge(taxon[[groupname,taxlev]],on=groupname, how='left')
    
    #Arrange plotting order by most highly represented classes 
    num_repr=pd.DataFrame(summary[taxlev].value_counts())
    num_repr=num_repr.reset_index().rename(columns={'count':'NumRepr'})
    num_repr['Label']=num_repr.apply(lambda row: row[taxlev] if row.NumRepr>=repr_thr else 'Other', axis=1)
    
    if groupname=='vOTU':
        num_repr['Label']=num_repr['Label'].apply(lambda row: row if row!='n.a.' else 'Unknown')
        unknown=num_repr.loc[num_repr['Label']=='Unknown','NumRepr'].sum()
        unclassified=num_repr.loc[num_repr['Label']=='Unclassified','NumRepr'].sum()
        
    other=num_repr.loc[num_repr['Label']=='Other','NumRepr'].sum()
    num_repr['NumRepr']=num_repr['NumRepr'].apply(lambda row: row if row>=repr_thr else other)
    
    order=num_repr.loc[num_repr['Label']!='Other']
    if groupname=='vOTU':
        order=order.loc[order['Label']!='Unknown']
        order=order.loc[order['Label']!='Unclassified']

    order['Label']=order.apply(lambda row: f'{row.Label} (n={row.NumRepr})', axis=1)
    order=order['Label'].tolist()
    if groupname=='vOTU':
        order.append(f'Unclassified (n={unclassified})')
        order.append(f'Unknown (n={unknown})')
    order.append(f'Other (n={other})')
    
    #Add number of representatives to the summary
    summary=summary.merge(num_repr[[taxlev,'Label','NumRepr']], on=taxlev, how='left')
    summary['Label']=summary.apply(lambda row: f'{row.Label} (n={row.NumRepr})', axis=1)
    
    fig,ax=plt.subplots(1,4,sharey=True, figsize=(10, 7))
    
    sb.boxplot(data=summary, x='Prevalence, %', y='Label', ax=ax[0],color=colors[groupname],order=order)
    sb.boxplot(data=summary, x='genome', y='Label', ax=ax[1], log_scale=True, color=colors[groupname],order=order)
    sb.boxplot(data=summary, x='Length, bp', y='Label', ax=ax[2], log_scale=True, color=colors[groupname],order=order)
    sb.boxplot(data=summary, x='Completeness', y='Label', ax=ax[3],color=colors[groupname],order=order)
    
    ax[0].set_ylabel('')
    ax[0].set_yticklabels(ax[0].get_yticklabels(),fontstyle='italic')
    ax[1].set_xlabel('Number of genomes')
    ax[2].set_xlim(genome_lims)
    ax[2].set_xlabel('Genome length, bp')
    ax[3].set_xlabel('Genome completeness, %')
    
    return summary

mag_summary=plot_stats(mags.loc[mags.index.isin(baseline)],mag_lengths,mag_comp, mag_taxon,mag_numgen, 'MAG', 'class', 10, [3*10**5,10**7])
plt.tight_layout()
plt.savefig('/'.join([wdir,prefix,'MAG_prevalence_length_completeness.pdf']))
mag_summary.to_csv('/'.join([wdir,prefix,'MAG_prevalence_length_completeness.csv']), sep='\t', index=False)

votu_summary=plot_stats(votus,votu_lengths,vir_comp,vir_taxon,votu_numgen,'vOTU', 'Family', 10, [10**3, 10**6])
plt.tight_layout()
plt.savefig('/'.join([wdir,prefix,'vOTU_prevalence_length_completeness.pdf']))
votu_summary.to_csv('/'.join([wdir,prefix,'Phages_prevalence_length_completeness.csv']), sep='\t', index=False)

#Make summary for each class
mag_summary=pd.read_csv('/'.join([wdir,prefix,'MAG_prevalence_length_completeness.csv']), sep='\t')

def sum_by_taxon(summary,taxlev,metric):
    
    Descr=pd.DataFrame()
    taxon=summary[taxlev].unique().tolist()
    for t in taxon:
        s=summary.loc[summary[taxlev]==t]
        sd=pd.DataFrame(s[metric].describe()).T.reset_index()
        sd[taxlev]=t
        Descr=pd.concat([Descr,sd])
    
    return Descr

Metrics=pd.DataFrame()
for m in ['Length, bp','Completeness']:
    Descr=sum_by_taxon(mag_summary, 'class', m)
    Descr['metric']=m
    Metrics=pd.concat([Metrics,Descr])

Metrics.to_csv('/'.join([wdir,prefix,'MAG_prevalence_length_completeness_metrics.csv']), sep='\t',index=False)

votu_summary=pd.read_csv('/'.join([wdir,prefix,'Phages_prevalence_length_completeness.csv']), sep='\t')
vMetrics=pd.DataFrame()
for m in ['Length, bp','Completeness']:
    Descr=sum_by_taxon(votu_summary, 'Family', m)
    Descr['metric']=m
    vMetrics=pd.concat([vMetrics,Descr])

vMetrics.to_csv('/'.join([wdir,prefix,'Phages_prevalence_length_completeness_metrics.csv']), sep='\t',index=False)

#-------------------------------------------------------------------------------------------------------------
#Make a summary plot for plasmids

colors={'MAG':'#4385BF', 'vOTU':'#5B7571', 'pOTU':'#CEBA7D'}

samples=pyr.read_r('/'.join([wdir, 'participant_data/sample_meta.Rds']))[None]
baseline=samples['sample_id'].tolist()


potus=pd.read_csv('/'.join([wdir, 'datasets/pOTUs_relab.tsv']),sep='\t')
potus=potus.set_index('sample_id')
potu_taxon=pd.read_csv('/'.join([wdir, 'datasets/pOTU_taxonomy_IMGPR_0.9.txt']),sep=',')
potu_taxon['host_taxonomy']=potu_taxon['host_taxonomy'].fillna('Unclassified')
potu_taxon['Hit_family']=potu_taxon['host_taxonomy'].apply(lambda row: row.split('f__')[1] if 'f__' in row else (row if row=='Unclassified' else 'Unknown'))
potu_taxon['Hit_family']=potu_taxon['Hit_family'].apply(lambda row: row.split(';')[0] if ';' in row else row)
potu_taxon=potu_taxon.rename(columns={'Query':'pOTU'})
potu_taxon=potu_taxon.query('LongestMatchPairwiseId>=90')
potu_taxon['Query_Hit_ratio']=potu_taxon['QueryLength']/potu_taxon['Hit_length']
potu_lengths=pd.read_csv('/'.join([wdir,'datasets/pOTUs_lengths.csv']))


potu_prev=potus.loc[potus.index.isin(baseline)].copy()
potu_prev[potu_prev>0]=1 #plasmid needs to be covered on at least 75% of its length by default
potu_summary=pd.DataFrame(potu_prev.sum(axis=0)/potu_prev.shape[0]*100, columns=['Prevalence, %'])
potu_summary=potu_summary.reset_index()
potu_summary=potu_summary.rename(columns={'index':'pOTU'})
potu_summary=potu_summary.merge(potu_numgen, on='pOTU', how='left')
potu_summary=potu_summary.merge(potu_lengths, on='pOTU', how='left')
potu_summary=potu_summary.merge(potu_taxon[['pOTU','LongestMatchPairwiseId','Query_Hit_ratio','Hit_family']],on='pOTU',how='left')
potu_summary.loc[potu_summary['Hit_family'].isna(), 'Hit_family'] = 'No reference'


### Compare the prevalence between MAGs and plasmids belonging to that family
mag_summary=pd.read_csv('/'.join([wdir,prefix,'MAG_prevalence_length_completeness.csv']), sep='\t')
mag_taxon=pd.read_csv('/'.join([wdir,'datasets/MAG_taxonomy_full.tsv']), sep=',')
mag_summary=mag_summary.merge(mag_taxon[['MAG','family']], on='MAG', how='left')

fams=pd.DataFrame(potu_summary['Hit_family'].unique().tolist())
fams.columns=['Hit_family']
for ix,f in fams.iterrows():
    ms=mag_summary.loc[mag_summary['family']==f.Hit_family]
    if len(ms)>0:
        fams.at[ix,'NumMAGs']=len(ms)
        #Calculate quartiles 
        qm=ms['Prevalence, %'].describe()
        qm=qm.loc[['min','25%','50%','75%','max']]
        p=potu_summary.loc[potu_summary['Hit_family']==f.Hit_family]
        qp=p['Prevalence, %'].describe()
        qp=qp.loc[['min','25%','50%','75%','max']]
        cor,p=pearsonr(qp.values,qm.values)
        fams.at[ix,'Corr']=cor
        fams.at[ix,'pval']=p
    else:
        fams.at[ix,'NumMAGs']=np.nan
        fams.at[ix,'Corr']=np.nan
        fams.at[ix,'pval']=np.nan

potu_summary=potu_summary.merge(fams[['Hit_family','Corr','pval']], on='Hit_family', how='left')

#Make labels for the plot

num_repr=pd.DataFrame(potu_summary['Hit_family'].value_counts()).reset_index()
num_repr=num_repr.rename(columns={'count':'NumRepr'})
num_repr['Label']=num_repr.apply(lambda row: row['Hit_family'] if row.NumRepr>=10 else 'Other', axis=1)

unclassified=num_repr.loc[num_repr['Label']=='Unclassified','NumRepr'].sum()
unknown=num_repr.loc[num_repr['Label']=='Unknown','NumRepr'].sum()
noref=num_repr.loc[num_repr['Label']=='No reference','NumRepr'].sum()
other=num_repr.loc[num_repr['Label']=='Other','NumRepr'].sum()

num_repr['NumRepr']=num_repr['NumRepr'].apply(lambda row: row if row>=10 else other)

order=num_repr.loc[num_repr['Label']!='Other']
order=order.loc[order['Label']!='Unknown']
order=order.loc[order['Label']!='No reference']
order=order.loc[order['Label']!='Unclassified']


order['Label']=order.apply(lambda row: f'{row.Label} (n={row.NumRepr})', axis=1)
order=order['Label'].tolist()

order.append(f'Other (n={other})')
order.append(f'Unclassified (n={unclassified})')
order.append(f'Unknown (n={unknown})')
order.append(f'No reference (n={noref})')


#Add number of representatives to the summary
potu_summary=potu_summary.merge(num_repr[['Hit_family','Label','NumRepr']], on='Hit_family', how='left')
potu_summary['Label']=potu_summary.apply(lambda row: f'{row.Label} (n={row.NumRepr})', axis=1)

potu_summary=pd.read_csv('/'.join([wdir,'results',prefix,'pOTU_prevalence_query_hit_lengthratio.csv']), sep='\t')

#Add info on mobility
potumob=pd.read_csv('/'.join([wdir, 'datasets/plasmids_dereplicated_0.9_MOBtyper.txt']), sep='\t')
potumob=potumob.rename(columns={'sample_id':'pOTU'})
potu_summary=potu_summary.merge(potumob[['pOTU','predicted_mobility']], on='pOTU',how='left')

mob_summary=potu_summary.groupby(['Label','predicted_mobility']).size().unstack(fill_value=0)
perc_mobility=mob_summary.div(mob_summary.sum(axis=1),axis=0)*100
perc_mobility=perc_mobility.reset_index()
perc_mobility['mob+conj']=100-perc_mobility['non-mobilizable']


fig,ax=plt.subplots(1,5,sharey=True, figsize=(10, 7))

sb.boxplot(data=potu_summary, x='Prevalence, %', y='Label', ax=ax[0],color=colors['pOTU'],order=order)
sb.boxplot(data=potu_summary, x='genome', y='Label', ax=ax[1],color=colors['pOTU'],order=order, log_scale=True)
sb.boxplot(data=potu_summary, x='Length, bp', y='Label', ax=ax[2], log_scale=True, color=colors['pOTU'],order=order)
sb.barplot(data=perc_mobility, x='mob+conj', y='Label', ax=ax[3], color=colors['pOTU'],order=order)
sb.boxplot(data=potu_summary, x='Query_Hit_ratio', y='Label', ax=ax[4],color=colors['pOTU'],order=order)

ax[0].set_ylabel('')
ax[0].set_yticklabels(ax[0].get_yticklabels(),fontstyle='italic')
ax[1].set_xlabel('Number of plasmids')
ax[2].set_xlabel('PTU length, bp')
ax[3].set_xlabel('Mobility, %')
ax[4].set_xlabel('PTU to reference plasmid')


plt.tight_layout()
plt.savefig('/'.join([wdir,'results',prefix,'pOTU_prevalence_length_query_hit_lengthratio_mobility.pdf']))
perc_mobility.to_csv('/'.join([wdir,'results',prefix,'pOTU_mobility_by_family.csv']), sep='\t')
potu_summary.to_csv('/'.join([wdir,'results',prefix,'pOTU_prevalence_query_hit_lengthratio_mobility.csv']), sep='\t', index=False)


