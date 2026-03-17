#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Find which CRISPR cassettes were detected in the final dataset samples and get some stats

"""

import pandas as pd
import seaborn as sb
import numpy as np
import pyreadr as pyr
import scipy.stats as stats
import matplotlib.pyplot as plt
from statsmodels.formula.api import ols
from statsmodels.stats.multitest import multipletests

crdir='PATH_TO/cctyper_snakemake_wf/data'
wdir='PATH_TO_MANUS_FOLDER' #set working directory

samples=pyr.read_r('/'.join([wdir, 'participant_data/sample_meta.Rds']))[None]
samples=samples.drop(columns=['richness', 'shannon', 'invsimpson'])
alpha_div=pd.read_csv('/'.join([wdir, 'datasets/MAGs_AlphaDiv.tsv']), sep='\t')

samples=samples.merge(alpha_div[['sample_id','kjonn','age_cat','age_invitation','senter','ObsSp','Shannon','InvSimpson']],on='sample_id',how='left')

#Add info on CRISPR proximity
samples['NumCRISPR_All']=None
samples['NumCRISPR_NearCas']=None
samples['NumCRISPR_Orphan']=None
samples['NumCRISPR_Putative']=None
samples['NumRepeats']=None

cassettes=pd.DataFrame() #make a dataframe for all samples
for i,sam in samples.iterrows():
    crfile=pd.read_csv('/'.join([crdir, 'CCtyper', sam.sample_id.replace('_','-'), 'crisprs_all_filtered.tab']), sep = '\t')
    cassettes=pd.concat([cassettes, crfile], ignore_index=True)
    samples.at[i,'NumCRISPR_All']=len(crfile)
    samples.at[i,'NumCRISPR_NearCas']=len(crfile.loc[crfile['Cas_proximity']=='near_cas'])
    samples.at[i,'NumCRISPR_Orphan']=len(crfile.loc[crfile['Cas_proximity']=='orphan'])
    samples.at[i,'NumCRISPR_Putative']=len(crfile.loc[crfile['Cas_proximity']=='putative'])
    samples.at[i,'NumRepeats']=crfile['N_repeats'].sum()
    
nfcassettes = pd.DataFrame()  # make a dataframe for all samples

samples['NumCRISPR_AllNF'] = None
samples['NumRepeatsNF'] = None

for i, sam in samples.iterrows():
    crall=pd.DataFrame()
    for t in ['near_cas','orphan','putative']:
        try:
            crfile = pd.read_csv(
                '/'.join([crdir, 'CCtyper', sam.sample_id.replace('_', '-'), f'crisprs_{t}.tab']), sep='\t')
            crfile['Cas_proximity']=t
            crall=pd.concat([crall,crfile])
            nfcassettes = pd.concat([nfcassettes, crfile], ignore_index=True)
        except:
            print(f'Sample {sam.sample_id} does not have {t} crisprs')
    samples.at[i, 'NumCRISPR_AllNF'] = len(crall)
    samples.at[i,'NumCRISPR_NearCas']=len(crall.loc[crall['Cas_proximity']=='near_cas'])
    samples.at[i,'NumCRISPR_Orphan']=len(crall.loc[crall['Cas_proximity']=='orphan'])
    samples.at[i,'NumCRISPR_Putative']=len(crall.loc[crall['Cas_proximity']=='putative'])
    samples.at[i, 'NumRepeatsNF'] = crall['N_repeats'].sum()

nfcassettes.to_csv('/'.join([wdir, 'datasets/Combined_crisprs_all_NOT_filtered.tab']), sep = '\t', index=False)

#Convert to integers
samples['NumCRISPR_All']=samples['NumCRISPR_All'].astype(int)
samples['NumCRISPR_NearCas']=samples['NumCRISPR_NearCas'].astype(int)
samples['NumCRISPR_Orphan']=samples['NumCRISPR_Orphan'].astype(int)
samples['NumCRISPR_Putative']=samples['NumCRISPR_Putative'].astype(int)
samples['NumRepeats']=samples['NumRepeats'].astype(int)

#Save data
samples.to_csv('/'.join([wdir,'datasets/data_by_samples_crispr.csv']), sep = '\t', index=False)

#Read the data
samples=pd.read_csv('/'.join([wdir,'datasets/data_by_samples_crispr.csv']),sep = '\t')
cassettes=pd.read_csv('/'.join([wdir, 'datasets/Combined_crisprs_all_filtered.tab']), sep = '\t')
nfcassettes=pd.read_csv('/'.join([wdir, 'datasets/Combined_crisprs_all_NOT_filtered.tab']), sep = '\t')


#Avg, median, std stats
samples['NumCRISPR_AllNF'].describe()

#Make a histogram
fig=sb.histplot(data=samples,x='NumCRISPR_All',color='#52a1b5')
plt.tight_layout()
fig.set(xlabel='Number of CRISPR cassettes per sample')
plt.savefig('/'.join([wdir, 'results/CRISPR_Cassettes_histogram.pdf']))
plt.close()

#Find how many spacers are in cassettes of 10 repeats and longer
over10=nfcassettes.query('N_repeats>=10')['N_repeats'].sum()
below10=nfcassettes.query('N_repeats<10')['N_repeats'].sum()
print(f"Total number of repeats: {nfcassettes['N_repeats'].sum()}")
print(f"Number of repeats in cassettes >=10: {over10}")
print(f"Number of repeats in cassettes <10: {below10}")
print(f"Fraction spacers filtered: {below10/(over10+below10)*100:.2f}")



#Find correlation between metadata and CRISPR

#technical
labels=['Number of cassettes','Total bases, QC','Processed reads', 'Number of contigs', 'N50']
tcorr=samples[['NumCRISPR_AllNF','Total_Bases_QC_ATLAS','reads_proc','n_contigs','N50']].corr()
mask = np.tril(np.ones_like(tcorr, dtype=bool))
fig=sb.heatmap(tcorr, annot=True,cmap='coolwarm', mask=mask,center=0.5, xticklabels=labels,yticklabels=labels)
fig.set_xticklabels(fig.get_xticklabels(),fontsize=9)
fig.set_yticklabels(fig.get_yticklabels(),fontsize=9)
plt.savefig('/'.join([wdir, 'results/Cassettes_technical_corr.pdf']))
plt.close()

#biological_microbiome
labels=['Number of cassettes','GC, %', 'Richness', 'InvSimpson', 'Shannon']
bcorr=samples[['NumCRISPR_AllNF','GC_perc','ObsSp','InvSimpson', 'Shannon']].corr()
fig=sb.heatmap(bcorr, annot=True,cmap='Blues',center=0.5, xticklabels=labels,yticklabels=labels)
fig.set_xticklabels(fig.get_xticklabels(),fontsize=9)
fig.set_yticklabels(fig.get_yticklabels(),fontsize=9)
plt.savefig('/'.join([wdir, 'results/Cassettes_biological_corr.pdf']))
plt.close()

#Make a scatter plot of number of repeats by sequencing depth
fig, axes = plt.subplots(2, 1, figsize=(5,12))
sb.regplot(data=samples, x='NumCRISPR_All',y='Total_Bases_QC_ATLAS',color='#5fc0bf',ax=axes[0])
#axes[0].set_yscale('log')
axes[0].set_ylim([10**9, 1.5*10**10])
axes[0].set_ylabel('Sequencing depth,bp')
axes[0].set_xlabel('Number of CRISPR cassettes')

sb.regplot(data=samples, x='NumRepeats',y='Total_Bases_QC_ATLAS',color='#5fc0bf',ax=axes[1])
#axes[1].set_yscale('log')
axes[1].set_ylim([10**9, 1.5*10**10])
axes[1].set_ylabel('Sequencing depth,bp')
axes[1].set_xlabel('Number of CRISPR repeats')
plt.savefig('/'.join([wdir, 'results/Num_cassettes_repeats_vs_depth.pdf']))
plt.close()

samples['Perc_NearCasNF']=samples['NumCRISPR_NearCas']/samples['NumCRISPR_All']*100

def cor_depth(samples,cols,target): 
    peardf=pd.DataFrame()
    for c in cols:
        pr, p = stats.pearsonr(samples[c],samples[target])
        k=pd.DataFrame.from_dict({'Column': [c], 'Pearsonr': [pr], 'Pval': [p]})
        peardf=pd.concat([peardf,k])
    
    _, padj, _, _ = multipletests(peardf['Pval'], method='fdr_bh')
    peardf['FDRp']=padj
    
    return peardf

cols=['NumCRISPR_AllNF','NumCRISPR_NearCas','NumCRISPR_Orphan','NumRepeats','Perc_NearCasNF']
peardf=cor_depth(samples,cols,'Total_Bases_QC_ATLAS')
peardf.to_csv('/'.join([wdir, 'results/Cassettes_Pearson_SeqDepth.tsv']),sep='\t',index=False)


cols=['n_contigs','N50','ObsSp','GC_perc']
peardf=cor_depth(samples,cols,'NumCRISPR_All')
peardf.to_csv('/'.join([wdir, 'results/Cassettes_Pearson_ObsSp_GC.tsv']),sep='\t',index=False)


#Make a summary of which CRISPR Types are detected in all samples
categ=['Unknown','I-A','I-B','I-C','I-D','I-E','I-F', 'I-F_T','I-G','III-A','III-B','III-C',
       'III-D','IV-A1','IV-A2','IV-A3','IV-D','II-A','II-B','II-C','V-A','V-F1','V-F2','VI-A','VI-B1','VI-B2','VI-C','VI-D']
cassettes['Prediction']=pd.Categorical(cassettes['Prediction'],categ)

fig=sb.jointplot(data=cassettes,x='N_repeats',y='Prediction',marginal_ticks=True,kind='hist', bins=28, color='#5fc0bf')
fig.ax_joint.set_xlabel('Number of repeats')
fig.ax_joint.set_ylabel('CRISPR type')
fig.ax_marg_x.hist(cassettes['N_repeats'],bins=28, color='#5fc0bf',edgecolor='black')
fig.ax_marg_x.set_yscale('log')
fig.ax_marg_y.set_xscale('log')
plt.savefig('/'.join([wdir, 'Cassettes_by_type_JointPlot.pdf']))
plt.close()

#Assign class to the samples

def assign_class(val):
    if 'I-' in val or 'III-' in val or 'IV-' in val:
        return('Class1')
    else:
        if 'Unknown' in val:
            return('Unknown')
        else:
            return('Class2')
            
nfcassettes['Class']=nfcassettes['Prediction'].apply(assign_class)
nfcassettes['Class'].value_counts()
nfcassettes['Prediction'].value_counts()

nfcassettes['N_repeats']=nfcassettes['N_repeats'].astype(int)
nfcassettes['N_repeats'].describe()

nfcassettes['Repeat_len'].describe()

nfcassettes['Type']=nfcassettes['Prediction'].apply(lambda row: row.split('-')[0])
nfcassettes['Type']=pd.Categorical(nfcassettes['Type'], categories=['I','III','IV','II','V','VI','Unknown'], ordered=True)

#Plot by type
fig=sb.jointplot(data=nfcassettes,x='N_repeats',y='Type',marginal_ticks=True,kind='hist', bins=28, color='#5fc0bf')
fig.ax_joint.set_xlabel('Number of repeats')
fig.ax_joint.set_ylabel('CRISPR type')
fig.ax_marg_x.hist(nfcassettes['N_repeats'],bins=35, color='#5fc0bf',edgecolor='black')
fig.ax_marg_x.set_yscale('log')
fig.ax_marg_y.set_xscale('log')
plt.savefig('/'.join([wdir, 'results/Cassettes_by_type_JointPlot.pdf']))
plt.close()


#Add sampleID to all_samples
cassettes['Sample']=cassettes['Contig'].apply(lambda row: row.split('_')[0])
last_column = cassettes.pop('Sample')
cassettes.insert(0, 'SampleID', last_column)

#Check which cassettes were detected inside MAGs
magcont=pyr.read_r('/'.join([wdir,'datasets/MAGs_contig_list.Rds']))[None]

MAGcont=pd.DataFrame()
for ix,mag in magcont.iterrows():
    c=mag.Contigs
    c=c.replace('[','')
    c=c.replace(']','')
    c=c.replace("'",'')
    c=c.split(', ')
    c=pd.DataFrame({'Contig':c})
    c['MAG']=mag.MAG
    MAGcont=pd.concat([MAGcont,c])  

nfcassettes=nfcassettes.merge(MAGcont,on='Contig', how='left')

nfcassettes['inMAG']=nfcassettes['Contig'].apply(lambda row: 'Yes' if row in MAGcont['Contig'].unique().tolist() else 'No')
nfcassettes['inMAG'].value_counts()
nfcassettes.to_csv('/'.join([wdir, 'datasets/Combined_crisprs_all_NOT_filtered.tab']), sep = '\t', index=False)

#Find if there is a correlation between length of contig and orphan or near cas CRISPR cassette

all_contigs=pd.read_csv('/'.join([wdir, 'datasets/All_contig_lengths.csv']),sep='\t')
nfcassettes=nfcassettes.merge(all_contigs,on='Contig',how='left')
nfcassettes=nfcassettes.rename(columns={'Length':'ContigLength'})
del all_contigs #remove variable because it's too large

nfcassettes.to_csv('/'.join([wdir, 'datasets/Combined_crisprs_all_NOT_filtered.tab']), sep = '\t', index=False)

#Plot proximity vs length (boxplots)
fig=sb.boxplot(data=nfcassettes,x='ContigLength',y='Cas_proximity',color='#5fc0bf')
fig.set_xscale('log')
fig.set_xlabel('Contig length, bp')
fig.set_ylabel('Cas proximity')
plt.savefig('/'.join([wdir, 'results/Cas_proximity_vs_length.pdf']))
plt.close()

#Check if there is significant difference in contig length between groups
def kruskal_group(data,cat_col,y):
    categ=data[cat_col].unique().tolist()
    for_stats={}
    for cat in categ:
        for_stats[cat]=data.loc[data[cat_col]==cat,y].tolist()

    h_st,pval=stats.kruskal(*list(for_stats.values())) # *indicates to treat each element in list as a group, the list will contain as many elements as there are dict entries, each entry is a separate list within a list

    return h_st,pval

prox=['near_cas','orphan','putative']
KruskalProx=pd.DataFrame()
for p in prox:
    test=[i for i in prox if i!=p]
    h_st,pval=kruskal_group(nfcassettes.loc[nfcassettes['Cas_proximity']!=p],'Cas_proximity','ContigLength')
    k=pd.DataFrame.from_dict({'Proximity': ['_vs_'.join(test)], 'Hstat': [h_st], 'Pval': [pval]})
    KruskalProx=pd.concat([KruskalProx,k])

_, padj, _, _ = multipletests(KruskalProx['Pval'], method='fdr_bh')
KruskalProx['FDRp']=padj
KruskalProx.to_csv('/'.join([wdir, 'Kruskal_CasProximity_ContigLength.tsv']),sep='\t',index=False)

h_st,pval=kruskal_group(nfcassettes,'Cas_proximity','ContigLength')


## Check difference between age_cat, sex, region and a/b use, adjusted for seq depth

meta=pd.read_csv('/'.join([wdir,'participant_data/screening_data.tsv']), sep='\t')
meta['beforeBL'] = pd.Categorical(meta['beforeBL'], categories=['Yes','No'], ordered=True)

samples=samples.merge(meta[['deltaker_id','beforeBL']], on='deltaker_id', how='left')

def diff_adjusted(samples, y, col, adj):
    
    model = ols(f'{y} ~ C({col}) + {adj}', data=samples).fit()
    
    print(model.summary())

    results=pd.DataFrame({'Y':[y], 'Group':[col], 'Adj':[adj], 'Rsq':[model.rsquared], 'RsqAdj':[model.rsquared_adj],
                          'Fstat':[model.fvalue],'Prob_F':[model.f_pvalue],'PvalGroupVar':[model.pvalues[1]],'PvalAgjVar':[model.pvalues[2]],
                          'PvalIntercept':[model.pvalues[0]]})
    
    return results

OLS_differ=pd.DataFrame()
for y in ['NumCRISPR_AllNF', 'NumRepeats']:
    for gr in ['beforeBL','age_cat','kjonn','senter']:
        olsres=diff_adjusted(samples, y, gr, 'Total_Bases_QC_ATLAS')
        OLS_differ=pd.concat([OLS_differ, olsres])
        
OLS_differ.to_csv('/'.join([wdir, 'OLS_cassettes_richness_SeqDepthAdj.tsv']),sep='\t',index=False)
    
                        
    