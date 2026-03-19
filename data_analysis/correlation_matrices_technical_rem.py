#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spearman correlation between bacteria and MGEs

"""

import pandas as pd
import pyreadr as pyr
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt

wdir='PATH_TO_MANUS_FOLDER' #set working directory

#-----------------------------------------------------------------------------------------------
##Load the data

#Metadata
meta=pd.read_csv('/'.join([wdir,'participant_data/screening_data.tsv']), sep='\t')
meta['beforeBL'] = pd.Categorical(meta['beforeBL'], categories=['Yes','No'], ordered=True)
meta['kjonn'] = pd.Categorical(meta['kjonn'], categories=['Male','Female'], ordered=True)
meta['age_cat'] = pd.Categorical(meta['age_cat'], categories=['50-59','60-69','>=70'], ordered=True)

samples=pyr.read_r('/'.join([wdir, 'participant_data/sample_meta.Rds']))[None]

#Relative abundances
mags=pd.read_csv('/'.join([wdir, 'datasets/MAGs_relab.tsv']),sep='\t').set_index('sample_id')
votus=pd.read_csv('/'.join([wdir, 'datasets/vOTUs_relab.tsv']),sep='\t').set_index('sample_id')
PTUs=pd.read_csv('/'.join([wdir, 'datasets/PTUs_relab.tsv']),sep='\t').set_index('sample_id')

#Taxonomy
mag_taxon=pd.read_csv('/'.join([wdir,'datasets/MAG_taxonomy_full.tsv']), sep=',')
votu_taxon=pd.read_csv('/'.join([wdir,'datasets/vOTUs_taxonomy.csv']),sep=',')
votu_taxon=votu_taxon.drop(columns='Unnamed: 0')
votu_taxon=votu_taxon.rename(columns={'Scaffold':'vOTU'})
PTU_taxon=pd.read_csv('/'.join([wdir,'datasets/PTU_taxonomy_IMGPR_0.95_modified.txt']))

#-----------------------------------------------------------------------------------------------
##Filter the data - remove technical correlations

#Filter based on the prevalence (keep only those that are >=20% prevalent in the population)

def filter_prev(relab,domain):
    prev=relab.copy()
    if domain!='MAG':
        prev[prev>0]=1 #keep all results over 0 since each virus/plasmid are covered on at least 75% of its length
    else:
        prev[prev>=0.001]=1 #set the threshold to 0.001 % (it makes 10000 bp out of a minimum threshold of 1Gb data per sample)
        prev[prev<0.001]=0
    
    tokeep=prev.sum(axis=0)
    tokeep=tokeep[tokeep>=20]
    prevtaxa=tokeep.index.tolist()
    
    relab=relab[prevtaxa]
    print(f'Number of {domain}s to keep: {len(tokeep)}')
    
    return relab

mags=filter_prev(mags,'MAG')
votus=filter_prev(votus,'vOTU')
PTUs=filter_prev(PTUs,'PTU')

#Filter based on cross-correlation to other taxa within a domain

def find_corr(data):
    dcorr=data.corr(method='spearman')
    dcorr_un=dcorr.unstack().reset_index()
    dcorr_un=dcorr_un.rename(columns={"level_0": "Tax1", "level_1": "Tax2", 0: "SpCorrCoef"})
    dcorr_un=dcorr_un.query('SpCorrCoef != 1') #remove self-comparison
    return dcorr_un, dcorr

def rem_technical(data):
    dcorr=find_corr(data)[0]
    dcorr=dcorr.query('SpCorrCoef >= 0.99')
    #Keep only unique pairs
    dcorr['Combination']=dcorr[['Tax1','Tax2']].apply(' '.join,axis=1)
    dcorr['Combination']=dcorr['Combination'].apply(lambda x: ''.join(sorted(x)))
    dcorr=dcorr.drop_duplicates(subset='Combination',keep='first')
    dcorr=dcorr.drop(columns='Combination')
    return dcorr

# Find technical correlation among MAGs
mc=rem_technical(mags) #get the list of correlating mags
mc=pd.merge(mc,mag_taxon[['MAG','species']],right_on='MAG',left_on='Tax1')#Find their taxonomy
mc=pd.merge(mc,mag_taxon[['MAG','species']],right_on='MAG',left_on='Tax2')#Find their taxonomy
mc=mc.rename(columns={'species_x':'Species1', 'species_y':'Species2'})
mc=mc.drop(columns=['MAG_x','MAG_y'])

mags=mags.drop(columns=mc['Tax2'])

#Find technical correlation among viruses
vc=rem_technical(votus)
vc=pd.merge(vc,votu_taxon[['vOTU','Family']],right_on='vOTU',left_on='Tax1')#Find their taxonomy
vc=pd.merge(vc,votu_taxon[['vOTU','Family']],right_on='vOTU',left_on='Tax2')#Find their taxonomy
vc=vc.rename(columns={'Family_x':'Family1', 'Family_y':'Family2'})
vc=vc.drop(columns=['vOTU_x','vOTU_y'])

votus=votus.drop(columns=vc['Tax2'])

#Find technical correlation among plasmids
pc=rem_technical(PTUs)
pc=pd.merge(pc,PTU_taxon[['PTU','Hit_family']],right_on='PTU',left_on='Tax1')#Find their taxonomy
pc=pd.merge(pc,PTU_taxon[['PTU','Hit_family']],right_on='PTU',left_on='Tax2')#Find their taxonomy
pc=pc.rename(columns={'Hit_family_x':'Family1', 'Hit_family_y':'Family2'})
pc=pc.drop(columns=['PTU_x','PTU_y'])

PTUs=PTUs.drop(columns=pc['Tax2'])

mc.to_csv('/'.join([wdir,'results/Cross_correlating_MAGs.csv']),index=False)
vc.to_csv('/'.join([wdir,'results/Cross_correlating_vOTUs.csv']),index=False)
pc.to_csv('/'.join([wdir,'results/Cross_correlating_PTUs.csv']),index=False)


##-------------------------------------------------------------------------
#Spearman correlation between domains

relabs={'MAGs':mags, 'vOTUs':votus, 'PTUs':PTUs}
taxonomy={'MAGs':mag_taxon[['MAG','species']],'vOTUs':votu_taxon[['vOTU','Family']],'PTUs':PTU_taxon[['PTU','Hit_family']]}

def correlate_domains(dom1,dom2, corlim, relabs=relabs):
    full=pd.concat([relabs[dom1],relabs[dom2]],axis=1,join='inner')

    #Calculate correlation
    print(f'Calculating Spearman correlation for {dom1} vs {dom2}')
    spcorr_un, spcorr=find_corr(full)
    spcorr=spcorr.reindex(index=relabs[dom1].columns.values,columns=relabs[dom2].columns.values)
    spcorr.to_csv('/'.join([wdir, f'results/Relab_Spearman_correlation_{dom1}_{dom2}_square.csv']),index=True)

    spcorr_un=spcorr_un.loc[spcorr_un['Tax1'].isin(relabs[dom1].columns.values)] #keep only dom1 as Tax1
    spcorr_un=spcorr_un.loc[spcorr_un['Tax2'].isin(relabs[dom2].columns.values)] #keep only dom2 as Tax2
    
    #add taxonomy data
    tax1=taxonomy[dom1]
    tax2=taxonomy[dom2]
    tax1.columns=['Tax1','Taxonomy1']
    tax2.columns=['Tax2','Taxonomy2']
    
    spcorr_un=spcorr_un.merge(tax1, on='Tax1', how='left')
    spcorr_un=spcorr_un.merge(tax2, on='Tax2', how='left')

    plt.plot(figsize=(4,6))
    sb.histplot(spcorr_un, x='SpCorrCoef')
    plt.tight_layout()
    spcorr_un.to_csv('/'.join([wdir,f'results/Relab_Spearman_correlation_{dom1}_{dom2}_melted.csv']),index=False)

    spcorr_un['AbsSpCor']=spcorr_un['SpCorrCoef'].apply(lambda x: abs(x))
    spcorr_un=spcorr_un.sort_values(by='AbsSpCor', ascending=False)
    
    spcorr_strong=spcorr_un.query(f'AbsSpCor>{corlim}')
    spcorr_strong=spcorr_strong.drop(columns='AbsSpCor')

    spcorr_strong.to_csv('/'.join([wdir,f'results/Relab_Spearman_correlation_{dom1}_{dom2}_AbsOver{str(corlim)}.csv']),index=False)

    return spcorr, spcorr_un, spcorr_strong

dompairs=[('MAGs','vOTUs'),('MAGs','PTUs')]

for d in dompairs:
    spcorr,spcorr_un, spcorr_strong=correlate_domains(d[0],d[1],0.8)

##-----------------------------------------------------------------------------------------------------
## Make a clustermap of domains that strongly correlate
def make_labels(data,col1,col2):
    hl=pd.DataFrame(data[col1].unique().tolist(),columns=[col1])
    hl=hl.merge(data[[col1,col2]],on=col1,how='left')
    hl=hl.drop_duplicates()
    hl['label']=hl[col1]+'_'+hl[col2]
    hl=hl.sort_values(by=col2)
    
    return hl

def make_clustermap(data,hm,hv):
    heat=data.reindex(index=hm['Tax1'],columns=hv['Tax2'])
    heat=heat.merge(hm[['Tax1','label']],left_index=True,right_on='Tax1')
    heat.set_index('label',inplace=True)
    heat=heat.drop(columns='Tax1')
    heat.columns=hv['label']
    ax=sb.clustermap(heat,cmap='Blues',linewidths=2)
    ax.ax_heatmap.set_xlabel('')
    ax.ax_heatmap.set_ylabel('')
    ax.ax_heatmap.set_title('')
    #ax.ax_heatmap.set_yticks(np.arange(0,len(hm),1))
    #ax.ax_heatmap.set_xticks(np.arange(0,len(hv),1))
    ax.ax_heatmap.set_yticklabels(ax.ax_heatmap.get_yticklabels(),fontdict={'fontstyle':'italic'})
    ax.ax_heatmap.set_xticklabels(ax.ax_heatmap.get_xticklabels(),fontdict={'fontstyle':'italic'})

dom1='MAGs'; dom2='PTUs'; corlim=0.8
spcorr=pd.read_csv('/'.join([wdir, f'results/Relab_Spearman_correlation_{dom1}_{dom2}_square.csv'])).set_index('Unnamed: 0')
spcorr_strong=pd.read_csv('/'.join([wdir,f'results/Relab_Spearman_correlation_{dom1}_{dom2}_AbsOver{str(corlim)}.csv']))

hm=make_labels(spcorr_strong,'Tax1','Taxonomy1')
hv=make_labels(spcorr_strong,'Tax2','Taxonomy2')
make_clustermap(spcorr,hm,hv)
plt.tight_layout()
plt.savefig('/'.join([wdir,'results/Relab_Spearman_correlation_{dom1}_{dom2}_AbsOver{str(corlim)}.pdf']))

##-----------------------------------------------------------------------------------------
#For plasmids, check negative correlations 
dom1='MAGs'; dom2='PTUs'
spcorr_un=pd.read_csv('/'.join([wdir, f'results/Relab_Spearman_correlation_{dom1}_{dom2}_melted.csv']))
sb.histplot(spcorr_un, x='SpCorrCoef')
spcorr_neg=spcorr_un.query('SpCorrCoef<0')

negstat=spcorr_neg['SpCorrCoef'].describe()
sb.histplot(spcorr_neg, x='SpCorrCoef')

spcorr_negtail=spcorr_un.query('SpCorrCoef<-0.1')
spcorr_negtail=spcorr_negtail.sort_values(by='SpCorrCoef',ascending=True)

#Check which families tend to have negative correlations
#Add family for the MAG

def check_tax_correspondence(corrdata, mag_taxon,vislim):
    corrdata=corrdata.merge(mag_taxon[['MAG','family']], left_on='Tax1', right_on='MAG', how='left')
    corrdata=corrdata.drop(columns={'MAG'})

    #Remove plasmids with unknown family delineation
    corrdata=corrdata.dropna(subset='Taxonomy2')
    corrdata=corrdata.loc[~corrdata['Taxonomy2'].isin(['Unknown', 'Unclassified'])]
    corrdata=corrdata.rename(columns={'Taxonomy2': 'PTU family', 'family': 'MAG family'})
    corrdata['MAG family']='MAG_'+corrdata['MAG family']
    corrdata['PTU family']='PTU_'+corrdata['PTU family']

    conttab=pd.crosstab(corrdata['MAG family'], corrdata['PTU family'])
    #Keep rows and columns where there are more than 20 pairs in total
    sumcol=conttab.sum(axis=0)
    vistab=conttab[sumcol[sumcol>vislim].index.tolist()]
    
    sumcol=vistab.sum(axis=1)
    vistab=vistab.reset_index()
    vistab=vistab.loc[vistab['MAG family'].isin(sumcol[sumcol>vislim].index.tolist())]
    vistab=vistab.set_index('MAG family')

    g=sb.clustermap(vistab, fmt='d', annot=None, cmap=sb.cubehelix_palette(as_cmap=True))
    # Modify tick labels for both the x and y axes
    g.ax_heatmap.set_xlabel('')
    g.ax_heatmap.set_ylabel('')

    for label in g.ax_heatmap.get_xticklabels():
        label.set_fontsize(8)  # Set font size for x-axis tick labels
        label.set_fontstyle('italic')  # Set font style to italic
    
    for label in g.ax_heatmap.get_yticklabels():
        label.set_fontsize(8)  # Set font size for y-axis tick labels
        label.set_fontstyle('italic')  # Set font style to italic

    plt.ylabel('')
    plt.xlabel('')
    
    return conttab

vislim=100
negconttab=check_tax_correspondence(spcorr_neg, mag_taxon, vislim)
plt.savefig('/'.join([wdir,'results/Negative_correlations_taxonomy_conttab_MAGs_PTUs.pdf']))

#Find number of plasmids that belong to dif classes, which MAGs they target
fortax=spcorr_negtail.drop_duplicates(subset='Tax2')
fortax['Taxonomy2'].value_counts()

#Check positive correlations
spcorr_pos=spcorr_un.query('SpCorrCoef>0')

posstat=spcorr_pos['SpCorrCoef'].describe()
sb.histplot(spcorr_pos, x='SpCorrCoef')

spcorr_postail=spcorr_un.query('SpCorrCoef>0.1')
vislim=500
posconttab=check_tax_correspondence(spcorr_postail, mag_taxon, vislim)
plt.savefig('/'.join([wdir,'results/Positive_correlations_taxonomy_conttab_MAGs_PTUs.pdf']))
