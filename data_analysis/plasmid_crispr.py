#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Find which MGEs have CRISPR elements

"""

import pandas as pd
from scipy.stats import binomtest

wdir='PATH_TO_MANUS_FOLDER'

spacers=pd.read_csv('/'.join([wdir,'datasets/spacers_manus_table.csv']), sep='\t')

plcrispr=pd.read_csv('/'.join([wdir,'results/CCTyper_plasmids/crisprs_all.tab']),sep='\t')
vircrispr=pd.read_csv('/'.join([wdir,'results/CCTyper_viruses/crisprs_all.tab']),sep='\t')

#Find number of plasmids with CRISPR cassettes (all)

plaslist=plcrispr['Contig'].unique().tolist()
virlist=vircrispr['Contig'].unique().tolist()

#Exclude those that are not detected in the manus

relab=pd.read_csv('/'.join([wdir, 'datasets/pOTUs_relab.tsv']),sep='\t')
vrelab=pd.read_csv('/'.join([wdir, 'datasets/vOTUs_relab.tsv']),sep='\t')

relab=relab.set_index('sample_id')
vrelab=relab.set_index('sample_id')


plaslist=[p for p in plaslist if p in relab.columns.tolist()]
virlist=[v for v in virlist if v in vrelab.columns.tolist()]

nrep=plcrispr['N_repeats'].describe()
nrep=vircrispr['N_repeats'].describe()


print('\n')
print(f"Number of viruses with CRISPR-Cas detected: {len(virlist)}")
print(f"Number of repeats per cassette: {nrep['mean'].round()}±{nrep['std'].round()}; [{nrep['min']}; {nrep['max']}]")

#Check what type of plasmids are those that contain CRISPR-Cas

mobtyper=pd.read_csv('/'.join([wdir,'datasets/plasmids_dereplicated_0.9_MOBtyper.txt']),sep='\t')
mobcr=mobtyper.loc[mobtyper['sample_id'].isin(plaslist)]
mobint=mobcr['predicted_mobility'].value_counts().reset_index()
mobtot=mobtyper['predicted_mobility'].value_counts()/len(mobtyper)
mobtot=mobtot.reset_index()
mobint=mobint.merge(mobtot,on='predicted_mobility')
mobint.columns=['predicted_mobility','Crispr','Total']
mobint['Sum']=len(mobcr)
mobint['binomp']=mobint.apply(lambda row: binomtest(row.Crispr, row.Sum, row.Total, alternative='two-sided').pvalue, axis=1)

#Check what type of viruses are those that contain CRISPR-Cas

virtypes=pd.read_csv('/'.join([wdir,'results/vOTUs_in_MAGs_contigs_clean.csv']),sep=',')
vircr=virtypes.loc[virtypes['vOTU'].isin(virlist)]
virint=vircr['IntegrationStatus'].value_counts().reset_index()
virtot=virtypes['IntegrationStatus'].value_counts()/len(virtypes)
virtot=virtot.reset_index()
virint=virint.merge(virtot,on='IntegrationStatus')
virint.columns=['IntegrationStatus','Crispr','Total']
virint['Sum']=len(vircr)
virint['binomp']=virint.apply(lambda row: binomtest(row.Crispr, row.Sum, row.Total, alternative='two-sided').pvalue, axis=1)
