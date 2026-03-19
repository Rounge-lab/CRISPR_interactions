#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Make files for cytoscape network visualization

"""

import pandas as pd
import seaborn as sb
import numpy as np

wdir='PATH_TO_MANUS_FOLDER'
spacers=pd.read_csv('/'.join([wdir,'datasets/spacers_manus_table_new.csv']), sep='\t')

mag_taxon=pd.read_csv('/'.join([wdir,'datasets/MAG_taxonomy_full.tsv']), sep=',')
vir_taxon=pd.read_csv('/'.join([wdir,'datasets/vOTUs_taxonomy.csv']),sep=',')
PTU_taxon=pd.read_csv('/'.join([wdir, 'datasets/PTU_taxonomy_IMGPR_0.9_modified.txt']),sep=',')

#Keep only spacers that are found in mags
spacers=spacers.dropna(subset='MAG')

mags=pd.read_csv('/'.join([wdir, 'datasets/MAGs_relab.tsv']),sep='\t').set_index('sample_id')
votus=pd.read_csv('/'.join([wdir, 'datasets/vOTUs_relab.tsv']),sep='\t').set_index('sample_id')
PTUs=pd.read_csv('/'.join([wdir, 'datasets/PTUs_relab.tsv']),sep='\t').set_index('sample_id')

spacers=spacers.loc[spacers['MAG'].isin(mags.columns.tolist())]

edges=spacers[['Cluster','MAG','vOTUs','PTUs']].drop_duplicates()
edges=pd.melt(edges, id_vars=['Cluster','MAG'], value_vars=['vOTUs','PTUs'], var_name='domain', value_name='taxon')
edges=edges.dropna(subset='taxon')
edges['taxon']=edges['taxon'].apply(lambda row: row.replace("'",''))
edges['taxon']=edges['taxon'].apply(lambda row: row.replace("[",''))
edges['taxon']=edges['taxon'].apply(lambda row: row.replace("]",''))
edges['taxon']=edges['taxon'].str.split(', ')
edges = edges.explode('taxon')

#remove those that target vOTUs/PTUs in CRCbiome, but not in manus:
mges=votus.columns.tolist()
mges.extend(PTUs.columns.tolist())
edges=edges.loc[edges['taxon'].isin(mges)]

#Remove edges to PTUs/vOTUs that contained crispr-cas cassettes
plcrispr=pd.read_csv('/'.join([wdir,'results/CCTyper_plasmids_rerun/crisprs_all.tab']),sep='\t')
vircrispr=pd.read_csv('/'.join([wdir,'results/CCTyper_viruses/crisprs_all.tab']),sep='\t')
edges=edges.loc[~edges['taxon'].isin(plcrispr['Contig'].unique().tolist())]
edges=edges.loc[~edges['taxon'].isin(vircrispr['Contig'].unique().tolist())]

edges['comb']=edges['MAG']+'_vs_'+edges['taxon']

edges=edges['comb'].value_counts()
edges=pd.DataFrame(edges.reset_index())
edges[['MAG', 'Taxon']]=edges['comb'].str.split('_vs_',expand=True)
edges=edges.rename(columns={'count': 'NumInteractions'})
edges=edges[['MAG','Taxon','NumInteractions']]
edges['Type']=edges['Taxon'].apply(lambda row: 'vOTU' if 'vOTU' in row else 'PTU')

sb.histplot(edges, x='NumInteractions')

#Add info on correlations between MAGs, PTUs and vOTUs

pcorr=pd.read_csv('/'.join([wdir, 'results/Relab_Spearman_correlation_MAGs_PTUs_melted.csv']))
vcorr=pd.read_csv('/'.join([wdir, 'results/Relab_Spearman_correlation_MAGs_vOTUs_melted.csv']))

def add_corr(edges,corr,letter):
    corr['comb']=corr['Tax1']+'_vs_'+ corr['Tax2']
    edges=edges.merge(corr[['comb','SpCorrCoef']],on='comb', how='left')
    edges[f'Correlation{letter}']=edges['SpCorrCoef'].apply(lambda row: 'positive' if row>0 else ('negative' if row<0 else row))
    edges=edges.drop(columns='SpCorrCoef')
    return edges

edges['comb']=edges['MAG']+'_vs_'+ edges['Taxon']
edges=add_corr(edges,pcorr,'p')
edges=add_corr(edges,vcorr,'v')
edges['Correlation']=edges.apply(lambda row: row.Correlationp if 'plasmid' in row.Taxon else row.Correlationv, axis=1)
edges=edges.drop(columns=['Correlationp','Correlationv','comb'])

edges.to_csv('/'.join([wdir,'results/Edges_network_crisprfree.csv']),index=False)

#Make a file for nodes

mag_taxon=pd.read_csv('/'.join([wdir,'datasets/MAG_taxonomy_full.tsv']), sep=',')
votu_taxon=pd.read_csv('/'.join([wdir,'datasets/vOTUs_taxonomy.csv']),sep=',')
votu_taxon=votu_taxon.drop(columns='Unnamed: 0')
votu_taxon=votu_taxon.rename(columns={'Scaffold':'vOTU'})
PTU_taxon=pd.read_csv('/'.join([wdir, 'datasets/PTU_taxonomy_IMGPR_0.9.txt']),sep=',')
PTU_taxon=PTU_taxon.query('LongestMatchPairwiseId>=90')
PTU_taxon['host_taxonomy']=PTU_taxon['host_taxonomy'].fillna('Unclassified')
PTU_taxon['Hit_family']=PTU_taxon['host_taxonomy'].apply(lambda row: row.split('f__')[1] if 'f__' in row else (row if row=='Unclassified' else 'Unknown'))
PTU_taxon['Hit_family']=PTU_taxon['Hit_family'].apply(lambda row: row.split(';')[0] if ';' in row else row)
PTU_taxon=PTU_taxon.rename(columns={'Query':'PTU'})

def make_nodes(edges, otherfilt, mag_taxon=mag_taxon):
    nodes=pd.DataFrame(pd.concat([edges['MAG'],edges['Taxon']]),columns=['Taxon'])
    nodes=nodes.drop_duplicates(keep='first')
    nodes['Domain']=nodes['Taxon'].apply(lambda row: 'MAG' if 'MAG' in row else ('vOTU' if 'vOTU' in row else 'PTU'))


    nodes=nodes.merge(mag_taxon[['MAG','family']], left_on='Taxon',right_on='MAG', how='left')
    nodes=nodes.drop(columns=['MAG'])

    nodes['Taxonomy']=nodes.apply(lambda row: row.family if row.Domain=='MAG' else None, axis=1)

    nodes=nodes[['Taxon','Domain','Taxonomy']]
    
    magtoplot=pd.DataFrame(nodes['Taxonomy'].value_counts()).reset_index()
    magtoplot.columns=['Taxonomy','Count']
    magtoplot['Label']=magtoplot.apply(lambda row: row.Taxonomy if row.Count>=otherfilt else 'Other', axis=1)

    nodes=nodes.merge(magtoplot[['Taxonomy','Label']],on='Taxonomy',how='left')
    nodes['Label']=nodes.apply(lambda row: 'vOTU' if 'vOTU' in row.Domain else ('PTU' if 'PTU' in row.Domain else row.Label), axis=1)
    
    return nodes

nodes=make_nodes(edges,10)

#Add info on virus/provirus and plasmid mobility
vir_comp=pd.read_csv('PATH_TO_VIRUSES/dereplication/checkV_summary.tsv', sep='\t')
vir_rename=pd.read_csv('PATH_TO_VIRUSES/dereplication/old_to_new_ids.tsv', sep='\t')
vir_comp=vir_comp.merge(vir_rename, left_on='virus_id',right_on='old_id', how='left')
vir_comp=vir_comp[['new_id','completeness','provirus']].rename(columns={'new_id':'Taxon','completeness':'Completeness','provirus':'status'}).dropna(subset='Taxon')
PTU_mob=pd.read_csv('/'.join([wdir, 'datasets/plasmids_dereplicated_0.9_MOBtyper.txt']),sep='\t').rename(columns={'sample_id':'Taxon','predicted_mobility':'status'})

status=pd.concat([vir_comp[['Taxon','status']],PTU_mob[['Taxon','status']]])
nodes=nodes.merge(status[['Taxon','status']], on='Taxon', how='left')
nodes.loc[nodes['status']=='Yes','status']='Provirus'
nodes.loc[nodes['status']=='No','status']='Virus'

nodes=nodes.drop_duplicates()

nodes.to_csv('/'.join([wdir,'results/Nodes_network_crisprfree.csv']),index=False)

#Add how many individuals each connection is detected in

numpairs=pd.read_csv('/'.join([wdir,'results/Interactions_counts_by_spacer_and_indiv.csv']))
numpairs=numpairs.rename(columns={'pair':'Pair'})

edges=edges.merge(numpairs,on='Pair',how='left')
edges['5orMore']=edges['NumInd'].apply(lambda row: 'Yes' if row>=5 else 'No')

edges.to_csv('/'.join([wdir,'results/Edges_network_crisprfree.csv']),index=False)

#Keep only those that have connection to more than one MAG

counts=edges['Taxon'].value_counts()
counts=counts[counts>1].index.tolist()
filt_edges=edges.loc[edges['Taxon'].isin(counts)]
filt_nodes=make_nodes(filt_edges,2)
filt_nodes=filt_nodes.merge(nodes[['Taxon','status']], on='Taxon', how='left')

filt_edges.to_csv('/'.join([wdir,'results/Edges_network_atleast2MAG_CF.csv']),index=False)
filt_nodes.to_csv('/'.join([wdir,'results/Nodes_network_atleast2MAG_CF.csv']),index=False)
