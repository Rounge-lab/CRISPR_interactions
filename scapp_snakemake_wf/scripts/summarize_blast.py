#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Filter BLAST results

"""

import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description='Rename dereplicated seqs')
parser.add_argument('-i', '--infile', help='file with BLAST hits', required=True)
parser.add_argument('-m', '--metadatafile', help='Info on PLSDB metadata', required=True)
parser.add_argument('-o', '--outfile', help='Summary_file', required=True)

args=parser.parse_args()

infile = args.infile
outfile = args.outfile
meta=args.metadatafile

metadata=pd.read_csv(meta, sep='\t')
metadata.columns=metadata.columns.str.replace('_NUCCORE','')

blast=pd.read_csv(infile, sep='\t', header=None)
blast.columns=['query','hit','id','len','qlen','mismatch','gapopen', 'evalue']

#Summarize blast table

def summarize(blast, metadata):
    
    sumtab=pd.DataFrame()
    queries=blast['query'].unique().tolist()
    for q in queries:
        hits=blast.loc[blast['query']==q]
        row=[q]
        row.extend([np.mean(hits['qlen'])]) 
        row.extend(hits['hit'].unique().tolist())    
        row.extend([max(hits['len'])])
        row.append(hits.loc[hits['len']==max(hits['len']),'id'].tolist()[0]) #in case ids are identical, extract the 1st value

        #Add all data together
        sumtab=sumtab.append(pd.Series(row, index=['Query','QueryLength','Best_Hit','LongestMatch','LongestMatchPairwiseId']), ignore_index=True)
        
    #Keep only those that are covered >= 20 %
    sumtab['Match_ratio']=sumtab['LongestMatch']/sumtab['QueryLength']*100
    sumtab=sumtab.query('Match_ratio>=20')
    
    #Add metadata
    if 'ACC' in metadata.columns.tolist(): #if PLSDB
        metadata=metadata.rename(columns={'ACC':'Best_Hit','Topology':'Hit_topology','Completeness': 'Hit_completeness','Length':'Hit_length'})

        sumtab=sumtab.merge(metadata[['Best_Hit','Hit_topology','Hit_completeness','Hit_length','Description','taxon_family_name',
                                  'taxon_order_name','taxon_class_name','taxon_phylum_name']],
                                    on='Best_Hit', how='left')
    
        sumtab.columns=sumtab.columns.str.replace('taxon_','Hit_')
        sumtab.columns=sumtab.columns.str.replace('_name','')

        sumtab['QueryHit_ratio']=sumtab['QueryLength']/sumtab['Hit_length']*100
        sumtab=sumtab.sort_values(by='Query',ascending=True)
        
    elif 'ptu' in metadata.columns.tolist(): #if IMGPR
        metadata=metadata.rename(columns={'plasmid_id':'Best_Hit','topology':'Hit_topology','Hit_completeness':'putatively_complete',
                                          'length':'Hit_length'})
        metadata['Best_Hit']=metadata.apply(lambda row: '|'.join([row.Best_Hit,str(row.taxon_oid),str(row.scaffold_oid)]), axis=1)
        metadata=metadata.drop(columns={'ptu','taxon_oid','scaffold_oid'})
        
        sumtab=sumtab.merge(metadata, on='Best_Hit', how='left')
    
    return sumtab

sumtab=summarize(blast,metadata)

print('Number of plasmids with hits: '+str(len(sumtab)))
print('Number of hit plasmids from database: '+str(len(sumtab['Best_Hit'].unique())))

sumtab.to_csv(outfile, index=False)