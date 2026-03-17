#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 16:11:00 2023

Filter BLAST results

@author: p1068-ekateria
"""

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Rename dereplicated seqs')
parser.add_argument('-i', '--infile', help='input file with blast results', required=True)
parser.add_argument('-m', '--minid', help='minimum id to the reference', required=True)
parser.add_argument('-r', '--ratio', help='spacer to reference length ratio', required=True)

args=parser.parse_args()

infile = args.infile
minid=int(args.minid)
ratio=int(args.ratio)

outfile=infile.replace('.txt','_minid'+str(minid)+'_lr'+str(ratio)+'.txt')

blast=pd.read_csv(infile, sep='\t', header=None)
blast.columns=['query','hit','id','len','qlen','mismatch', 'qseq', 'sseq', 'evalue']

#Filter based on the identity
blast['id']=blast['id'].astype(float)
blast=blast.loc[blast['id'] >= minid]

#Filter based on the length ratio
blast['Ratio']=blast['len']/blast['qlen']*100
blast=blast.loc[blast['Ratio']>=ratio]

blast.to_csv(outfile, index=False)

def summarize(blast):
    queries=blast['query'].unique().tolist()
    hits=blast['hit'].unique().tolist()

    sumtab=pd.DataFrame({'UniqueQueries':[len(queries)], 'UniqueHits':[len(hits)]})
    num_hits=pd.DataFrame(blast['query'].value_counts())
    return sumtab,num_hits

sumtab,num_hits=summarize(blast)

sumtab.to_csv(outfile.replace('.txt','_stats.txt'), index=False)
num_hits.to_csv(outfile.replace('.txt','_numhits.txt'), index=True)