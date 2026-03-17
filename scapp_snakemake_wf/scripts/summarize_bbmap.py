#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 09:47:38 2024

Make a summary table for plasmids coverage

@author: p1068-ekateria
"""

import pandas as pd
import os
import argparse

par=argparse.ArgumentParser(description='Make summary table for BBMap results, provide parameters which data to look for')

#Arguments
par.add_argument('-w', help='working directory with BBMap coverage files')
par.add_argument('-t', help='target info to summarize on, f.ex. covered percent or number of reads')
par.add_argument('-o', help='path for the output file')

args=par.parse_args()

wdir=args.w
target=args.t
out=args.o

samlist=os.listdir(wdir)
samlist=[s for s in samlist if 'covstats' in s]

for s in samlist:
    query=pd.read_csv('/'.join([wdir,s]),sep='\t')
    query['#ID']=query['#ID'].apply(lambda row: row.split(' ')[0])
    query=query[['#ID', target]]
    query=query.set_index('#ID')
    colname=s.replace('_covstats.txt','')
    query=query.rename(columns={target:colname})
    if 'summary' not in globals():
        summary=query.copy()
    else:
        summary=summary.merge(query, left_index=True,right_index=True)

summary=summary.T

summary.to_csv(out, index=True)
