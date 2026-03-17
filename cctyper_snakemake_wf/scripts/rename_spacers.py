#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 11:52:58 2023

Rename dereplicated sequences

@author: p1068-ekateria
"""
import argparse
from Bio import SeqIO
import pandas as pd

parser = argparse.ArgumentParser(description='Rename dereplicated seqs')
parser.add_argument('-i', '--infile', help='indicate file to be renamed', required=True)
parser.add_argument('-o', '--outfile', help='indicate new file', required=True)
parser.add_argument('-s', '--seqname', help='string to be used for renaming', required=True)

args=parser.parse_args()

infile = args.infile
outfile = args.outfile
seqname = args.seqname

def rename_sequences(infile, outfile, seqname):
    with open(infile, "r") as handle:
        records = list(SeqIO.parse(handle, "fasta"))

    old_ids=[]
    new_ids=[]
    for i, record in enumerate(records):
        old_ids.append(record.id)
        new_ids.append(seqname+'_'+str(i+1)+'_spacer_cctyper')
        record.id = f"{seqname}_{i+1}_spacer_cctyper"
        record.description = f"{seqname}_{i+1}_spacer_cctyper"
    
    maptab=pd.DataFrame({'NewID':new_ids, 'OldID':old_ids})
    
    with open(outfile, "w") as outhandle:
        SeqIO.write(records, outhandle, "fasta")

    maptab.to_csv(outfile.replace('.fasta','_map.csv'), index=False)
# Rename_sequences
rename_sequences(infile, outfile, seqname)