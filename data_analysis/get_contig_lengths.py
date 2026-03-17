#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 11:57:51 2023

Get sequence lengths for all contigs in all samples

@author: ekateria
"""
from Bio import SeqIO
import pandas as pd
import os

resfol='PATH_TO_MANUS_FOLDER/datasets'
contpath='PATH_TO_ASSEMBLY_CONTIGS'

samples=os.listdir(contpath)
samples=[sam for sam in samples if 'S-' in sam]

def get_seqlen(fasta):
    
    lengths=dict()
    nofile=[]
    # Open the file
    try:
        with open(fasta, "r") as fs:
            # Iterate over each record in the FASTA file
            for record in SeqIO.parse(fs, "fasta"):
                # Get the sequence ID and length
                lengths[record.id] = [len(record.seq)]  
    
    except:
        nofile.append(sam)
        
    return lengths,nofile

All_contigs=pd.DataFrame()

for sam in samples:
    print(sam+'_contigs.fasta')
    contlen,nofile=get_seqlen('/'.join([contpath,sam,sam+'_contigs.fasta']))
    contlen=pd.DataFrame(contlen).T.reset_index()
    All_contigs=pd.concat([All_contigs,contlen],axis=0,ignore_index=True)

All_contigs.columns=['Contig','Length']

All_contigs.to_csv('/'.join([resfol,'All_contig_lengths.csv']),sep='\t',index=False)
