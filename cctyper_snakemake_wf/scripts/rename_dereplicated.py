#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Rename dereplicated sequences

"""
import argparse
import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Rename dereplicated seqs')
parser.add_argument('-i', '--infile', help='indicate file to be renamed', required=True)
parser.add_argument('-o', '--outfile', help='indicate new file', required=True)
parser.add_argument('-t', '--tag', help='tags on the end of the seqnames', required=True)

args=parser.parse_args()

infile = args.infile
outfile = args.outfile
tag=args.tag

def rename_sequences(infile, outfile, tag):
    
    #Find which sequence is mapped with what
    mapped={}
    key=None
    refs={}
    with open ('.'.join([infile,'clstr']),'r') as cluster:
        for line in cluster:
            line=line.strip()
            if line.startswith('>'):
                key=line[1:].replace(' ','_')
                mapped[key]=[]
                refs[key]=[]
            elif key is not None:
                x=line.split('\t')[-1].split(', >')[1]
                seqn='_'.join(x.split('_')[:2])
                mapped[key].append(seqn)
                if '. *' in x: #if a reference
                    refs[key].append(seqn)
                       
    mapped=pd.DataFrame(mapped.items(), columns=['Cluster','Spacers'])
    mapped.set_index('Cluster',inplace=True)
    refs=pd.DataFrame(refs).transpose().reset_index()
    refs.columns=['Cluster','Representative']
    refs['Representative']=refs['Representative'].apply(lambda row: row +'_' + tag)
    
    with open(infile, "r") as ffile:
        records = list(SeqIO.parse(ffile, "fasta"))

    for i, record in enumerate(records):
        clname=refs.loc[refs['Representative']==record.id,'Cluster'].tolist()[0]
        record.id = f"{clname}"
        #record.description = f"{clname}"
    
    #Write files
    with open(outfile, "w") as outhandle:
        SeqIO.write(records, outhandle, "fasta")
        
    refs.to_csv(outfile.replace('.fasta','_refs.csv'), index=False)
    mapped.to_csv(outfile.replace('.fasta','_clustermap.csv'))

# Rename_sequences
rename_sequences(infile, outfile, tag)