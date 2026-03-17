#!/usr/bin/env python
"""
Filter CRISPRs based on the number of repeats, and add info on proximity to cas genes

"""
import argparse
import pandas as pd
import os
import pickle

parser = argparse.ArgumentParser(description='Filter CRISPRs', epilog='Remove CRISPRs with <10 repeats, add info on cas proximity')
parser.add_argument('-w', '--workfile', help='indicate sample folder where the files are', required=True)
args=parser.parse_args()

workfile=args.workfile
workfol=os.path.dirname(workfile)

#Read_files
file_all=pd.read_csv(workfile, sep = '\t')
file_all['Cas_proximity']=None

index=['putative', 'near_cas', 'orphan']
for ind in index:
    if os.path.isfile('/'.join([workfol, 'crisprs_'+ind+'.tab'])):
        toadd=pd.read_csv('/'.join([workfol, 'crisprs_'+ind+'.tab']), sep = '\t')['CRISPR'].tolist()
        #add location of CRISPR on the contig relative to the cas
        file_all.loc[file_all['CRISPR'].isin(toadd), 'Cas_proximity']=ind
        del toadd

#Filter CRISPRs that have less than 10 spacers
file_all=file_all.query('N_repeats>=10')

#Save the file
file_all.to_csv('/'.join([workfol, 'crisprs_all_filtered.tab']), sep = '\t', index=False)

with open('/'.join([workfol, 'crisprs_all_filtered.pkl']),'wb') as pickle_file:
    pickle.dump(file_all, pickle_file)

#Add info for blast files
blast=file_all['CRISPR'].tolist()
blast=['/spacers/' + element for element in blast]

with open('/'.join([workfol, 'spacers_to_blast.txt']),'w') as blfile:
    for element in blast:
        blfile.write(element + "\n")