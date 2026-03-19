#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Find spacers ntersect between databases

"""

import pandas as pd
from itertools import combinations
from upsetplot import UpSet
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO

bldir='PATH_TO/cctyper_snakemake_wf/data/BLASTRes_dereplicated'

dbs=['CRCbiome-vOTUs','CRCbiome-plasmids','SpacerDB','CRISPRCas-db','Inphared','PLSDB', 'IMGPR']

def read_file(db,bldir=bldir):
    f=pd.read_csv('/'.join([bldir,f'blastres_{db}_minid99_lr95.txt']),sep=',')
    f['Database']=db
    return f

blastres=pd.DataFrame()
for db in dbs:
    file=read_file(db)
    blastres=pd.concat([blastres,file])


#Load the list of spacer clusters from the manuscript dataset 

clusters=pd.read_csv('PATH_TO_MANUS_FOLDER/datasets/spacers_manus_table.csv',sep='\t')

#Keep only those clusters that are in the manuscript database
blastres=blastres.loc[blastres['query'].isin(clusters['Cluster'].unique().tolist())]

numclus=pd.DataFrame()
for db in dbs:
    x=len(blastres.loc[blastres['Database']==db,'query'].unique().tolist())
    print(f'Number of clusters in {db}: {str(x)}')
    n=pd.DataFrame({'Db1':[db],'IntersectClusters':[x]})
    numclus=pd.concat([numclus,n])
    
for db in dbs:
    x=len(blastres.loc[blastres['Database']==db,'hit'].unique().tolist())
    print(f'Number of hits in {db}: {str(x)}')
    
for db in dbs:
    hits_range=pd.read_csv('/'.join([bldir,f'blastres_{db}_minid99_lr95_numhits.txt']),sep=',')
    hits_range=hits_range.loc[hits_range['Unnamed: 0'].isin(clusters['Cluster'].unique().tolist())]
    x=hits_range['query'].max()
    print(f'Maximum number of hits per cluster in {db}: {str(x)}')
    
def intersect(l1,l2):
    common=[l for l in l1 if l in l2]
    return common

print('Intersecting 2')
IntDatabase=pd.DataFrame()
for x,y in list(combinations(dbs,2)):
    xl=blastres.loc[blastres['Database']==x,'query'].unique().tolist()
    yl=blastres.loc[blastres['Database']==y,'query'].unique().tolist()
    com=intersect(xl,yl)
    db=pd.DataFrame.from_dict({'Db1':[x],'Db2':[y],'IntersectClusters':[len(com)]})
    IntDatabase=pd.concat([IntDatabase,db])
    
IntDatabase.to_csv('PATH_TO_MANUS_FOLDER/results/Spacers_Intersect_databases.csv',index=False, sep='\t')

final=IntDatabase.copy()

print('Intersecting 3')
IntDatabase=pd.DataFrame()
for x,y,z in list(combinations(dbs,3)):
    xl=blastres.loc[blastres['Database']==x,'query'].unique().tolist()
    yl=blastres.loc[blastres['Database']==y,'query'].unique().tolist()
    com=intersect(xl,yl)
    zl=blastres.loc[blastres['Database']==z,'query'].unique().tolist()
    com2= intersect(com,zl)
    db=pd.DataFrame.from_dict({'Db1':[x],'Db2':[y],'Db3':[z],'IntersectClusters':[len(com2)]})
    IntDatabase=pd.concat([IntDatabase,db])
    
final=pd.concat([final,IntDatabase])

final.to_csv('PATH_TO_MANUS_FOLDER/results/Spacers_Intersect_databases_all_3intONLY.csv',index=False, sep='\t')
print('Done with intersection')

##Make an UpSet plot
print('Making UpSet plot')
final=final.reset_index()
final=final.drop(columns=['index'])
for ix,f in final.iterrows():
    for db in dbs:
        if db in final.loc[ix,['Db1','Db2','Db3']].tolist():
            final.at[ix,db]=True
        else:
            final.at[ix,db]=False
           
multi_index = pd.MultiIndex.from_arrays([final[dbs[0]],final[dbs[1]],final[dbs[2]],final[dbs[3]],final[dbs[4]],final[dbs[5]],final[dbs[6]]], names=dbs)
final=pd.Series(final['IntersectClusters'].values, index=multi_index)

upset=UpSet(final,orientation='vertical',show_counts='{:d}',facecolor='#5fc0bf',element_size=15,sort_by='cardinality',totals_plot_elements=5, intersection_plot_elements=6)

upset.plot()
plt.savefig('PATH_TO_MANUS_FOLDER/results/Spacers_Intersect_databases_UpSet_vertical_3only.pdf')
plt.savefig('PATH_TO_MANUS_FOLDER/results/Spacers_Intersect_databases_UpSet_vertical_3only.jpg')

upset=UpSet(final,show_counts='{:d}',facecolor='#5fc0bf',element_size=None,sort_by='cardinality',totals_plot_elements=5, intersection_plot_elements=6)
upset.plot()
plt.savefig('PATH_TO_MANUS_FOLDER/results/Spacers_Intersect_databases_UpSet_horizontal_3only.pdf')
plt.savefig('PATH_TO_MANUS_FOLDER/results/Spacers_Intersect_databases_UpSet_horizontal_3only.jpg')

#Check the lengths of databases
dbfol='PATH_TO_MANUS_FOLDER/databases'
dblist={'CRISPRCasdb':['spacer_id.fsa'], 'Inphared':['Inphared01052023.fa'],'IMGPR':['IMGPR_nucl.fna'], 'PLSDB': ['plsdb.fna'],
        'CRCbiome-vOTUs_final':['repr_viral_seqs.fasta'],'CRCbiome-plasmids90':['plasmids_dereplicated_0.9.fasta'], 'SpacerDB':['nr_spacers_hq-all_25-05-10.fna']}

dblen={key: None for key in dblist.keys()}

for db in dblist:
    num_seqs=0
    for seq in SeqIO.parse('/'.join([dbfol,db,dblist[db][0]]),'fasta'):
        num_seqs +=1
    dblen[db]=num_seqs
    print(f'Number of seqs in {db}: {num_seqs}')
    
