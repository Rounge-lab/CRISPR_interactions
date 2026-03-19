#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Rarefy the databases to assess the saturation of the spacers target identification

"""

import pandas as pd
import seaborn as sb
from Bio import SeqIO
import random
import os
import matplotlib.pyplot as plt

bldir='PATH_TO/cctyper_snakemake_wf/data/BLASTRes_dereplicated'
wdir='PATH_TO_MANUS_FOLDER/results/db_rarefaction'

dbfol='PATH_TO_DATABASES/databases'
dblist={'CRISPRCasdb':['spacer_id.fsa'], 'Inphared':['Inphared01052023.fa'],'IMGPR':['IMGPR_nucl.fna'], 'PLSDB': ['plsdb.fna'],
        'CRCbiome-vOTUs_final':['repr_viral_seqs.fasta'],'CRCbiome-plasmids90':['plasmids_dereplicated_0.9.fasta'], 'SpacerDB':['nr_spacers_hq-all_25-05-10.fna']}

clustersloc='PATH_TO_MANUS_FOLDER/datasets/spacers_manus_table.csv'
relabloc='PATH_TO_MANUS_FOLDER/datasets'

#Extract headers of the fasta db entries
def get_refs(db,wdir=wdir):
    print(f'Extracting sequence IDs for {db}')
    os.makedirs('/'.join([wdir,db]),exist_ok=True)
    seqids=[]
    for seq in SeqIO.parse('/'.join([dbfol,db,dblist[db][0]]),'fasta'):
        seqids.append(seq.id)
        
    print(f"{db} is parsed, {len(seqids)} headers are now written to '/db_rarefaction/{db}/Db_seqIDs.csv'...")
    seqids=pd.DataFrame({'SeqID':seqids})
    seqids.to_csv('/'.join([wdir,db,f'{db}_seqIDs.csv']),sep='\t',index=False)
    print('Ids recorded')
    
for db in dblist:
    get_refs(db)

#Read the file and keep only those clusters that are in samples from the manus
clusters=pd.read_csv(clustersloc,sep='\t')
def read_file(db,bldir=bldir,clusters=clusters):
    f=pd.read_csv('/'.join([bldir,f'blastres_{db}_minid99_lr95.txt']),sep=',')
    f=f.loc[f['query'].isin(clusters['Cluster'].unique().tolist())]
    hits=f['hit'].tolist()
    return f,hits

#Rarefy database entries and check how many hits we would get with each db size

dbs={'CRCbiome-vOTUs':['CRCbiome-vOTUs_final'],'CRCbiome-plasmids':['CRCbiome-plasmids90'], 'Inphared':['Inphared'],
     'PLSDB':['PLSDB'], 'IMGPR':['IMGPR'],'CRISPRCas-db':['CRISPRCasdb'],'SpacerDB':['SpacerDB']}

random.seed=45 #set a seed to generate 10 random seeds for rarefaction
seeds=[random.randint(1,1000) for _ in range(10)]
   
def rarefy_db(db,blastres,hits,qhits,step,wdir=wdir,dbs=dbs, seeds=seeds):
    print(f'Rarefying data for {db}')
    rartab={'Db':[],'Seed':[],'NumEntries':[],'NumHits':[],'NumSpacers':[]}
    for s in seeds:
        print(f'going through seed {s}')
        random.seed=s
        c=1
        if db!='SpacerDB':
            lim=len(qhits)
        else:
            lim=10**6
        while c*step<=lim:
              rhits=random.sample(qhits,c*step)
              found=list(set(hits)&set(rhits)) #check how many records are detected in clusters
              rartab['Db'].append(db)
              rartab['Seed'].append(s)
              rartab['NumEntries'].append(c*step) #Db size
              rartab['NumHits'].append(len(found)) #get number of unique DB entries that spacers map to
              spacers=blastres.loc[blastres['hit'].isin(found)] #get list of spacers that match the list
              rartab['NumSpacers'].append(len(spacers['query'].unique().tolist())) #get number of clusters mapping to these hits
              c+=1
    rartab=pd.DataFrame(rartab)
    return rartab

Rartab=pd.DataFrame()
for db in dbs:
    blastres,hits=read_file(db)
    qhits=pd.read_csv('/'.join([wdir,dbs[db][0],f'{dbs[db][0]}_seqIDs.csv']))
    qhits=qhits['SeqID'].tolist()
    rartab=rarefy_db(db,blastres,hits,qhits,1000)
    rartab.to_csv('/'.join([wdir,dbs[db][0],f'{dbs[db][0]}_rarefaction.csv']),sep='\t',index=False)
    Rartab=pd.concat([Rartab,rartab])
    
Rartab.to_csv('/'.join([wdir,'Rarefaction_all_NumSpacers.csv']),sep='\t',index=False)

Rartab=pd.read_csv('/'.join([wdir,'Rarefaction_all_NumSpacers.csv']),sep='\t')

#Visualize data

ax=sb.lineplot(Rartab.query('NumEntries<100000'),x='NumEntries',y='NumSpacers', estimator='mean', ci='sd', marker='*',
             hue='Db', palette=['#6b519d','#4d4f50','#ff5757','#ff914d','#00bf63','#4b99cc','#ffde59'])

plt.legend(title='Database')
ax.set(xlabel='Database size, number of entries', ylabel='CRISPR spacers with homology to database entries')
plt.savefig('/'.join([wdir, 'Rarefaction_Databases.pdf']))


#Split CRCbiome-vOTUs to virus and provirus based on fraction of viral genomes identified as such
virprov=pd.read_csv(f"{wdir.replace('db_rarefaction','')}vOTUs_in_MAGs_contigs_clean.csv", sep=',')
checkv=pd.read_csv(f"{wdir.replace('results/db_rarefaction','datasets')}/viral_contigs_list_organized_onlymanus_votus.csv")

for ix,v in virprov.iterrows():
    query=checkv.loc[checkv['new_id']==v.vOTU]
    virprov.loc[ix,'ProvFraction']=len(query.loc[query['checkV']=='Yes'])/len(query)
    
virprov['ProvFraction'].describe()

virprov['IntegrationStatus']=virprov['ProvFraction'].apply(lambda row: 'Always' if row==1 else ('Never' if row==0 else 'Mixed'))
virprov.to_csv(f"{wdir.replace('db_rarefaction','')}vOTUs_in_MAGs_contigs_clean.csv", sep=',', index=False)

col='IntegrationStatus'
Rartab_vir=pd.DataFrame()
for st in ['Always','Never','Mixed']:
    blastres,hits=read_file('CRCbiome-vOTUs')
    qhits=virprov.loc[virprov[col]==st]
    qhits=qhits['vOTU'].tolist()
    rartab=rarefy_db('CRCbiome-vOTUs',blastres,hits,qhits,500)
    rartab[col]=st
    Rartab_vir=pd.concat([Rartab_vir,rartab])
    
Rartab_vir.to_csv('/'.join([wdir,f'Rarefaction_vOTUs_{col}_status.csv']),sep='\t',index=False)

ax=sb.lineplot(Rartab_vir,x='NumEntries',y='NumSpacers', estimator='mean', errorbar='sd', marker='*',
             hue=col, palette=['#6b519d','#00bf63','#af9a5b'])
plt.savefig('/'.join([wdir, f'Rarefaction_vOTUs_by_{col}.pdf']))

#Split CRCbiome-PTUs to mobilizable, non-mobilizable and conjugative
#Change step to a 100 because of the small number of mobilizable and conjugative plasmids
mobplas=pd.read_csv(f"{wdir.replace('db_rarefaction','')}PTUs_mobility_ARG_summary.tsv", sep='\t')

col='Mobility'
Rartab_plas=pd.DataFrame()
for st in ['mobilizable','conjugative','non-mobilizable']:
    blastres,hits=read_file('CRCbiome-plasmids')
    qhits=mobplas.loc[mobplas[col]==st]
    qhits=qhits['PTU'].tolist()
    rartab=rarefy_db('CRCbiome-plasmids',blastres,hits,qhits,20)
    rartab[col]=st
    Rartab_plas=pd.concat([Rartab_plas,rartab])
    
Rartab_plas.to_csv('/'.join([wdir,f'Rarefaction_PTUs_{col}_status.csv']),sep='\t',index=False)

ax=sb.lineplot(Rartab_plas,x='NumEntries',y='NumSpacers', estimator='mean', errorbar='sd', marker='*',
             hue=col, palette=[ '#af9a5b','#6b519d','#00bf63',], hue_order=['non-mobilizable','mobilizable','conjugative'])
ax.set(xlim=[0,800],ylim=[0,2000])
plt.savefig('/'.join([wdir, f'Rarefaction_PTUs_by_{col}.pdf']))

#----------------------------------------------------------------------------------------------------------
#Split data randomly, remove all vOTUs/PTUs that are only present in one group and 
# search the spacers of this group in the database
#----------------------------------------------------------------------------------------------------------

databases={'CRCbiome-vOTUs':[1000],'CRCbiome-plasmids':[1000]} #database: size for rarefaction step

random.seed=22
newseeds=[random.randint(1,1000) for _ in range(10)]

#get a random set of samples #check
Rartab=pd.DataFrame()
for db in databases:
       
    blastres,hits=read_file(db)
    
    if 'vOTU' in db:
        relab=pd.read_csv('/'.join([relabloc,'vOTUs_relab.tsv']),sep='\t').set_index('sample_id')
    else:
        relab=pd.read_csv('/'.join([relabloc,'PTUs_relab.tsv']),sep='\t').set_index('sample_id')
        
    prev=relab.copy()
    prev[prev>0]=1
    allsums=prev.sum(axis=0)
    
    pick=1
    for ns in newseeds:
        random.seed=ns #set a seed to generate 10 random seeds for rarefaction
        id_a=random.sample(relab.index.tolist(),517) #groupA: group where individual vOTUs are excluded, groupB: all other individuals
        #find prevalence of taxa in group A
        asums=prev.loc[id_a].sum(axis=0)
        #Find all vOTUs that are only in these samples (their prevalence in whole dataset is equal their prevalence in groupA)
        torem=allsums-asums
        torem=torem[torem==0].index.tolist()
        qhits=[q for q in prev.columns.tolist() if q not in torem] #keep only vOTUs from other samples (groupB) as queries in the database
        
        #Keep only spacers that are detected in the randomly picked samples
        aclust=clusters.loc[clusters['Sample'].isin(id_a),'Cluster'].unique().tolist()
        ablast=blastres.loc[blastres['query'].isin(aclust)]
        ahits=ablast['hit'].tolist()
        
        rartab=rarefy_db(db,ablast,ahits,qhits,databases[db][0])
        rartab['PopSubsample']=f'Subsample_{pick}'
        Rartab=pd.concat([Rartab,rartab])
        pick+=1
    
Rartab.to_csv('/'.join([wdir,'Rarefaction_PopulationSampling.csv']),sep='\t',index=False)

ax=sb.lineplot(Rartab,x='NumEntries',y='NumSpacers', estimator='mean', errorbar='sd', marker='*',
             hue='Db', palette=['#6b519d','#4d4f50'])
    

##Make a complete figure with data from full databases as well

Fulltab=pd.read_csv('/'.join([wdir,'Rarefaction_all_NumSpacers.csv']),sep='\t')

Rartab['Db']=Rartab['Db'].str.replace('vOTUs','vOTUs_IndivExc')
Rartab['Db']=Rartab['Db'].str.replace('plasmids','plasmids_IndivExc')

Fulltab=pd.concat([Fulltab,Rartab])

ax=sb.lineplot(Fulltab.query('NumEntries<100000'),x='NumEntries',y='NumSpacers', estimator='mean', ci='sd', marker='*',
             hue='Db', palette=['#6b519d','#4d4f50','#ff5757','#ff914d','#00bf63','#4b99cc','#ffde59'])

plt.legend(title='Database')
ax.set(xlabel='Database size, number of entries', ylabel='CRISPR spacers with homology to database entries')
plt.savefig('/'.join([wdir, 'Rarefaction_Databases_with_PopulationSubsampling.pdf']))