#!/usr/bin/env python3
import os,sys
import pandas as pd
import numpy as np
import operator
from ete3 import NCBITaxa
ncbi = NCBITaxa()
import warnings
warnings.filterwarnings("ignore")
warnings.simplefilter("ignore")
#ncbi.update_taxonomy_database()

infile=sys.argv[1]

def taxid2speid(taxid):
    speid='0'
    try:
        lineage=ncbi.get_lineage(taxid)
        for i in lineage:
            ranki=ncbi.get_rank([i])
            if ranki[i]=='species':
                speid=i
    except:
        pass
    return(speid)


df=pd.read_table(infile,sep='\t')
#print (df)
df=df[df['queryLength']>=1000]
#print (df)
#print(len(np.unique(df['readID'].values)))
unitaxid=np.unique(df['taxID'].values)
#print (len(unitaxid))

tidvalue=dict()
tax2spe=dict()
for i in unitaxid:
    #print (i)
    dfset=df[df['taxID']==i]
    #print (dfset)
    tidvalue[i]=sum(1/dfset['numMatches'])
    #print (sum(1/dfset['numMatches']))
    if i==0:
        tax2spe['0']=[0]
        continue
    else:
        sid=taxid2speid(i)
        if sid not in tax2spe.keys():
            tax2spe[sid]=[i]
        else:
            tax2spe[sid].append(i)
#print (tax2spe)
#print (tidvalue)

sidvalue=dict()
allsum=0
for sid in tax2spe.keys():
    taxid=tax2spe[sid]
    if len(taxid)>1:
        value=0
        for i in taxid:
            tempv=tidvalue[i]
            value=value+tempv
    else:
        value=tidvalue[taxid[0]]
    sidvalue[sid]=value
    allsum=allsum+value
#print (sidvalue)
#print (allsum)

sortd=sorted(sidvalue.items(),key=operator.itemgetter(1),reverse=True)
#print (sortd)
speid=[]
spename=[]
readcount=[]
for i,j in sortd:
    speid.append(i)
    if i=='0':
        spename.append('Unclassified')
    else:
        sname=''.join(ncbi.get_taxid_translator([i]).values())
        spename.append(sname)
    readcount.append(j)

for i in range(0,5):
    if spename[i]!='Unclassified':
        print ('{0} with {1} reads ({2}%)'.format(spename[i],readcount[i],readcount[i]/allsum*100))

fw=open(infile.replace('_cenout.txt','_cenout1.txt'),'w')
fw.write('taxID\tSpecies\n')
for i,j in zip(speid[0:5],spename[0:5]):
    fw.write('{0}\t{1}\n'.format(i,j))
fw.close()
#print (speid[0:5])
#print (spename[0:5])

