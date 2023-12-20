#!/usr/bin/env python3
import os,sys
import pandas as pd
import numpy as np
import operator
from ete3 import NCBITaxa
ncbi = NCBITaxa()
#ncbi.update_taxonomy_database()

infile=sys.argv[1]

def taxid2speid(taxid):
    speid='0'
    lineage=ncbi.get_lineage(taxid)
    for i in lineage:
        ranki=ncbi.get_rank([i])
        if ranki[i]=='species':
            speid=i
    return(speid)



df=pd.read_table(infile,sep='\t')
print (df)
#print(len(np.unique(df['readID'].values)))
unitaxid=np.unique(df['taxID'].values)
#print (len(unitaxid))
#dfref=df[df['numMatches']==1]
dfref=df.drop_duplicates('readID')
dfref['taxID']=dfref['taxID'].astype('int')
dfref['seqID_taxID']=dfref['seqID']+':'+dfref['taxID'].astype('str')
print (dfref)

refs=dfref['seqID_taxID'].values
print (refs)
refcount=dict()
for i in refs:
    if i not in refcount.keys():
        refcount[i]=1
    else:
        refcount[i]+=1
refcount=sorted(refcount.items(),key=operator.itemgetter(1),reverse=True)
print (refcount[:10])
speout=[]
accout=[]
spename=[]

fw=open(infile.replace('_cenout.txt','_names.txt'),'w')
fw1=open(infile.replace('_cenout.txt','_cenout2.txt'),'w')
for i in range(0,10):
    temp=refcount[i]
    #print (temp)
    acci=temp[0].split(':')[0]
    taxi=temp[0].split(':')[1]
    #print (taxi)
    #print (ncbi.get_taxid_translator([taxi]))
    if taxi=='0': continue
    else:
        spei=taxid2speid(taxi)
        #print (spei)
    if len(speout)==0:
        speout.append(spei)
        accout.append(acci)
        tmp=''.join(ncbi.get_taxid_translator([spei]).values())
        spename.append(tmp)
        #print (acci)
        fw.write(acci+'\n')
        fw1.write('{0}\t{1}\t{2}\n'.format(spei,acci,tmp))
    else:
        if spei not in speout:
            #print (acci)
            fw.write(acci+'\n')
            speout.append(spei)
            accout.append(acci)
            tmp=''.join(ncbi.get_taxid_translator([spei]).values())
            spename.append(tmp)
            fw1.write('{0}\t{1}\t{2}\n'.format(spei,acci,tmp))


        else:
            continue
        if len(speout)>=5:
            break

#print (speout)
#print (spename)
fw.close()
fw1.close()

'''
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
print (sidvalue)
print (allsum)


sortd=sorted(sidvalue.items(),key=operator.itemgetter(1),reverse=True)
print (sortd)
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
    print ('{0} with {1} reads ({2}%)'.format(spename[i],readcount[i],readcount[i]/allsum*100))

#print (speid[0:5])
#print (spename[0:5])
'''
