#!/usr/bin/env python3
import os,sys
import pandas as pd
import numpy as np
import subprocess

inpaf=sys.argv[1]
thr=int(sys.argv[2])
refnames=sys.argv[3]

def coverage(tstart,tend):
    pos=[]
    for m,n in zip(tstart,tend):
        temp=[x for x in range(m,n)]
        #print (temp)
        pos.extend(temp)
        #print (pos)
    #print (len(pos)/len(np.unique((pos))))
    depth=len(pos)/len(np.unique((pos)))
    covpos=len(np.unique((pos)))
    return(depth,covpos)


df=pd.read_table(inpaf,names=['Qname','Qlen','Qstart','Qend','strand','Tname','Tlen','Tstart','Tend','Nmatch','Alen','MQ1','MQ2','MQ3','MQ4','MQ5','MQ6','MQ7']) 
print (df)

df1=df[df['Qlen']>=thr]
print(df1)
df2=df1[df1['Alen']/df1['Qlen']>=0.8]
print(df2)
print (sum(df2['Alen'].values))

#refs=np.unique(df2['Tname'].values)
#print (refs)
refdf=pd.read_table(refnames,names=['Ref'])
refs=refdf['Ref'].values

refcov=dict()
for i in refs:
    dfset=df2[df2['Tname']==i].drop_duplicates('Qname')
    #print (dfset)
    if len(dfset)==0:
        continue
    if len(dfset['Tlen'])>0:
        gsize=dfset['Tlen'].values[0]
    else:
        gsize=dfset['Tlen'].values
    #print (gsize)
    tstart=dfset['Tstart'].values
    tend=dfset['Tend'].values
    qnames=dfset['Qname'].values
    depth,covpos=coverage(tstart,tend)
    #print (depth)
    #print (covpos/gsize)
    cov=covpos/gsize
    if cov>=0.75:
        reffile=refnames.replace('_names.txt','_Dref.fasta')
        comm='seqkit grep -p "{0}" {1} >> {2}'.format(i,refnames.replace('_names.txt','_refs.fasta'),reffile)
        #print (comm)
        subprocess.getoutput(comm)
        print (i)
        print (depth)
        print (cov*100)
        df2=df2[~df2['Qname'].isin(qnames)]

