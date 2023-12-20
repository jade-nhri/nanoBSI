#!/usr/bin/env python3
import pandas as pd
import sys,os
import json

infile=sys.argv[1]
outfile=sys.argv[2]
df=pd.read_table(infile,keep_default_na=False)
print (df)

gene2note=dict()
for i,j in zip(df['Gene_accession no.'],df['Notes']):
    #print (i)
    #print ('{0}\t{1}'.format(i,j))
    temp=i.split('_')
    if len(temp)==3:
        gene=temp[0]+'_'+temp[2]
    if len(temp)==4:
        gene=temp[0]+'_'+temp[2]+'_'+temp[3]
    if gene not in gene2note.keys():
        gene2note[gene]=j
    else:
        print(gene)
#print (gene2note)
gene2note['Point_Point']='PointFinder'
gene2note['Intrinsic_Intrinsic']='Intrinsic'

print (len(gene2note.keys()))
with open(outfile,'w') as ofile:
    json.dump(gene2note,ofile,ensure_ascii=False)
    ofile.write('\n')
