#!/usr/bin/env python3
import os,sys
import pandas as pd

cwd=os.getcwd()
mydirs=[x for x in os.listdir() if '_res' in x and 'barcode' in x]
print (mydirs)
count=0
for i in sorted(mydirs):
    os.chdir(i)
    if not os.path.exists('ResFinder_results_tab.txt'):
        continue
    count+=1
    df=pd.read_table('ResFinder_results_tab.txt')
    print (df)
    print (count)
    df['#File']='{0}'.format(i)
    if count==1:
         dfout=df
    else:
        dfout=pd.concat((dfout,df),axis=0)
        #dfout=dfout.append(df,ignore_index=True)
    os.chdir('..')
print (os.getcwd())
print (dfout)
dfout.to_csv('Resfinder_resistance_report.csv',index=False)
