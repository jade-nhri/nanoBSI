#!/usr/bin/env python3 
import os,sys 
import pandas as pd 
cwd=os.getcwd() 
mydirs=[x for x in os.listdir() if '_spe' in x and os.path.isdir(x)] 
print (mydirs) 

count=0 
for i in sorted(mydirs):
    if os.path.exists(os.path.join(i,'resfinderc1')):
        os.chdir(i+'/resfinderc1')
        if not os.path.exists('ResFinder_results_tab.txt'):
            continue
        count+=1
        df=pd.read_table('ResFinder_results_tab.txt')
        print (df)
        print (count)
        df['#File']='{0}_c1'.format(i)
        if count==1:
             dfout=df
        else:
            dfout=pd.concat((dfout,df),axis=0)
        os.chdir('../..')
    if os.path.exists(os.path.join(i,'resfinderc2')):
        os.chdir(i+'/resfinderc2')
        if not os.path.exists('ResFinder_results_tab.txt'):
            continue
        count+=1
        df=pd.read_table('ResFinder_results_tab.txt')
        print (df)
        print (count)
        df['#File']='{0}_c2'.format(i)
        if count==1:
             dfout=df
        else:
            dfout=pd.concat((dfout,df),axis=0)
            #dfout=dfout.append(df,ignore_index=True)
        os.chdir('../..')
        continue

    if os.path.exists(os.path.join(i,'resfinderc')):
        os.chdir(i+'/resfinderc')


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
    os.chdir('../..') 
print (os.getcwd()) 
print (dfout) 
dfout.to_csv('Resfinder_resistance_report.csv',index=False)
