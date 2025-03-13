#!/usr/bin/env python3
#-*-coding utf-8-*-
import pandas as pd
import numpy as np
import os,json,sys
import argparse
import webbrowser

fileVersion="20250313"
parser = argparse.ArgumentParser()
parser.add_argument("-i",
                    "--inputDirectory",
                    type=str,
                    help="input directory path")

parser.add_argument("-d",
                    "--DB",
                    type=str,
                    default="/opt/nanoBSI/antibiotic_name.json",
                    help="input DB path")
parser.add_argument("-D",
                    "--drug",
                    type=str,
                    default="/opt/nanoBSI/drug_list.xlsx",
                    help="drug_list file")      

parser.add_argument("-v",
                    "--version",
                    action='version',
                 
                    version=f'version:\t{fileVersion}')
args = parser.parse_args()
print (args.inputDirectory)
#print (args.DB)



#sys.exit()
column_name=['Barcode','Species name','Credibility']

df_drug=pd.read_excel(args.drug)
drug_list=df_drug['drug_list'].tolist()

column_name.extend(drug_list)



#讀取之前做好的json
# load json
f=open(f'{args.DB}')
dict_name=json.load(f)

#to load the dictionary about note of gene
f=open(f'/opt/gene2note.json')
gene2note=json.load(f)


#read barcode & species  data
df_species=pd.read_csv(f'{args.inputDirectory}/Species_identification_wconf.txt',names=['barcode','species','Credibility'],sep='\t')

#read Resfinder file
df_resistance=pd.read_csv(f'{args.inputDirectory}/Resfinder_resistance_report.csv')

#抓取要用的兩列，但要進行copy,不然之後會出warring
#select ['Resistance gene','Position in contig','Phenotype','Accession no.','#File'] field and copy to new DataFrame
df_tmp1=df_resistance.iloc[:,[0,6,7,8,9]].copy()
#print (df_tmp1)

#read Phenotype & Class data
df_Resfinder_DB=pd.read_excel('/opt/nanoBSI/organizeResfinderDb.xlsx')

#改成需要的名稱 
#change we want to use columns's name
df_tmp1.columns=['Resistance gene','Position in contig','Phenotype','Accession','barcode']
#print (df_tmp1)

#把barcode的文字處理成可以和檔案合併的文字
#split df_tmp1['barcode']  string for merge with df_species
#df_tmp1['barcode']=df_tmp1['barcode'].str.split('_',expand=True)[0]

#把兩個文件合併起來，方便之後利用
#merge  df_species and  merge by 'barcode'
df_tmp2=df_species.merge(df_tmp1,on='barcode',how='left')
#print (df_tmp2)

#更改欄位名稱
#change columns name
df_fin=pd.DataFrame(columns=['barcode','Species','Phenotype'])
#print (df_tmp2)
df_tmp2['Gene_accession']=df_tmp2['Resistance gene']+'_'+df_tmp2['Accession']
#print (df_tmp2)


#####這意思應該是
#print (df_tmp2['Position in contig'][0])
if 'NA..NA'!=df_tmp2['Position in contig'][0]:
    df_tmp2=df_tmp2.drop_duplicates(subset=['barcode','Position in contig','Phenotype'])
    #print (df_tmp2)

uniBarcode=df_species['barcode'].drop_duplicates().to_list()
uniBarcode=sorted(uniBarcode)

#顯示出 pdf有的Phenotype
#show all pdf's Phenotype 
for barcode in uniBarcode:
    #print (barcode)
    havePhenotype=False
    #df_tmpSpecies=df_tmp2[df_tmp2['species']==speciesName]
    df_tmpBarcode=df_tmp2[df_tmp2['barcode']==barcode].dropna(subset=['Phenotype'])
    if  not df_tmpBarcode.empty :
        tmpSpeciesName=df_tmpBarcode.iloc[0,1]
        list_uniSpeciesPhenotype=set(df_tmpBarcode['Phenotype'].str.split(', ').sum())
    else :
        print (f'{barcode} doesn\'t have data in Resfinder_resistance_report')
        df_tmpBarcode=df_species[df_species['barcode']==barcode]
        tmpSpeciesName=df_tmpBarcode.iloc[0,1]
        df_fin.loc[len(df_fin.index)]=[barcode,tmpSpeciesName,'']
        continue
    try:
        for phenotype in list_uniSpeciesPhenotype:
            if phenotype in dict_name[tmpSpeciesName] :
                df_fin.loc[len(df_fin.index)]=[barcode,tmpSpeciesName,phenotype]
                havePhenotype=True
###判斷例外，之後可能用function來處理，這樣比較不會讓程式碼太長
###user  try: & except 

        if 'Ampicillin' in list_uniSpeciesPhenotype and 'Ampicillin+Sulbactam' in dict_name[tmpSpeciesName]:
            df_fin.loc[len(df_fin.index)]=[barcode,tmpSpeciesName,'Ampicillin+Sulbactam']
            havePhenotype=True
        if 'Sulfamethoxazole' in list_uniSpeciesPhenotype and 'Trimethoprim' in list_uniSpeciesPhenotype:
            df_fin.loc[len(df_fin.index)]=[barcode,tmpSpeciesName,'Sulfamethoxazole']
            df_fin.loc[len(df_fin.index)]=[barcode,tmpSpeciesName,'Trimethoprim']
            df_fin.loc[len(df_fin.index)]=[barcode,tmpSpeciesName,'Trimethoprim+Sulfamethoxazole']
            havePhenotype=True
        if not havePhenotype:
            df_fin.loc[len(df_fin.index)]=[barcode,tmpSpeciesName,'']
    except KeyError  as e:
        print ('{} not in dict'.format(e))
        for phenotype in list_uniSpeciesPhenotype:
            if phenotype in dict_name['noFindSpecies'] :
                df_fin.loc[len(df_fin.index)]=[barcode,tmpSpeciesName,phenotype]
                havePhenotype=True
        if 'Ampicillin' in list_uniSpeciesPhenotype :
            df_fin.loc[len(df_fin.index)]=[barcode,tmpSpeciesName,'Ampicillin+Sulbactam']
            
            havePhenotype=True
        if 'Sulfamethoxazole' in list_uniSpeciesPhenotype and 'Trimethoprim' in list_uniSpeciesPhenotype:
            df_fin.loc[len(df_fin.index)]=[barcode,tmpSpeciesName,'Sulfamethoxazole']
            df_fin.loc[len(df_fin.index)]=[barcode,tmpSpeciesName,'Trimethoprim']
            df_fin.loc[len(df_fin.index)]=[barcode,tmpSpeciesName,'Trimethoprim+Sulfamethoxazole']                  
        if not havePhenotype:
            df_fin.loc[len(df_fin.index)]=[barcode,tmpSpeciesName,'']

    except :
        print ('Other Error')



df_fin=df_fin.merge(df_Resfinder_DB,on='Phenotype',how='left')

#print (uniBarcode)

count=0
dfout=pd.DataFrame(columns=['barcode','Credibility','Species','Phenotype','Class','Note'])
#print (df_species)
for barcode in uniBarcode:
    #print (barcode)
    dftmpout=pd.DataFrame(columns=['barcode','Credibility','Species','Phenotype','Class','Note'])
    
    subdf=df_fin[df_fin['barcode']==barcode]
    spe=subdf['Species'].values[0]
    #print (spe)
    subdf2=df_tmp2[df_tmp2['barcode']==barcode]
    #print (subdf)
    #print (subdf2)
    pheno=subdf['Phenotype'].values
    clas=subdf['Class'].values
    Credibility=subdf2['Credibility'].values[0]
    
    #print (pheno)
    if pheno[0] != '':
        for i,c in zip(pheno,clas):
#            print (i)
            if i =='Ampicillin+Sulbactam' :
                i='Ampicillin'
            for j,k in zip(subdf2['Phenotype'],subdf2['Gene_accession']):
                #print('{0}\t{1}'.format(j,k))
                temp=j.split(',')
                for m in temp:
                    m=m.replace(' ','')
                    if m==i:
                        count+=1
                        if count==1:
                            #print ('{0}\t{1}\t{2}\t{3}\t{4}'.format(barcode,spe,i,c,k))
                            dftmpout=pd.DataFrame([[barcode,Credibility,spe,i,c,str(k+':'+gene2note[k])]],columns=['barcode','Credibility','Species','Phenotype','Class','Note'])
                        else:
                            dfnew=pd.DataFrame([[barcode,Credibility,spe,i,c,str(k+':'+gene2note[k])]],columns=['barcode','Credibility','Species','Phenotype','Class','Note'])
                            #print ('{0}\t{1}\t{2}\t{3}\t{4}'.format(barcode,spe,i,c,k))
    
                            dftmpout=pd.concat([dftmpout,dfnew],ignore_index=True)
                    #else:
                        #print (f'####need check:\t{barcode} \t {m}')

        #print (dfout)
    #if barcode not in dfout['barcode'].values:
    else:
        count+=1
        tmpdf=pd.DataFrame([[barcode,Credibility,spe,'','','']],columns=['barcode','Credibility','Species','Phenotype','Class','Note'])
        dftmpout=pd.concat([dftmpout,tmpdf],ignore_index=True)
        #print (dfout)
    detailPhenotypeList=dftmpout['Phenotype'].to_list()
    #print(detailPhenotypeList)
    subdf_noDetail=subdf[~subdf['Phenotype'].isin(detailPhenotypeList)]
    #print(subdf_noDetail)
    dftmpout=pd.concat([dftmpout,subdf_noDetail])    
    dftmpout['Credibility']=Credibility
    
    
    dfout=pd.concat([dfout,dftmpout])
    
dfoutout=dfout.drop_duplicates(subset=['barcode','Phenotype','Note'])
#print (dfoutout)


df_fin=df_fin.sort_values(by=['barcode','Species','Class','Phenotype'])

df_fin.to_csv(f'{args.inputDirectory}/Resistance_prediction.csv',index=False)
df_fin.to_html(f'{args.inputDirectory}/Resistance_prediction.html',index=False)


dfoutout.to_csv(f'{args.inputDirectory}/Resistance_prediction_details.csv',index=False)
dfoutout.to_html(f'{args.inputDirectory}/Resistance_prediction_details.html',index=False)


#print (os.path.abspath(f'{args.inputDirectory}/Finish.html'))
#fpath=os.path.abspath(f'{args.inputDirectory}/Finish.html')
#webbrowser.open_new_tab('file://'+fpath)
######################



#讀取我的表
#copy df_fin DataFrame
df_me=df_fin.copy()

#讀取uniBarcode資訊
# get  uniBarcode list
uniBarcode=df_me['barcode'].drop_duplicates()


df_fin=pd.DataFrame(columns=column_name)
tmpdict=dict()
for Barcode in uniBarcode:
    tmpdict=dict()
    df_tmp=df_me[df_me['barcode']==Barcode]
    tmpSpecies=df_tmp.iloc[0,1]
    uniPhenotype=df_tmp['Phenotype'].drop_duplicates().to_list()


    tmpdict['Barcode']=Barcode
    tmpdict['Species name']=tmpSpecies

    for colName in drug_list:
        try :
            if colName in uniPhenotype:
                tmpdict[colName]='R'
            elif colName in dict_name[tmpSpecies]:
                tmpdict[colName]='S'
        except KeyError  as e:
            if colName in uniPhenotype:
                tmpdict[colName]='R'
            else :
                tmpdict[colName]='S'
#    for phenotype in uniPhenotype:
#        if phenotype in column_name:

    df_tmp=pd.DataFrame(data=[tmpdict])
    df_fin=pd.concat([df_fin,df_tmp])
#df_fin.fillna('R',inplace=True)
df_fin.to_csv(f'{args.inputDirectory}/Resistance_4comparison.csv',index=False)

