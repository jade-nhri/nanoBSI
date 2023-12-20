#!/usr/bin/env python3
import os,sys
import subprocess
import datetime
import pandas as pd

inlist=sys.argv[1]
infa=sys.argv[2]
phenfile='/opt/resfinder_db/phenotypes.txt'
classfile='/opt/resfinder_db/antibiotic_classes.txt'
cur_datetime=datetime.datetime.now()
for_datetime=cur_datetime.strftime("%Y%m%d_%H%M%S")
print (for_datetime)

comm='cp {0} {1}'.format(phenfile,phenfile.replace('.txt','_'+for_datetime+'.txt'))
print (comm)
subprocess.getoutput(comm)
comm='cp {0} {1}'.format(classfile,classfile.replace('.txt','_'+for_datetime+'.txt'))
print (comm)
subprocess.getoutput(comm)


d=dict()
f=open(infa)
def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

with f as fp:
    for name, seq in read_fasta(fp):
        name=name.replace('>','')
        d[name]=seq
f.close()
print (d)

listd=dict()
def read_infile(fp):
    for line in fp:
        line = line.rstrip()
        text=line.split('\t')
        print (text)
        listd[text[0]]=text[1:]
    return (listd)

f=open(inlist)
with f as fp:
    ld=read_infile(fp)
    print (ld)
f.close()


df=pd.read_table(phenfile,keep_default_na=False)
print (df)

for i in ld.keys():
    process=ld[i]
    cond=len(process)
    if cond==5:
        print (process[0])
        fsa=os.path.join('/opt/resfinder_db',process[0])
        if os.path.exists(fsa) and i in d.keys():
            comm='cp {0} {1}'.format(fsa,fsa.replace('.fsa','_'+for_datetime+'.fsa'))
            subprocess.getoutput(comm)
            fw=open(fsa,'a')
            fw.write('>'+i+'\n')
            fw.write(d[i]+'\n')
            fw.close()
            new_row=pd.DataFrame({'Gene_accession no.':[i],'Class':[process[1]],'Phenotype':[process[2]],'PMID':[process[3]],'Mechanism of resistance':[''],'Notes':[process[4]],'Required_gene':['']})
            #new_row=pd.DataFrame({'Gene_accession no.':[i],'Class':[process[1]],'Phenotype':[process[2]],'PMID':[process[3]],'Notes':[process[4]]})
            #print (new_row)
            df=pd.concat([df,new_row],ignore_index=True)
            #print (df)
    
    if cond==2:
        #print (i)
        stri=i+'_'
        #print (stri)
        s1=df['Gene_accession no.']
        #print (s1)
        a=s1[s1.str.contains(stri,regex=True)].index.tolist()
        #print (a)
        if os.path.exists(classfile):
            comm="grep '{0}' {1} -n".format(process[0],classfile)
            #print (comm)
            stdout=subprocess.getoutput(comm)
            #print (stdout)
            linen=stdout.split(':')[0]
            #print (linen)
            comm="sed -i '{0}s/$/\t{1}/' {2}".format(linen,process[1],classfile)
            #print (comm)
            subprocess.getoutput(comm)


        for j in a:
            #print (j)
            tmp=df.loc[j].at['Phenotype']
            df.at[j,'Phenotype']=process[1]+', '+tmp
            #print (df.loc[j])


    if cond==1:
        #print (process[0])
        #print (i)
        a=df[df['Gene_accession no.']==i].index.tolist()
        #print (a)
        tmp=df.loc[a[0]].at['Phenotype']
        df.at[a[0],'Phenotype']=tmp+', '+process[0]
        #print (df.loc[a[0]])


df.to_csv(phenfile,sep='\t',index=False)
