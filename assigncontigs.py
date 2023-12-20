#!/usr/bin/env python3
import os,sys
import subprocess
import pandas as pd
import numpy as np

assembly=sys.argv[1]
reffile=sys.argv[2]
outdir=sys.argv[3]
if not os.path.exists(outdir):
    os.mkdir(outdir)


d=dict()
f=open(assembly)
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
        readid=name.replace('>','')
        d[readid]=seq
f.close()


assemblyfile=os.path.split(reffile)[-1]
print (assemblyfile) #Dref
outpaf=os.path.join(outdir,assemblyfile.replace('_Dref.fasta','_temp.paf'))
outpaf1=os.path.join(outdir,assemblyfile.replace('_Dref.fasta','_plasmid_temp.paf'))
comm='minimap2 {0} {1} > {2}'.format(reffile,assembly,outpaf)
subprocess.getoutput(comm)

###20230609
speinfo=pd.read_table(reffile.replace('_Dref.fasta','_cenout2.txt'),names=['taxID','Tname','Species'])
print (speinfo)

def assignspe(acclist):
    specount=dict()
    a2s=dict()
    for i,j in zip(speinfo['Tname'].values,speinfo['Species'].values):
        a2s[i]=j
    if len(acclist)==1:
        speout=a2s[acclist[0]]
    else:
        uniacc=np.unique(acclist)
        #print (uniacc)
        count=1
        for i in uniacc:
            tmpc=acclist.count(i)
            speout=a2s[i]
            if count==1:
                maxc=tmpc
            if count>1 and tmpc>maxc:
                speout=a2s[i]
            specount[a2s[i]]=tmpc
    return(speout)
    


#info=pd.read_table(assembly.replace('.fasta','_info.txt'),names=['seq_name','length','cov.','circ. ','repeat','mult.','alt_group','graph_path'])
info=pd.read_table(assembly.replace('.fasta','_info.txt'))
print (info)
circontigs=info[info['circ.']=='Y']['#seq_name'].values
print (circontigs)

df=pd.read_table(outpaf,names=['Qname','Qlen','Qstart','Qend','strand','Tname','Tlen','Tstart','Tend','Nmatch','Alen','MQ1','MQ2','MQ3','MQ4','MQ5','MQ6','MQ7'])
print (df)

ref=np.unique(df['Tname'].values)
print (ref)

if len(ref)==1:
    dfset1=df.drop_duplicates('Qname')
    contigids=dfset1['Qname'].values
    #print (contigids)
    fw=open(os.path.join(outdir,'cluster.fasta'),'w')
    fw1=open(os.path.join(outdir,'cluster.txt'),'w')
    fw1.write('Cluster: {0}\n'.format(assignspe(ref)))
    fw1.close()
    for i in contigids:
        if i in circontigs:
            fw.write('>{0} Circular\n'.format(i))
        else:
            fw.write('>{0}\n'.format(i))
        fw.write('{0}\n'.format(d[i]))
    fw.close()
    fw=open(os.path.join(outdir,'others.fasta'),'w')
    for i in d.keys():
        if i not in contigids:
            if i in circontigs:
                fw.write('>{0} Circular\n'.format(i))
            else:
                fw.write('>{0}\n'.format(i))
            fw.write('{0}\n'.format(d[i]))
    fw.close()

    comm='minimap2 /opt/plasmid.mmi {0} > {1}'.format(os.path.join(outdir,'others.fasta'),outpaf1)
    print (comm)
    subprocess.getoutput(comm)
    df1=pd.read_table(outpaf1,names=['Qname','Qlen','Qstart','Qend','strand','Tname','Tlen','Tstart','Tend','Nmatch','Alen','MQ1','MQ2','MQ3','MQ4','MQ5','MQ6','MQ7'])
    print (df1)
    print (len(df1))
    if len(df1)>0:
        plasids=df1.drop_duplicates('Qname')['Qname'].values
        #print (plasids)
        contigsnew=np.append(contigids,plasids)
        #print (contigsnew)
        comm='mv {0} {1}'.format(os.path.join(outdir,'cluster.fasta'),os.path.join(outdir,'cluster_old.fasta'))
        subprocess.getoutput(comm)
        comm='mv {0} {1}'.format(os.path.join(outdir,'others.fasta'),os.path.join(outdir,'others_old.fasta'))
        subprocess.getoutput(comm)
        fw=open(os.path.join(outdir,'cluster.fasta'),'w')
        for i in contigsnew:
            if i in circontigs:
                fw.write('>{0} Circular\n'.format(i))
            else:
                fw.write('>{0}\n'.format(i))
            fw.write('{0}\n'.format(d[i]))
        fw.close()
        fw=open(os.path.join(outdir,'others.fasta'),'w')
        for i in d.keys():
            if i not in contigsnew:
                if i in circontigs:
                    fw.write('>{0} Circular\n'.format(i))
                else:
                    fw.write('>{0}\n'.format(i))
                fw.write('{0}\n'.format(d[i]))
        fw.close()



else:
    meancov=[]
    refs=ref.tolist()
    print (refs)
    cluster=dict()
    for i in range(0,len(ref)):
        cluster[i]=dict()
        meancov.append(0)
    #print (cluster)
    for i,j,k,l in zip(info['#seq_name'],info['length'],info['cov.'],info['circ.']):
        #print ('{0}\t{1}\t{2}\t{3}'.format(i,j,k,l))
        dfset=df[df['Qname']==i]
        #print (dfset)
        dfset1=dfset.drop_duplicates('Qname')
        #print (dfset1)
        if len(dfset1)==0:
            continue
        else:
            index=refs.index(dfset1['Tname'].values)
            #print (index)
            if j>=10000 and l=='N':
                for s in range(0,len(ref)):
                    if index==s:
                        if not cluster[s].keys():
                            cluster[s]['contig']=[i]
                            cluster[s]['cov']=[k]
                            cluster[s]['refs']=[refs[s]]
                        else:
                            meancov[s]=np.mean(cluster[s]['cov'])
                            #print (meancov)
                            cluster[s]['contig'].append(i)
                            cluster[s]['cov'].append(k)
                            cluster[s]['refs'].append(refs[s])
            if j<10000 and l=='N':
                for s in range(0,len(ref)):
                    if k<=1.6*meancov[s] and k>0.6*meancov[s]:
                        meancov[s]=np.mean(cluster[s]['cov'])
                        #print (meancov)
                        cluster[s]['contig'].append(i)
                        cluster[s]['cov'].append(k)
                        cluster[s]['refs'].append(refs[s])
                    if index==s and abs(k-meancov[s])<=10 and i not in cluster[s]['contig']:
                        cluster[s]['contig'].append(i)
                        cluster[s]['cov'].append(k)
                        cluster[s]['refs'].append(refs[s])
    #print (cluster)
   

    contigsnew=[]
    for i in range(0,len(ref)):
        cluster[i]['Ref']=assignspe(cluster[i]['refs'])
        fw=open(os.path.join(outdir,'cluster{0}.fasta'.format(i+1)),'w')
        fw1=open(os.path.join(outdir,'cluster{0}.txt'.format(i+1)),'w')
        fw1.write('Cluster{0}: {1}\n'.format(i+1,cluster[i]['Ref']))
        for j,k,l in zip(cluster[i]['contig'],cluster[i]['cov'],cluster[i]['refs']):
            contigsnew.append(j)
            if j in circontigs:
                fw.write('>{0} Circular\n'.format(j))
                fw1.write('{0}\t{1}\t{2}\n'.format(j,k,l))
            else:
                fw.write('>{0}\n'.format(j))
                fw1.write('{0}\t{1}\t{2}\n'.format(j,k,l))
            fw.write('{0}\n'.format(d[j]))
        fw.close()
        fw1.close()
    #print (len(contigsnew))
    fw=open(os.path.join(outdir,'others.fasta'),'w')
    for i in d.keys():
        if i not in contigsnew:
            if i in circontigs:
                fw.write('>{0} Circular\n'.format(i))
            else:
                fw.write('>{0}\n'.format(i))
            fw.write('{0}\n'.format(d[i]))
    fw.close()
    print (cluster)

    comm='minimap2 /opt/plasmid.mmi {0} > {1}'.format(os.path.join(outdir,'others.fasta'),outpaf1)
    print (comm)
    subprocess.getoutput(comm)
    df1=pd.read_table(outpaf1,names=['Qname','Qlen','Qstart','Qend','strand','Tname','Tlen','Tstart','Tend','Nmatch','Alen','MQ1','MQ2','MQ3','MQ4','MQ5','MQ6','MQ7'])
    print (df1)
    print (len(df1))
    contignew1=[]
    if len(df1)>0:
        plasids=df1.drop_duplicates('Qname')['Qname'].values
        print (plasids)
        print (len(plasids))
        comm='mv {0} {1}'.format(os.path.join(outdir,'cluster1.fasta'),os.path.join(outdir,'cluster1_old.fasta'))
        subprocess.getoutput(comm)
        comm='mv {0} {1}'.format(os.path.join(outdir,'cluster2.fasta'),os.path.join(outdir,'cluster2_old.fasta'))
        subprocess.getoutput(comm)
        comm='mv {0} {1}'.format(os.path.join(outdir,'others.fasta'),os.path.join(outdir,'others_old.fasta'))
        subprocess.getoutput(comm)
        for i in range(0,len(ref)):
            fw=open(os.path.join(outdir,'cluster{0}.fasta'.format(i+1)),'w')
            for j in cluster[i]['contig']:
                contignew1.append(j)
                if j in circontigs:
                    fw.write('>{0} Circular\n'.format(j))
                else:
                    fw.write('>{0}\n'.format(j))
                fw.write('{0}\n'.format(d[j]))
            for k in plasids:
                contignew1.append(k)
                if k in circontigs:
                    fw.write('>{0} Circular\n'.format(k))
                else:
                    fw.write('>{0}\n'.format(k))
                fw.write('{0}\n'.format(d[k]))
            for l in circontigs:    #20230718
                if l not in cluster[i]['contig'] and l not in plasids:
                    contignew1.append(l)
                    fw.write('>{0} Circular\n'.format(l))
                    fw.write('{0}\n'.format(d[l]))
            fw.close()
        fw=open(os.path.join(outdir,'others.fasta'),'w')
        for i in d.keys():
            if i not in contignew1:
                if i in circontigs:
                    fw.write('>{0} Circular\n'.format(i))
                else:
                    fw.write('>{0}\n'.format(i))
                fw.write('{0}\n'.format(d[i]))
        fw.close()

    
