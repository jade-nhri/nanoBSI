#!/usr/bin/env python3
import os,sys
import pandas as pd
import subprocess

indir=os.path.abspath(sys.argv[1])
cwd=os.path.abspath(os.getcwd())
#print (cwd)

os.chdir(indir)
oldspeid=pd.read_table('Species_identification.txt',names=['barcode','Species'])
#print (oldspeid)

osped=dict()
for i,j in zip(oldspeid['barcode'].values,oldspeid['Species']):
    osped[i]=j
#print (osped)

mydirs=[]
for i in osped.keys():
    mydirs.append(i+'_spe')
#print (mydirs)

fw=open(os.path.join(indir,'Species_identification_wconf.txt'),'w')
for i in sorted(mydirs):
    os.chdir(i)
    bc=i.replace('_spe','')
    if not os.path.exists('cluster.txt') and not os.path.exists('cluster1.txt'):
        #print ('{0}\t{1}\t{2}'.format(i,osped[bc],'L'))
        fw.write('{0}\t{1}\t{2}\n'.format(i,osped[bc],'L'))
        os.chdir('..')
        continue
    myfiles=[x for x in os.listdir() if 'cluster' in x and '.txt' in x and  'old' not in x]
    #print (myfiles)
    for j in sorted(myfiles):
        csize=0
        rsize=0
        comm='grep "^Cluster" {0}'.format(j)
        #print (comm)
        stdout=subprocess.getoutput(comm)
        #print (stdout)
        nspe=stdout.split(': ')[-1]
        #print (j)
        if j=='cluster.txt':
            csize=os.path.getsize('cluster.fasta')
            #print (csize)
            if nspe==osped[bc]:
                #print ('{0}\t{1}\t{2}'.format(i,nspe,'H'))
                rsize=os.path.getsize('../{0}_Dref.fasta'.format(bc))
                #print (rsize)
                if csize>=rsize*0.9:
                    #print ('rule1')
                    fw.write('{0}\t{1}\t{2}\n'.format(i,nspe,'HC'))
                    continue
                if csize<rsize*0.7:
                    fw.write('{0}\t{1}\t{2}\n'.format(i,nspe,'HIC'))
                    continue
                else:
                    #print ('rule3')
                    fw.write('{0}\t{1}\t{2}\n'.format(i,nspe,'H'))
            else:
                #print ('{0}\t{1}\t{2}'.format(i,nspe,'M'))
                fw.write('{0}\t{1}\t{2}\n'.format(i,nspe,'M'))
        else:
            tmp=j.replace('luster','').replace('.txt','')
            if nspe==osped[bc]:
                #print ('{0}\t{1}\t{2}'.format(i,nspe,'H'))
                fw.write('{0}_{3}\t{1}\t{2}\n'.format(i,nspe,'H',tmp))
            else:
                #print ('{0}\t{1}\t{2}'.format(i,nspe,'M'))
                fw.write('{0}_{3}\t{1}\t{2}\n'.format(i,nspe,'M',tmp))

    os.chdir('..')
os.chdir(cwd)
fw.close()


