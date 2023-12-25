#!/usr/bin/env python3
import os,sys
import subprocess
import multiprocessing as mp
import operator
import pandas as pd
from ete3 import NCBITaxa
ncbi = NCBITaxa()
import warnings
warnings.filterwarnings("ignore")
warnings.simplefilter("ignore")
database_path = ncbi.dbfile
if not os.path.exists(database_path):
    ncbi.update_taxonomy_database()

threads=8
multip=24

indir=sys.argv[1]
indir=os.path.abspath(indir)
print (indir)
outdir=sys.argv[2]
outdir=os.path.abspath(outdir)
if os.path.exists(outdir):
    comm='rm {0} -rf'.format(outdir)
    subprocess.getoutput(comm)
if not os.path.exists(outdir):
    os.makedirs(outdir)



def catfiles(bc):
    #print ('Run concatenating files')
    comm='cat {0}/{2}/*.fastq.gz > {1}/{2}.fastq.gz'.format(indir,outdir,bc)
    #print (comm)
    subprocess.getoutput(comm)
    comm='cat {0}/{1}.fastq.gz | seqkit stats'.format(outdir,bc)
    stdout=subprocess.getoutput(comm)
    print ('{0}\n{1}'.format(bc,stdout))
    runcen(bc)


def runcen(bc):
    predictspe=''
    psped={'Campylobacter coli':'Campylobacter','Campylobacter jejuni':'Campylobacter','Enterococcus faecalis':'Enterococcus faecalis','Enterococcus faecium':'Enterococcus faecium','Escherichia coli':'Escherichia coli','Helicobacter pylori':'Helicobacter pylori','Klebsiella pneumoniae':'Klebsiella','Klebsiella quasipneumoniae':'Klebsiella','Klebsiella aerogenes':'Klebsiella','Klebsiella oxytoca':'Klebsiella','Klebsiella variicola':'Klebsiella','Mycobacterium tuberculosis':'Mycobacterium tuberculosis','Neisseria gonorrhoeae':'Neisseria gonorrhoeae','Plasmodium falciparum':'Plasmodium falciparum','Salmonella enterica':'Salmonella','Staphylococcus aureus':'Staphylococcus aureus','Staphylococcus capitis':'Staphylococcus aureus','Staphylococcus epidermidis':'Staphylococcus aureus','Staphylococcus caprae':'Staphylococcus aureus','Staphylococcus haemolyticus':'Staphylococcus aureus','Staphylococcus hominis':'Staphylococcus aureus','Staphylococcus cohnii':'Staphylococcus aureus','Staphylococcus pettenkoferi':'Staphylococcus aureus','Staphylococcus argenteus':'Staphylococcus aureus'}
    intrinsic_res={'Enterococcus faecium':['Penicillin'],'Enterococcus lactis':['Penicillin'],'Klebsiella aerogenes':['Cefoxitin'],'Bacteroides fragilis':['Penicillin'],'Bacteroides thetaiotaomicron':['Penicillin'],'Serratia marcescens':['Ampicillin','Cefazolin','Cefoxitin'],'Enterobacter cloacae complex sp.':['Ampicillin','Cefoxitin'],'Enterobacter asburiae':['Ampicillin','Cefoxitin'],'Enterobacter cloacae':['Cefoxitin'],'Enterobacter hormaechei':['Ampicillin','Cefoxitin'],'Ralstonia mannitolilytica':['Amikacin','Gentamicin'],'Leptotrichia wadei':['Amikacin'],'Staphylococcus pettenkoferi':['Pencillin','Oxacillin']}

    comm='centrifuge -x /opt/cbsidb -U {0}/{1}.fastq.gz -S {0}/{1}_cenout.txt -p {2} --min-hitlen 60'.format(outdir,bc,threads)
    #print (comm)
    subprocess.getoutput(comm)
    comm='getspe.py {0}/{1}_cenout.txt'.format(outdir,bc)
    #print (comm)
    stdout=subprocess.getoutput(comm)
    #print (stdout)
    with open(os.path.join(outdir,'{0}_speout.txt'.format(bc)),'w') as f:
        print (stdout, file=f)
    predictspe=stdout.split(' with')[0]
    #print (predictspe)
    if predictspe in intrinsic_res.keys():
        print(predictspe)

    fw=open(os.path.join(outdir,'Species_identification.txt'),'a')
    fw.write('{0}\t{1}\n'.format(bc,predictspe))
    fw.close()
    print ('{0} is {1}'.format(bc,predictspe))
    spe=predictspe.replace(' ','_')
    if predictspe in psped.keys():
        comm='python3 -m resfinder -k /opt/kma/kma -o {0}/{1}_res -l 0.6 -t 0.9 --acquired -ifq {0}/{1}.fastq.gz --nanopore -c --species "{2}"'.format(outdir,bc,psped[predictspe])
    else:
        comm='python3 -m resfinder -k /opt/kma/kma -o {0}/{1}_res -l 0.6 -t 0.9 --acquired -ifq {0}/{1}.fastq.gz --nanopore'.format(outdir,bc)
    #print (comm)
    stdout=subprocess.getoutput(comm)

    pfile='{0}/{1}_res/PointFinder_results.txt'.format(outdir,bc)
    #print (pfile)
    rfile='{0}/{1}_res/ResFinder_results_tab.txt'.format(outdir,bc)
    #print (rfile)

    if os.path.exists(pfile):
        #print (bc)
        pdf=pd.read_table(pfile)
        #print (pdf)
        pdf.rename(columns={'Mutation':'Contig','Resistance':'Phenotype'},inplace=True)
        pdf=pdf[['Contig','Phenotype']]
        pdf['Resistance gene']='Point'
        pdf['Accession no.']='Point'
        #print (pdf)
        rdf=pd.read_table(rfile)
        #print (rdf)
        rdfnew=pd.concat((rdf,pdf),axis=0)
        #print (rdfnew)
        subprocess.getoutput('mv {0} {1}'.format(rfile,rfile.replace('.txt','_old.txt')))
        rdfnew.to_csv(rfile,index=False,sep='\t')
    if predictspe in intrinsic_res.keys():
        ind={'Phenotype':intrinsic_res[predictspe]}
        indf=pd.DataFrame(ind)
        indf['Resistance gene']='Intrinsic'
        indf['Accession no.']='Intrinsic'
        rdf=pd.read_table(rfile)
        rdfnew=pd.concat((rdf,indf),axis=0)
        subprocess.getoutput('mv {0} {1}'.format(rfile,rfile.replace('.txt','_iold.txt')))
        rdfnew.to_csv(rfile,index=False,sep='\t')

def getdirsize(indir):
    total_size=0
    for dirpath, dirnames, filenames in os.walk(indir):
        for f in filenames:
            fp=os.path.join(dirpath,f)
            if not os.path.islink(fp):
                total_size+=os.path.getsize(fp)
    return (total_size)


os.chdir(indir)
mydirs=[x for x in os.listdir(indir) if os.path.isdir(x) and 'barcode' in x]
print (mydirs)

dsize=dict()
for i in mydirs:
    dsize[i]=getdirsize(i)
sortedbc=sorted(dsize.items(),key=operator.itemgetter(1),reverse=True)
print (sortedbc)

def main():
    po=mp.Pool(multip)
    for i,j in sortedbc:
        #print (i)
        if j>1000000:
            po.apply_async(catfiles,[i])
            #po.apply_async(runresfinder,[i])

    po.close()
    po.join()

    os.chdir(outdir)
    comm='getres_reads.py'    #Input folder *_res
    print (comm)
    subprocess.getoutput(comm)


    comm='getresults_reads.py -i {0}'.format(outdir)
    print (comm)
    subprocess.getoutput(comm)

if __name__=='__main__':
    main()
