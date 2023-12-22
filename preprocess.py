#!/usr/bin/env python3
import os, sys
import subprocess

print ("To download Centrifuge index...") 
comm='wget -O /opt/cbsidb.zip "https://www.dropbox.com/scl/fi/g65ndnvw1gf927n3idr2h/cbsidb.zip?rlkey=nr1q6d0b6200ir2b84mhg04ir&dl=0"'
subprocess.getoutput(comm)
comm='unzip /opt/cbsidb.zip -d /opt/' 
subprocess.getoutput(comm)

print ("To download database sequences...") 
comm='wget -O /opt/csequences.fna "https://www.dropbox.com/scl/fi/8kft1bc0xtlu05spthkie/csequences.fna?rlkey=p09mkkuzxtgcnrn2nj0o441se&dl=0"' 
subprocess.getoutput(comm)

print ("To download plasmid sequences...") 
comm='wget -O /opt/plasmid.mmi "https://www.dropbox.com/scl/fi/97vs73gjt3pngx8hf55jt/plasmid.mmi?rlkey=dy1onif2gr4wozmiz7c0i5tp3&dl=0"' 
subprocess.getoutput(comm)

print ("To update ResFinder database...") 
comm='updateres.py /opt/nanoBSI/newres.txt /opt/nanoBSI/newres.fasta'
subprocess.getoutput(comm)

print ("To generate gene2note.json...") 
comm='getnotes.py /opt/resfinder_db/phenotypes.txt /opt/gene2note.json'
subprocess.getoutput(comm)

print ("Preprocess was done.")
