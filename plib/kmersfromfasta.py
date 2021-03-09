# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 17:25:21 2021

@author: thomas.hauling
"""

#kmersfromfasta

#import pandas as pd
import glob

# get all fasta files
def getfilenames():
    fastatypes = ('*.fa', '*.fasta', '*.fna')
    foundfasta = []
    for f in fastatypes:
        foundfasta.extend(glob.glob(f))
    return foundfasta

#files = getfilenames()
#if len(getfilenames()) == 0: print('no fasta file(s) found')


def makekmers(kmerlen, fastaoutput):
    '''processes fasta files into k-mers'''
    fastainput = getfilenames()
    allf = []
    for f in fastainput:
        with open(f) as allfasta:
            #next(allfasta)  ## use this to remove FIRST line, ie fasta header
            for line in allfasta:
                allf.append(line.rstrip())
            
            allfstr = '\n'.join(allf)
           # allfstrold = allfstr
### remove ALL line breaks '\n' and carriage returns '\r' from sequence portion BEFORE using split
            allfstr = ''.join(allf)

#### make kmers ### run this overnight
#window = 40
#kmers = [allfstr[i:i+window] for i in range(len(allfstr)-2)] ## why -2?
    kmers = (allfstr[i:i+kmerlen] for i in range(len(allfstr)-2)) ## why -2?
#kmers =[i for i in kmers if len(i) == 30]
#kmers =[i for i in kmers if len(i) == window]
    kmers =(i for i in kmers if len(i) == kmerlen)

    purekmers=[]
#s=['>', ' ', '.', 'N']    ### add expression to test if any char are NOT in basel=['A','T','C','G']
    for k in kmers:
      if all(c.isalpha() and c.isupper() for c in k):
        purekmers.append(k)
    else:
        pass
        
        # if '>' in k:
        #     pass
        # elif ' ' in k:
        #     pass
        # elif '.' in k:
        #     pass
        # elif 'N' in k:
        #     pass
        # elif any(c.islower() for c in k):
        #     pass
        # elif any(c.isdigit() for c in k):
        #     pass
        # elif any(c.isspace() for c in k):
        #     pass
        # elif any(c.isnumeric() for c in k):
        #     pass
        # elif all(c.isalpha() for c in k):
        #     purekmers.append(k)
        # # return
        # else:
        #     purekmers.append(k)
            
            
    kmers_uni = list(dict.fromkeys(purekmers)) #reduces kmer list by 75ish % 
##118 Mil vs 437 Mil

## save all kmers
    fastafile = open(fastaoutput, "w")
    counter = 0
#for fa in purekmers:
    for fa in kmers_uni:
#for fa in short_kmers:
        fastafile.writelines(['>' + str(counter) + '\n', fa + '\n'])
        counter = counter +1
    fastafile.close()
    
#makekmers(fastanames, 32, 'kmers32.fa')
