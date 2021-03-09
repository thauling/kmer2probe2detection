# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 16:52:42 2021

@author: thomas.hauling
"""


################################ dictionary created ########################
## save and recreate gene library
#import sys
import pandas as pd

#sys.path.append(r'C:\Users\Thomas\Documents\Python_and_R\plib') 

## load already created refseq:symbol dict d:
def ref2sym_chunk(alignout):     
    ddf=pd.read_csv('d.csv', sep='\t', header=0)
#print(ddf == ddf2)
    d = pd.Series(ddf.symbol.values,index=ddf.refseq).to_dict()

##### function to find matching gene symbols for mmseq2 output chunks  
    symbolout=[]

    def ref2sym(alignchunk):
        for mm in alignchunk[1]: 
            try: 
                symbolout.append(''.join(d[mm]))
        
            except:
            #pass #'nan' appended only at start of symbolout - why? 
                symbolout.append('na')
            #print(symbolout)
            #return symbolout #not required?

#    def ref2sym_chunk(alignout): 
    chunks = pd.read_csv(alignout, sep='\t', header=None, chunksize=1000000)
#print(chunks)
    symbolout=[]

    for c in chunks:
        print(c[1])
        ref2sym(c)

### write only symbols to file (to be save)
#f=open('/home/octopod/symbolout_start.txt','w') #create new text file, will be overwritten ('w' option), use 'a+' to create and appen
#f=open('/media/octopod/2TBext4/kmers40_symbols.txt','a+') # create text file and APPEND with data if already exisiting
    f=open(alignout + '_symbols.txt','w') 
    for ele in symbolout:
        f.write(ele+'\n')
    f.close()