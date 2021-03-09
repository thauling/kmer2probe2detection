# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 16:52:42 2021

@author: thomas.hauling
"""


import pandas as pd
import numpy as np
import sqlite3 
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC
from seqfold import dg

## load already created refseq:symbol dict d:
def ref2sym_chunk(alignout):     
    ddf=pd.read_csv('d.csv', sep='\t', header=0)
#print(ddf == ddf2)
    d = pd.Series(ddf.symbol.values,index=ddf.refseq).to_dict()

##### function to find matching gene symbols for mmseq2 output chunks  
    symbolout=[]

    def ref2sym(alignchunk):
        '''finds gene symbol for every RefSeq Nr in alignment-out (mmseqs) file'''
        for mm in alignchunk[1]: 
            try: 
                symbolout.append(''.join(d[mm]))
        
            except:
            #pass #'nan' appended only at start of symbolout - why? 
                symbolout.append('na')
    chunks = pd.read_csv(alignout, sep='\t', header=None, chunksize=1000000)
    symbolout=[]

    for c in chunks:
        ref2sym(c)
### write only symbols to file (to be save)
    f=open(alignout + '_symbols.txt','w') 
    for ele in symbolout:
        f.write(ele+'\n')
    f.close()
    
def out2db(mmseqout, symbols, dbout):
    '''removes repetetive motifs, calculates Tm, %GC, dG (currently set to 0) and transfers results to db'''
    con = sqlite3.connect(dbout)
## write txt to db

    chsize = 100000 ##chsize 1000000000 1617 sec chsize=1000000 1664 sec (4242 sec with modin) , chsize=10000 3398 sec
    rc1 = pd.read_csv (mmseqout, sep = '\t', header=None, chunksize=chsize)
    rc2 = pd.read_csv (symbols, sep = '\t', header=None, chunksize=chsize)

#def out2db(mmseqout, symbols, db):
#for df1, df2 in zip(rc1,rc2):
    for df1, df2 in zip(rc1,rc2):
        #print(df1, df2)
            df3=pd.concat([df1,df2], axis=1, ignore_index=True)
        #print(df3)
        ### sort df by RefSeq Acc. and position (start or end) ## first step to calculate distances between kmers/ RefSeq Acc.
            df3.sort_values(by=[1,3], ascending=[True, True],inplace=True) ## first SORT according to kmer id 

            df3=df3.drop_duplicates() ## QC - in case there are entire duplicate rows
        #df3b=df3 ## might be better to drop duplicates AFTER overlap determination, anyways 
            df3 = df3[df3[7].str.contains('CCCC') != True] 
            df3 = df3[df3[7].str.contains('GGGGG') != True] 
            df3 = df3[df3[7].str.contains('AAAAAA') != True] 
            df3 = df3[df3[7].str.contains('TTTTTT') != True] 
            df3 = pd.DataFrame({'Symbol':df3[8],'RefSeq':df3[1],'Start':df3[3],'End':df3[4],'Bitscore':df3[6],'Sequence':df3[7]})
            df3['Tm'] = df3['Sequence'].apply(mt.Tm_NN)
            df3['GC'] = df3['Sequence'].apply(GC)
            #df3['DeltaG'] = df3['Sequence'].apply(dg) # takes FOREVER
            df3['DeltaG'] = pd.DataFrame(np.zeros(len(df3))) #mock dG!
            #print(df3)
            df3.to_sql('mmseqt', con, if_exists='append', index = False)
    con.commit()  ## save changes to db
    con.close()   