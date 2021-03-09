# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 08:15:06 2021

@author: Thomas
"""

# read kmers from SQlite DB and other formats
# for probe building etc
import os
import pandas as pd
import numpy as np
import sqlite3 
from seqfold import dg
#import sqlalchemy

cwd = os.getcwd()
os.chdir(r'C:\Users\Thomas\Documents\Python_and_R')
#dbfile='test2_db.db'
def readdb(dbfile):
    con = sqlite3.connect(dbfile)
    df = pd.read_sql_query("SELECT * from mmseqt", con) #.groupby('symbol').count()['RNA seq']
    con.close()
    print(df.head())
    print(df.tail())
    return df
    
def readparquet(pqfile):
    df = pd.read_parquet(path=pqfile, engine='pyarrow')
    print(df.head())
    print(df.tail())
    return df

## query database, retrieve all entries of symbol=x, bitscore=y
dbfile = 'dbout'
symbolfile = 'mysymbols.txt' 
symbol = 'symbol'   
gene = 'Meis2'
genes = ('Meis2', 'Robo2') # this creates a tuple
score = 'bitscore'
# returns df of db rows for specified genes
def getkmers(dbfile, symbolfile): ##symbols, scores have to be tuples
    allrows = pd.DataFrame()
    try:
        con = sqlite3.connect(dbfile)
    except:
        print('could not connect to db')
    cur = con.cursor()
    f = open(symbolfile, 'rt')
    genes = tuple(f.read().split())
    
    for g in genes:
        cur.execute("SELECT * FROM mmseqt WHERE %s=?" % (symbol,), (g,)) #this works!
        #cur.execute("SELECT * FROM mmseqt WHERE %s=? AND WHERE %s=73" % (symbol, score), (g,)) #
        rows = pd.DataFrame(cur.fetchall()) ## reads db as list
        #print(rows)
        allrows = allrows.append(rows) # why is this not growing??
    #df = pd.read_sql_query("SELECT * FROM mmseqt", con) #, symbols[0]) # reads db as df
    con.close()
    return allrows 
    #seems to work! #QC: remove duplicate kmers if present

myrows = getkmers(dbfile,symbolfile) # test func
def filterkmers(myrows, minscore=73, mintm=55, maxtm=65, mingc=45, maxgc=55, mindg=-0.5, maxdg=0.5):
    selection = myrows[(myrows[4] >= minscore) & (myrows[6].between(mintm,maxtm)) & (myrows[7].between(mingc,maxgc))]#.any())]
    ## calulate dG on selection
    selection['deltaG'] = selection[5].apply(dg) ##assumes 37degC I guess
    selection = selection[selection['deltaG'].between(mindg,maxdg)]
    return selection

myselection = filterkmers(myrows, mindg=-10, maxdg=10) ## this works but setting with copy warning

def getkmerspergene(myselection, genenumbers, kmerdistance=20):
    nonoverlapAll = pd.DataFrame()
   # allprobes= pd.DataFrame()
   # numbersdict = pd.Series(numbersfile[1].values,index=numbersfile[0]).to_dict()
    refseq=list(np.unique(myselection[1])) ##make non-redundant list of all RefSeq Accessions
    #rs=refseq[0]
    ## apply simple distance filter, has to be done PER RefSeq! check previous script
    #kmerdistance=20 ##fro debugging
    for rs in refseq: #[0]:
             myrefseq=myselection[myselection[1] == rs] 
             y = np.array(myrefseq)
             z = myrefseq[4][0::kmerdistance] ##remove [0::distance] to select all
             posdict = dict(zip(z,y))
             nonoverlapRs = pd.DataFrame(posdict.values()) ##converts posdict to dataframe 
             nonoverlapAll = nonoverlapAll.append(nonoverlapRs)     
             grouped = nonoverlapAll.groupby(0)[5].apply(list) #.agg([5])
             nonoverlapAll.columns=['Gene','RefSeq','Start','End','Bitscore','Sequence','Tm','GC','DeltaG']
    return grouped
## now get n kmers per gene, starting from first kmer
   # for no in set(nonoverlapAll[0]): # use this to iterate over genes
def getnumberedkmers(groupednonoverlap, numbersfile):
    allprobes= pd.DataFrame()
    numbersdict = pd.Series(numbersfile[1].values,index=numbersfile[0]).to_dict()
    for no in groupednonoverlap.items(): # or iterate of items
        gene = (no[0])
        try:
            probenumber = numbersdict.get(gene)
            print(probenumber)
            print(no[1][:probenumber]) #voila! 
            probeseqs = (no[1][:probenumber])
            probes =([gene], [probenumber], [probeseqs])
            probesdf = pd.DataFrame(probes).T
            print(probes)
            allprobes = allprobes.append(probesdf) ##works but should transpose
            allprobes.columns=['Gene','Number of probes','Sequence'] #rememebr this! to assign names to df columns!
            #nonoverlapAll.columns=['Gene','RefSeq','Start','End','Bitscore','Sequence','Tm','GC','DeltaG']
    #team.columns =['Name', 'Code', 'Age', 'Weight'] 
        except:
            pass
    return allprobes #fix this

## need to write function to get specified number of probes if available! #numbersdict
   # selector = myselection[3][0::kmerdistance]
    #kmerspergene = myselection[0].isin 
    #return spacedselection

numbersfile = pd.read_csv('mygenenumbers.txt', sep='\t', header=None) #[1]## provided as csv file (generated based on known scRNAseq or other expression values)

aallprobes = getkmerspergene(myselection, numbersfile, kmerdistance=1) ## fix this


# write iterator to maximize probe numbers ??

