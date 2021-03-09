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
from collections import Counter
#import sqlalchemy

# cwd = os.getcwd()
# os.chdir(r'C:\Users\Thomas\Documents\Python_and_R')



# # ## query database, retrieve all entries of symbol=x, bitscore=y
# dbfile = 'dbout'
# symbolfile = 'mysymbols.txt' 
 

###########################################################################

# returns df of db rows for specified genes
def getkmersfromsymbols(dbfile, symbolfile): ##symbols, scores have to be tuples
    ''' retrieves kemrs from SQLite3 DB based on gene/ symbol list supplied by user '''
    
    symbol = 'symbol'  ## quirky hack to get cur.execute query line to work
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

# returns df of db rows for specified genes
def getkmers(dbfile, symbolfile): ##symbols, scores have to be tuples
    ''' retrieves kemrs from SQLite3 DB based on gene/ symbol list supplied by user '''
    
    symbol = 'symbol'  ## quirky hack to get cur.execute query line to work
    allrows = pd.DataFrame()
    try:
        con = sqlite3.connect(dbfile)
    except:
        print('could not connect to db')
    cur = con.cursor()
    # f = open(symbolfile, 'rt')
    # genes = tuple(f.read().split())
    genes = tuple(pd.read_csv(symbolfile, sep='\t', header=0)['Symbol'])
    for g in genes:
        cur.execute("SELECT * FROM mmseqt WHERE %s=?" % (symbol,), (g,)) #this works!
        #cur.execute("SELECT * FROM mmseqt WHERE %s=? AND WHERE %s=73" % (symbol, score), (g,)) #
        rows = pd.DataFrame(cur.fetchall()) ## reads db as list
        #print(rows)
        allrows = allrows.append(rows) # why is this not growing??
    #df = pd.read_sql_query("SELECT * FROM mmseqt", con) #, symbols[0]) # reads db as df
    con.close()
    allrows.columns=['Symbol','RefSeq','Start','End','Bitscore','Sequence','Tm','GC','DeltaG']
    return allrows 
    #seems to work! #QC: remove duplicate kmers if present

# myrows = getkmers(dbfile,symbolfile) # test func
#def filterkmers(myrows, minscore=73, mintm=55, maxtm=65, mingc=45, maxgc=55, mindg=-0.5, maxdg=0.5):
def filterkmers(myrows, minscore=1, mintm=1, maxtm=100, mingc=1, maxgc=100, mindg=-100, maxdg=100): ##make very permissive
    ''' applies threshold filters for multiple DNA/RNA parameters on db selection
        bitscore, Tm, GC-content, and calculates deltaG Gibbs free energy '''
        
    selection = myrows[(myrows['Bitscore'] >= minscore) & (myrows['Tm'].between(mintm,maxtm)) & (myrows['GC'].between(mingc,maxgc))]#.any())]
    ## calulate dG on selection
    selection.drop(['DeltaG'], axis=1)
    selection['DeltaG'] = selection['Sequence'].apply(dg) ##assumes 37degC I guess
    selection = selection[selection['DeltaG'].between(mindg,maxdg)]
    selection.columns=['Symbol','RefSeq','Start','End','Bitscore','Sequence','Tm','GC','DeltaG']
    return selection

# myfilteredkmers = filterkmers(myrows, mindg=-10, maxdg=10) ## this works but setting with copy warning
# myfilteredkmers_copy = myfilteredkmers
###### seems to work until here !!
# print(id(myfilteredkmers_copy) == id(myfilteredkmers)) ##spoiler: they are THE SAME object (same ID)
# print(id(myfilteredkmers))
#### not so relevant right now ########################
#myselection.reset_index(inplace=True)
#myselectiondict = myselection[0] = defaultdict(int)
# symcounts = Counter(myfilteredkmers[0]) # to get number of kmers per gene/ symbol dict!
######## not so relevant end ############################

def getspacedkmers(filteredkmers, kmerdistance=1): ## distance defaults to 1 to select all
    ''' applies distance filter on filtered kmers from previous step '''
    
    #filteredkmers = myfilteredkmers ##debugging
    filteredkmers.sort_values(['RefSeq','End'], inplace=True)
    allspaced = pd.DataFrame()
   # allprobes= pd.DataFrame()
   # numbersdict = pd.Series(numbersfile[1].values,index=numbersfile[0]).to_dict()
    uniquerefseqs=list(np.unique(filteredkmers['RefSeq'])) ##make non-redundant list of all RefSeq Accessions
    #rs=refseq[0]
    ## apply simple distance filter, has to be done PER RefSeq! check previous script
    #kmerdistance=20 ##fro debugging
    #newstr = 'start'
    for rs in uniquerefseqs:
        #newstr += rs
        myrefseq=filteredkmers[filteredkmers['RefSeq'] == rs] 
        #y = np.array(myrefseq)
        ## below still a constrcutino site:
        spaced = myrefseq[np.r_[True,np.diff(myrefseq['End'])>=kmerdistance]] 
        allspaced = allspaced.append(spaced)
        ## np.r_ is the way to go, but might discard all if dist too large
    allspaced.columns=['Symbol','RefSeq','Start','End','Bitscore','Sequence','Tm','GC','DeltaG']
    return allspaced

#spaced_all = pd.DataFrame()
## (generated based on known scRNAseq or other expression values)
# myspacedkmers = getspacedkmers(myfilteredkmers) ## now it seems to work??!!

# numbersfile = pd.read_csv('mygenenumbers.txt', sep='\t', header=None) #
# spacedkmers = myspacedkmers
## seems to work BUT with caveat - might discard too many...


################################################################################
## now get n kmers per gene, starting from first kmer
   # for no in set(nonoverlapAll[0]): # use this to iterate over genes
def getnumberedkmers(spacedkmers, numbersfile):
    ''' gets user-speciefied number of kmers from previous (filtered, spaced) kmer selection '''
    numbersfile = pd.read_csv(numbersfile, sep='\t', header=0)
    wanteddict = pd.Series(numbersfile['Number'].values,index=numbersfile['Symbol']).to_dict() ## user supplied numbers to dict
    havedict = Counter(spacedkmers['Symbol']) ## creates a dictionary obj[(symbol(key): kmercount(value)]
    #print(havedict.get('Robo2')) ## and it works
    enough = pd.DataFrame()
    notenough = pd.DataFrame()
    try:
        for wa in wanteddict.keys(): # or iterate of items
        #print(wa)
        #print(wanteddict.get(wa))
            if wanteddict.get(wa) <= havedict.get(wa):
                print('sufficient kmers available for ' + wa)
                wanted1 = spacedkmers[spacedkmers['Symbol'] == wa][:(wanteddict.get(wa))] #selects n probes according to wanteddict, starting from pos 0
                enough = enough.append(wanted1)
            else:
                print('only ' + havedict.get(wa) + 'kmers available for ' + wa)
                wanted2 = spacedkmers[spacedkmers['Symbol'] == wa][:(havedict.get(wa))] #selects n probes according to wanteddict, starting from pos 0
                notenough = notenough.append(wanted2)
 
    except:
        pass
    return enough, notenough

#numbersfile = pd.read_csv(r'C:\Users\Thomas\Documents\Python_and_R\mygenenumbers.txt', sep='\t', header=0) #[1]## provided as csv file 
# numbersfile = pd.read_csv('20neurogenes_numbers2.txt', sep='\t', header=0) #[1]## provided as csv file #make sure files contain expected number of columns..
#mynumberedkmers, mynumberedkmers_notenough = getnumberedkmers(myspacedkmers, numbersfile)

## write out selected kmers

def saveselectedkmers(enoughkmers, notenoughkmers, enoughname='enoughkmers', notenoughname='notenoughkmers'):
    ''' writes selected kmers from getnumberedkmers to file'''
    enoughkmers.to_csv(f'{enoughname}.csv', index=False)
    notenoughkmers.to_csv(f'{notenoughname}.csv', index=False)
    
    return enoughname, notenoughname #print('selected kmers saved')

# saveselectedkmers(mynumberedkmers, mynumberedkmers_notenough)


