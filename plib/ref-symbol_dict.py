# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 14:30:18 2021

@author: thomas.hauling
"""


import pandas as pd
import mygene 

mg = mygene.MyGeneInfo()

markerlist = pd.read_csv(r'C:\Users\thomas.hauling\WORK\padlock_design\whole transcriptome\MRK_List1.rpt', sep='\t')
genelist = markerlist[markerlist['Marker Type']=='Gene']
genelist2 = list(set(genelist['Marker Symbol']))#[0:1000]

## get RefSeq IDs based on gene symbol
#symbols = mg.querymany(genelist['Marker Symbol'], scopes='symbol', species='mouse', as_dataframe=True)
symbols = mg.querymany(genelist2, scopes='symbol', species='mouse', fields='refseq.rna', as_dataframe=True)
symbols = symbols['refseq.rna'].dropna()
symbols = symbols.to_frame()


df = symbols[['refseq.rna']]
df.reset_index(inplace=True)
## Dimitris code - does the job 

d = {}
for index, row in df.iterrows():
    rna = row['refseq.rna']
    rna = [rna] if isinstance(rna, str) else rna
    for el in rna:
        if el in d.keys():
            d[el].append(row['query'])
        else:
            d[el] = [row['query']]

duplicates = [x for x in d.items() if len(x[1]) > 1]



dkeys = list(d.keys())
dvalues = list(d.values())
dvalcopy = dvalues

## in case there are double values per key in dict:
for sl in dvalues:
    print(len(sl))
    if len(sl)>1:
        try:
            del sl[0]
        except:
            pass
        
    #del sl[0]# if len(sl) >1
print(dvalues == dvalcopy)
dvalues_f = [i for sublist in dvalues for i in sublist]
ddf = pd.DataFrame(list(zip(dkeys, dvalues_f)), columns=['refseq', 'symbol'])
ddf.to_csv('d.csv', sep='\t', index=False)

ddf2=pd.read_csv('d.csv', sep='\t', header=0)
print(ddf == ddf2)
d2 = pd.Series(ddf2.symbol.values,index=ddf2.refseq).to_dict()
## seems to works!
#ddf2.set_index('refseq', inplace=True)
#d2 = ddf2.to_dict(orient='index')
