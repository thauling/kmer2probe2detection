# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 11:58:23 2021

@author: Thomas
"""
import pandas as pd
# actual probe constructor, takes approved, selected kmers as input

## options:
# cDNA directed padlocks
# RNA directed padlocks
# RNA directed smFISH probes

#input:
# kmers
# barcode library, barcodes / gene

#output:
# barcoded probes
# gene: barcode(s) dict (csv format)
# print('specify probetype: 1 cDNA padlocks 2 RNA padlocks 3 smFISH')
# probetype = input('Type 1,2 or 3 \n')

# def getselectedkmers(probetype, selectedkmers=0, barcodes=0): ##user provides selectedkmers and barcodes as csv
#     '''build DNA/RNA molecular probes sequences from kmers'''
    
#     print(f'probe type {probetype} has been selected')
#     return probetype


# #selectedkmers = 
# mynewprobes = getselectedkmers(probetype) 
def revcomp(kmers=pd.DataFrame()):
    '''generates reverse complementary seqs from DNA/RNA, returns pd.DF'''
    rcdict = {'A':'T','T':'A','G':'C','C':'G'}
        #seq = 'TTGCCACT'
    seqrclist = []
    for seq in kmers['Sequence']:
        #rcseq = "".join(rcdict.get(base, base) for base in seq)[::-1]
        seqrclist.append("".join(rcdict.get(base, base) for base in seq)[::-1])
    return pd.DataFrame(seqrclist, columns=['Sequence rc'])


def cdnapadlocks(kmers, barcodes, barcodedict=pd.DataFrame()): #two options: list of barcodes OR barcode-gene LUT
    '''create cDNA directed padlocks '''
    assignedbc = []
    seqhalf = (kmers['Sequence'].str.len()/2)[1].astype(int) 
    threep = (kmers['Sequence'].str[:seqhalf])
    fivep = (kmers['Sequence'].str[seqhalf:])
    #barcode = barcodes[:kmers.shape[0]] # create sereis with as many barcodes as there are kmers
    barcodeset = set(barcodes['Barcode'])
    symbolset = set(kmers['Symbol'])
    if barcodedict.empty == 1:
        symbolbarcodedict = dict(zip(symbolset,barcodeset)) ## works but picks RANDOM barcodes
    else:
        symbolbarcodedict = dict(zip(barcodedict['Symbol'],barcodedict['Barcode']))
    for s in kmers['Symbol']:
        assignedbc.append(symbolbarcodedict.get(s))
        #print(assignedbc)
    kmers['Barcode'] = pd.DataFrame(assignedbc)     
    kmers['Probe'] = fivep + kmers['Barcode'] + threep
    kmers.to_csv('cdnapadlock.csv', index=False)    
    return kmers

def rnapadlocks(kmers, barcodes, barcodedict=pd.DataFrame()): ## same as cDNA padlocks but needs revcomp function
    '''create RNA directed padlocks '''
    assignedbc = []
    kmers['Sequence rc'] = revcomp(kmers)
    seqhalf = (kmers['Sequence rc'].str.len()/2)[1].astype(int) 
    threep = (kmers['Sequence rc'].str[:seqhalf])
    fivep = (kmers['Sequence rc'].str[seqhalf:])
    #barcode = barcodes[:kmers.shape[0]] # create sereis with as many barcodes as there are kmers
    barcodeset = set(barcodes['Barcode'])
    symbolset = set(kmers['Symbol'])
    if barcodedict.empty == 1:
        symbolbarcodedict = dict(zip(symbolset,barcodeset)) ## works but picks RANDOM barcodes
    else:
        symbolbarcodedict = dict(zip(barcodedict['Symbol'],barcodedict['Barcode']))
    for s in kmers['Symbol']:
        assignedbc.append(symbolbarcodedict.get(s))
        #print(assignedbc)
    kmers['Barcode'] = pd.DataFrame(assignedbc)     
    kmers['Probe'] = fivep + kmers['Barcode'] + threep
    kmers.to_csv('rnapadlock.csv', index=False)    
    return kmers

def rnasmfish(kmers, barcodes, barcodedict=pd.DataFrame()): ## same as RNA padlocks but different string concatenation in assembly step
    '''create RNA directed (sm)FISH probes'''
    assignedbc = []
    kmers['Sequence rc'] = revcomp(kmers)
    #seqhalf = (kmers['Sequence rc'].str.len()/2)[1].astype(int) 
    #threep = (kmers['Sequence rc'].str[:seqhalf])
    #fivep = (kmers['Sequence rc'].str[seqhalf:])
    #barcode = barcodes[:kmers.shape[0]] # create sereis with as many barcodes as there are kmers
    barcodeset = set(barcodes['Barcode'])
    symbolset = set(kmers['Symbol'])
    if barcodedict.empty == 1:
        symbolbarcodedict = dict(zip(symbolset,barcodeset)) ## works but picks RANDOM barcodes
    else:
        symbolbarcodedict = dict(zip(barcodedict['Symbol'],barcodedict['Barcode']))
    for s in kmers['Symbol']:
        assignedbc.append(symbolbarcodedict.get(s))
        #print(assignedbc)
    kmers['Barcode'] = pd.DataFrame(assignedbc)     
    kmers['Probe'] = kmers['Barcode'] + kmers['Sequence rc'] + kmers['Barcode']
    kmers.to_csv('rnasmfish.csv', index=False)    
    return kmers
    


