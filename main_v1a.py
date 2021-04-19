#! /usr/bin/env/ python3  ##shebang line to allow execution from commandline on unix based systems, will be ignored on win

#import flask
import os
import glob
import sys
import pandas as pd

sys.path.append('.') #modify sys path at runtime, adds current path that contains plib
sys.path.append('.\plib') #modify sys path at runtime, add path that contains plib

import params

sys.path.append(params.inoutdir) #modify sys path at runtime, add path that contains plib

from plib import kmersfromfasta
from plib import kmeralign    ## outputs mmseqs commands, kmers.fa plus commands need to be parsed to a UNIX-based OS (eg. WSL2 UBUNTU)
from plib import align2db
from plib import selectfromdb 
from plib import probeconstructor
from plib import detectionconstructor


### user input ####

validinput1 = ['makekmerdb', 'makeprobes'] 
validinput2 = ['cdnapadlock', 'rnapadlock', 'rnasmfish'] 
if (len(sys.argv) < 2 or sys.argv[1] not in validinput1): # or len(sys.argv) > 3):
    print('specify makekmerdb or makeprobes plus probetype')
if sys.argv[1] == validinput1[1] and len(sys.argv) < 3: # or sys.argv[2] not in validinput2:
    print('specify cdnapadlock, rnapadlock or rnasmfish')   
# elif sys.argv[2] not in validinput2:
#         print('specify cdnapadlock, rnapadlock or rnasmfish')
if sys.argv[1] == validinput1[0] and len(sys.argv) > 2:
    print(f'{validinput1[0]} does not take additional arguments')     

else: 
    print(f'running {sys.argv[0]} in {sys.argv[1]} mode')
#### user input end ################################################################################################### 
    ## makekmerdb chosen ##
    if sys.argv[1] == validinput1[0]:
        os.chdir(params.inoutdir)
        print('this can take a while')
        #kmerlength = params.kmerlength
        kmerfasta = f'kmers{params.kmerlength}.fa'
##### generate kmers from fasta
        #os.chdir(params.outdir)
        kmersfromfasta.makekmers(params.kmerlength, kmerfasta)
        print('k-mer fastafile created')
        # ### 3) mmseqs command: specify output text file for mmseqs command
        instructout = f'kmers{params.kmerlength}mmseqs_instruct.txt'
        print('instruction file for mmseqs generated')
##### generate mmseqs command for alignment
        kmeralign.kmeralign(kmerfasta, instructout) 
##### put mmseqsout to database ## need to test if file exists! print 'create mmseqs-out file'if file not found
        mmseqout = f'{kmerfasta[:-3]}_outdb'
        symbols = f'{mmseqout}_symbols.txt' 
        dbout = f'kmers{params.kmerlength}_sqlite.db'
        if os.access(mmseqout, os.R_OK) == True:
            print(f'{mmseqout} found. Calculating Tm, GC content, deltaG and parsing results to database')
            align2db.ref2sym_chunk(mmseqout) ## note: expects d.csv and mmseqout in os.cwd
            print(f'{symbols} created')
            align2db.out2db(mmseqout, symbols, dbout)
            print(f'{mmseqout} transferred to databsse')
        else: print('mmseqs-out- file not found - run mmseqs first')
        
####### makeprobes chosen  ################################################################################################        
    if sys.argv[1] == validinput1[1]:
        # ### 5) select genes from DB
        dbout = f'kmers{params.kmerlength}_sqlite.db'
        os.chdir(params.inoutdir)
        kmerselection = selectfromdb.getkmers(dbout, params.genenumberlist) # symbols = symbols to select from DB, provided as csv
        kmerselection = selectfromdb.filterkmers(kmerselection)
        kmerselection = selectfromdb.getspacedkmers(kmerselection, params.spacing) 
        kmerselection, kmerselection_notenough = selectfromdb.getnumberedkmers(kmerselection, params.genenumberlist) # specify k-mer number per gene symbol (tab-separated file, header MUST contain Symbol, Number)
        ##save kmers as csv
        enoughname, notenoughname = selectfromdb.saveselectedkmers(kmerselection, kmerselection_notenough) # weird output, too tired
        print('selected kmers saved')
        kmerselection = pd.read_csv(f'{enoughname}.csv', header=0)
        barcodelist = pd.read_csv(params.barcodelist, header=0)

        if sys.argv[1] == validinput1[1] and sys.argv[2] == validinput2[0]:
             probes = probeconstructor.cdnapadlocks(kmerselection, barcodelist)
             #print(f'{sys.argv[1]} done')
        if sys.argv[1] == validinput1[1] and sys.argv[2] == validinput2[1]:
             probes = probeconstructor.rnapadlocks(kmerselection, barcodelist)
             #print(f'{sys.argv[1]} done')
        if sys.argv[1] == validinput1[1] and sys.argv[2] == validinput2[2]:
             probes = probeconstructor.rnasmfish(kmerselection, barcodelist)
        print(f'{sys.argv[1]} {sys.argv[2]} done')
#     ### 7) construct detection assay: specify csv files that constain sequences of probes and dye probes
        mykmers = pd.read_csv(f'{sys.argv[2]}.csv', header=0) #use length of this to get one barcode per kmer 
        mydyes = pd.read_csv(params.dyelist, header=0)
        uniquesymbols = mykmers['Symbol'].unique() #use length of this to get one barcode per gene symbol
# ### create detection oligos for selected probes (currently creates as many codes as there are barcodes, does not remove duplicates)
        mybridgeprobes, mycode = detectionconstructor.makedetection(params.dyes,params.rounds,len(uniquesymbols)) #only works for colors up to 7, needs at least as many rounds as colors. can return duplicates if too few rounds/ dyes are chosen. returns 1004 columns but should return same as length of probes.csv returns error if len=1
        print('bridge probes done')




## notes:
        ## add gene-symbol and barcode column to bridge probes output
        ## solve weird, 'gappy' output after numberedkmers has been applied 




