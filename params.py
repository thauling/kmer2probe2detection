# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 09:44:02 2021

@author: Thomas
"""

## config file. specify k-mer and probedesign variables

### input
# directory that contains main.py
maindir = r'C:\Users\Thomas\Documents\Python_and_R'   
# directory for input and output
inoutdir = r'C:\Users\Thomas\Documents\Python_and_R\Python_scripts\kmer2probe2detection\inoutdata' 
# file that contains list of gene symbols for probe construction plus number of desired probes per gene (TAB-DELIMITED), column labels: Symbol, Number                                                    
genenumberlist = r'C:\Users\Thomas\Documents\Python_and_R\Python_scripts\kmer2probe2detection\inoutdata\neurogenes_numbers.txt'    
# file that contains barcodes for probe construction                                                                                                                        
barcodelist = r'C:\Users\Thomas\Documents\Python_and_R\Python_scripts\kmer2probe2detection\inoutdata\barcodes_twentymer_55.csv'    
# file that contains dye-probe sequences                                  
dyelist = r'C:\Users\Thomas\Documents\Python_and_R\inputdata\dyes.csv'  


### output 
# name of output alignment file                 
alignout = ''   
# length of k-mers in bp                                                      
kmerlength = 40   
# minimum Bitscore, min/max Tm, min/max GC percentage, min/max Gibbs Free Energy deltaG
minscore, mintm, maxtm, mingc, maxgc, mindg, maxdg = 60,30,90,40,60,-1,1           
# spacing between probes on bp                                                 
spacing = 20    
# name of k-mer selection from db
kmerselection = r'C:\Users\Thomas\Documents\Python_and_R\outputdata\neurokmers.csv'
# number of dyes for bridge assay (has to be an odd prime power)                                                                
dyes = 7         
# number of imaging rounds                                                              
rounds = 8                                                                    


### 

