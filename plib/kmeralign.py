# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 17:44:49 2021

@author: thomas.hauling
"""

# kmeralign

#inputfile = 'mmseqs_instruct.txt'
#kmerfasta = 'kmers30'

def kmeralign(kmerfasta, instructout, targetdb='rnadb'):
    mmseqs_input = (f'''mmseqs createdb {kmerfasta} {kmerfasta[:-3]}_db \n mmseqs search {kmerfasta} {kmerfasta[:-3]}_db {targetdb} {kmerfasta[:-3]}_outdb tmp --search-type 3 \n mmseqs convertalis {kmerfasta[:-3]}_db {targetdb} {kmerfasta[:-3]}_outdb --format-output "query,target,evalue,bits,qseq" --search-type 3''')  
    print(mmseqs_input, file=open(instructout, 'w'))
#with open(inputname, "w") as text_file:
 #   text_file.write(mmseqs_input)


#eg:
# 1. create query and search DBs from fasta files (kmers and transcriptome/ genome fasta from NCBI)
# mmseqs createdb query.fa queryDB 
# 2. search query DB  against nucleotide DB
# mmseqs search queryDB targetDB query_outdb tmp --search-type 3
# 3. then convert output DB to blastn-like output txt file
# mmseqs convertalis queryDB targetDB result.tsv --format-output "query,target,evalue,bits,qseq" --search-type 3
# does the job: query ID, acc Nr, scores, actual kmer sequence of query

#### results in output textfile that contains 
##### kmer (query) ID, RefSeq Acc. Nr of matching gene, e-value, bit-score, kmer seq