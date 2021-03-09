# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 16:07:50 2021

@author: thomas.hauling
"""

from Bio.Blast.Applications import NcbiblastnCommandline ## this works! outputs xml file


### BLAST module (use this for kmers <30 bp)
def blastn(kmerfasta, blastnout):
    query=kmerfasta
    out=blastnout
    blastdbpath="C:/Worky/DB/RefSeq/mouse/mouse.transcriptome"
    task="blastn-short"
    cpucores=8
    evalue=0.1
    outfmt='6 qseqid sseqid nident sstart send evalue bitscore'
    cline = NcbiblastnCommandline(task=task, query=query, db=blastdbpath, out=out, outfmt=outfmt, num_threads=cpucores, evalue=evalue)
    cline
    stout, stderr = cline()

### this is still taking TOO long - >18h - days (despite multicore option) - is only an option for shorter gene lists
## combine with kmer seqeunce to create same output format as with mmseq2