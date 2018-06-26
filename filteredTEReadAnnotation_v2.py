#!/usr/bin/python
# programmer : nshah 
# usage: To be used in order to annotate a GTF file for presence of transcripts that begin in transposable elements. 

from __future__ import division #integer division has to be //, now all / are floating point
import fileinput
import cPickle as pickle
import sys
import subprocess as sp
import re
import math
import copy
import tabix
import os 

def annotateTE(chrTE, startTE, endTE, strandCandidate, geneCandidate, numCandidate, fileBase):
    
    filename = '/bar/nshah/links/hlee.stringtie/filterreadspe/{}.tab.bgz'.format(fileBase)

    #This checks the size. If the paired end file is empty, then this likely means that it is s single-end file

    statinfo = os.stat(filename)
    sizepe = statinfo.st_size
    filetype = 'pe'

    # If the paired end index file is empty, it means the file is single end and that needs to be used
    #if sizepe == 28:
    #    filename = '/scratch/twlab/nshah/filter12045/{}.tabse.bgz'.format(fileBase)
    #    filetype = 'se'

    rmsk = tabix.open(filename)
    
    res = []
    tmp = ''
        
    try:
        tmp=rmsk.query(chrTE,startTE,endTE)
    except:
        tmp = ''
        
    for i in tmp:
        res.append(i)
    
    if res:
        subfam = res
    else:
        subfam = []
            
    # There should only be 1 Subfamily that is returned. If none is returned, then the transcript is not fruther annotate
    if subfam == []:
        returncounts = [0, 0, 0, filetype]
    else:
        #There are cases where TE coordinates overlap. In that case, I do not want all the reads 
        subfamcandidate = filter(lambda x: x[1] == str(startTE), subfam)
        
        countreads = len(subfamcandidate)/2
        countreadsstart = len(filter(lambda x: x[7] == strandCandidate and x[8] == geneCandidate and x[9] == 'exon' and int(x[10]) >= int(numCandidate) and x[4] == 'start', subfamcandidate))
        countreadsend = len(filter(lambda x: x[7] == strandCandidate and x[8] == geneCandidate and x[9] == 'exon' and x[4] == 'end', subfamcandidate))
        returncounts = [countreads, countreadsstart,  countreadsend, filetype]
         
    return returncounts

with open(sys.argv[1],"r") as candidatelist:
    
    # The file name of the peak file has information on the chromosome, start, end, and strand of the peak
    for eachline in candidatelist:
        
        temp = eachline.strip().split("\t")
        filename = temp[40]
        resultcounts = annotateTE(temp[31], int(temp[32]), int(temp[33]), temp[37], temp[15], temp[20], temp[-1] + '--' + filename)
        resultcounts = [str(x) for x in resultcounts]
        outputnew = '\t'.join([temp[-1]] + [filename] + resultcounts)
        print outputnew
        
        
        
        


        
        
        
    
                
                
                