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

rmsk = tabix.open('/home/nakul.m.shah/reference/rmskhg38.bed6.gz')
plustab = tabix.open('/home/nakul.m.shah/reference/genecode_plus.tab.bgz')
minustab = tabix.open('/home/nakul.m.shah/reference/genecode_minus.tab.bgz')

peaktotalsize = 0 #Will be added up from the bed/peak file

TE={}
TEsub={}
TEdic={}

# This is a file form Daofeng that has the class, family, and subfamily. I will replace this with the rmsk subfam dictionary I made earlier.
# Format: Subfamily     Class   Family
with open("/home/nakul.m.shah/reference/TE.lst","r") as fin:
        for eachline in fin:
                temp = eachline.strip().split('\t')
                TE[temp[0]]=temp[1]
                TEsub[temp[0]]=temp[2]
                TEdic[temp[2]]=temp[1]

def annotateTE(start_t, chr_t, TE):
    res = []
    tmp = ''
        
    try:
        tmp=rmsk.query(chr_t,start_t-1,start_t+1)
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
        transcriptTE = ["None", "None", "None", "None", "None", "None"]
    else:    
        transcriptTE = subfam[0]
        subfamTE = transcriptTE[3]
        
        # I want to specifically look at TE for now, but maybe others later.
        if subfamTE not in TE.keys():
            transcriptTE = ["None", "None", "None", "None", "None", "None"]
    
    return transcriptTE

def annotateloc(start_t, chr_t, strand_t):
    
    if strand_t == "+":
        tabq = plustab
    else:
        tabq = minustab
                
    res = []
    tmp = ''
        
    try:
        tmp=tabq.query(chr_t,start_t,start_t+1)
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
        transcriptanno = ["None", "None", "None"]
    elif len(subfam) == 1:    
        transcriptanno = subfam[0][3:6]
    else:
        typelist = [sublist[4] for sublist in subfam]
        if 'exon' in typelist:
                transcriptanno = subfam[0][3:6]
        else:
                transcriptanno = subfam[0][3:6]        
        
    
    return transcriptanno



allowedchr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

      
# The file name of the peak file has information on the chromosome, start, end, and strand of the peak
for eachline in fileinput.input():
    
    finalout = "None"
    finalout1b = "None"
    finalout2 = "None"
    finalout2b = "None"
    
    temp = eachline.strip().split("\t")
    
    # Obtain pertinent transcript information
    qualityread = int(temp[4])
    
    #If the first or second mate are not aligned properly then this read is disregarded
    #Also, it must be a uniquelly mapped read (255 in STAR or 60 in hisat)
    if qualityread < 60:
            continue
    
    chr_read = temp[0]
    
    if chr_read not in allowedchr:
            continue

    start_read = int(temp[1])
    end_read = int(temp[2])
    startTE = annotateTE(start_read, chr_read, TE)
    endTE = annotateTE(end_read, chr_read, TE)
    
    #print([start_read, end_read, startTE, endTE, temp[6]])
    
    if (startTE[0] == "None" and endTE[0] == "None"):
        continue

    #If both are TE-derived, I need to determine how I should annotate them.
    if (startTE[0] != "None" and endTE[0] != "None"):
        
        # Plus strand annotation of read
        
        startanno = annotateloc(start_read, chr_read, '+')
        endanno = annotateloc(end_read, chr_read, '+')
        
        #If both ends of the reads are in exons or if neither one is in exons then I do not want to annotate a + strand position for the read
        if not ((startanno[1] == "exon" and endanno[1] == "exon")):
                if (endanno[1] == "exon"):
                        finalout = startTE[0:4] + ['start'] + [str(start_read)] + [str(end_read)] + ["+"] + endanno
                elif (startanno[1] == "exon"):
                        finalout = endTE[0:4] + ['end'] + [str(start_read)] + [str(end_read)] + ["+"] + startanno
                else:
                        finalout = startTE[0:4] + ['start'] + [str(start_read)] + [str(end_read)] + ["+"] + endanno
                        finalout1b = endTE[0:4] + ['end'] + [str(start_read)] + [str(end_read)] + ["+"] + startanno
            
            
        # Minus strand annotation of read
        
        startanno2 = annotateloc(start_read, chr_read, '-')
        endanno2 = annotateloc(end_read, chr_read, '-')
        
        #If both ends of the reads are in exons or if neither one is in exons then I do not want to annotate a + strand position for the read
        if not ((startanno2[1] == "exon" and endanno2[1] == "exon")):
                if (endanno2[1] == "exon"):
                        finalout2 = startTE[0:4] + ['end'] + [str(start_read)] + [str(end_read)] + ["-"] + endanno2
                elif (startanno2[1] == "exon"):
                        finalout2 = endTE[0:4] + ['start'] + [str(start_read)] + [str(end_read)] + ["-"] + startanno2
                else:
                        finalout2 = startTE[0:4] + ['end'] + [str(start_read)] + [str(end_read)] + ["-"] + endanno2
                        finalout2b = endTE[0:4] + ['start'] + [str(start_read)] + [str(end_read)] + ["-"] + startanno2
        
    #If only one is TE-derived then based on if it is the start (5' + strand) of the read or the end (3- + strand) of the read
    else:
        if (startTE[0] != "None"):
                finalout = startTE[0:4] + ['start'] + [str(start_read)] + [str(end_read)] + ["+"] + annotateloc(end_read, chr_read, '+')
                finalout2 = startTE[0:4] + ['end'] + [str(start_read)] + [str(end_read)] + ["-"] + annotateloc(end_read, chr_read, '-')
        else:
                finalout = endTE[0:4] + ['end'] + [str(start_read)] + [str(end_read)] + ["+"] + annotateloc(start_read, chr_read, '+')
                finalout2 = endTE[0:4] + ['start'] + [str(start_read)] + [str(end_read)] + ["-"] + annotateloc(start_read, chr_read, '-')
 
    if (finalout != "None"):
        print('\t'.join(finalout))
        
    if (finalout1b != "None"):
        print('\t'.join(finalout1b))

    if (finalout2 != "None"):
        print('\t'.join(finalout2))

    if (finalout2b != "None"):
        print('\t'.join(finalout2b))
        
        
        
    
                
                
                
        
        
    
                
                
                
