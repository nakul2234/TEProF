#! /usr/bin/env python
# programmer : nshah 
# usage: To be used to get summarized FPKM fractions for each transcript from a stringtie output file

import sys
import os
import glob
import numpy as np

filename = sys.argv[1]
filteredcontent = [i.strip().split('\t') for i in open(filename).readlines()[1:]]
referencecontent = np.loadtxt(filename, dtype=object, delimiter='\t', skiprows=1)
foutfractot=open('{}_frac_tot'.format(sys.argv[1]),'w')
foutfracmax=open('{}_frac_max'.format(sys.argv[1]),'w')
fouttot=open('{}_tot'.format(sys.argv[1]),'w')
for line in filteredcontent:
	geneline = line[9]
	transcriptidline = line[5]
	fpkmline = line[-1]
	subset = referencecontent[np.where(referencecontent[:,9] == geneline)]
	fpkmlist = subset[:,-1].astype(np.float)
	
	#If there is no expression of gene then no division will be performed
	if (max(fpkmlist) == 0.0):
		print >> foutfractot, '\t'.join([transcriptidline] + [str(0.0)])
		print >> foutfracmax, '\t'.join([transcriptidline] + [str(0.0)])
		print >> fouttot, '\t'.join([transcriptidline] + [str(0.0)])
	
	#If at least one isoform is expressed then there will be a collapsing
	else:
		print >> foutfractot, '\t'.join([transcriptidline] + [str(float(fpkmline)/sum(fpkmlist))])
		print >> foutfracmax, '\t'.join([transcriptidline] + [str(float(fpkmline)/max(fpkmlist))])
		print >> fouttot, '\t'.join([transcriptidline] + [str(sum(fpkmlist))])
	
	