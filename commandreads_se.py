#!/usr/bin/python
# programmer : nshah 
# usage: To be used in order to get read information from max files

import sys

# The distance that reads can be to overlap
distbuffer = 1000

bamfoldername = sys.argv[2]

filesuffix = sys.argv[3]
# This is the rest of the command that will create a searchable tabix file for downstream analysis


fout=open('filterreadcommands_se.txt','w')

with open(sys.argv[1],"r") as filteredfile:
	for eachline in filteredfile:
		temp = eachline.strip().split("\t")
		elementsvec = temp[24]
		location = temp[6]
		chromosome = temp[31]
		filename = temp[40]
		uniqid = temp[47]

		if location == 'intron':
			elementvecint = [int(x) for x in elementsvec.split(',')]
			startloc = min(elementvecint) - distbuffer
			endloc = max(elementvecint) + distbuffer
			commandroot = 'samtools view -b -h ' + bamfoldername + '/' + filename + filesuffix + ' ' + chromosome + ':' + str(startloc) + '-' + str(endloc)
		else:
			strand = temp[37]
			elementvecint = [int(x) for x in elementsvec.split(',')]
			if strand == '+':
				startloc = int(temp[32]) - distbuffer
				endloc = max(elementvecint) + distbuffer
			else:
				startloc = min(elementvecint) - distbuffer
				endloc = int(temp[33]) + distbuffer
			commandroot = 'samtools view -b -h ' + bamfoldername + '/' + filename + filesuffix + ' ' + chromosome + ':' + str(startloc) + '-' + str(endloc)


		commandfinish = '| bedtools bamtobed -i stdin | rmsk_annotate_bed_se.py | sort -k1,1 -k2,2n | bgzip > ./filterreadsse/' + uniqid + '--' + filename + '.tabse.bgz; tabix -p bed -f ./filterreadsse/' + uniqid + '--' + filename + '.tabse.bgz'
		commandall = commandroot + commandfinish
		print >> fout, commandall

			




