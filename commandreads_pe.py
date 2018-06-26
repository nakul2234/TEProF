#!/usr/bin/python
# programmer : nshah 
# usage: To be used in order to get read information from max files

import sys

# The distance that reads can be to overlap
distbuffer = 1000

bamfoldername = sys.argv[2]

filesuffix = sys.argv[3]
# This is the rest of the command that will create a searchable tabix file for downstream analysis

#filenameraw = sys.argv[1].split("/")[-1]
#filename = filenameraw.split(".")[0]

fout=open('filterreadcommands_pe_combined.txt','w')

with open(sys.argv[1],"r") as filteredfile:
	for eachline in filteredfile:
		temp = eachline.strip().split("\t")
		elementsvec = temp[24]
		location = temp[6]
		chromosome = temp[31]
		# R code that generates uniq id- paste(res2$subfamTE,res2$startTE,res2$gene1,res2$exonintron1, res2$number1, res2$gene2, res2$exonintron2, res2$number2,res2$transcriptstart2,sep = "_")
		filename = temp[40]
		locationsofuniqid = [34, 32, 2, 6, 7, 15, 19, 20, 22]
		uniqid = '_'.join([temp[i] for i in locationsofuniqid])

		if location == 'intron':
			elementvecint = [int(x) for x in elementsvec.split(',')]
			startloc = min(elementvecint) - distbuffer
			endloc = max(elementvecint) + distbuffer
			commandroot = 'samtools view -h -f 2 -F 4 ' + bamfoldername + '/' + filename + filesuffix + ' ' + chromosome + ':' + str(startloc) + '-' + str(endloc)
		else:
			strand = temp[37]
			elementvecint = [int(x) for x in elementsvec.split(',')]
			if strand == '+':
				startloc = int(temp[32]) - distbuffer
				endloc = max(elementvecint) + distbuffer
			else:
				startloc = min(elementvecint) - distbuffer
				endloc = int(temp[33]) + distbuffer
			commandroot = 'samtools view -h -f 2 -F 4 ' + bamfoldername + '/' + filename + filesuffix + ' ' + chromosome + ':' + str(startloc) + '-' + str(endloc)


		commandfinish = ' | samtools sort -n - | bedtools bamtobed -bedpe -i stdin | rmsk_annotate_bedpe.py | sort -k1,1 -k2,2n | bgzip > ./filterreadspe/' + uniqid + '--' + filename + '.tab.bgz; tabix -p bed -f ./filterreadspe/' + uniqid + '--' + filename + '.tab.bgz'
		commandall = commandroot + commandfinish
		print >> fout, commandall

			




