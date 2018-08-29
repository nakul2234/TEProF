#!/usr/bin/python
# programmer : nshah 
# usage: To be used in order to get read information from max files

import sys

distbuffer = 1000

folderdic = {}
with open("/scratch/twlab/nshah/TCGA_IDs_folder","r") as foldertcga:
	for eachline in foldertcga:
		temp = eachline.strip().split("\t")
		folderdic[temp[0]] = temp[1]

# This is the rest of the command that will create a searchable tabix file for downstream analysis


fout=open('filterreadcommands.txt','w')

with open(sys.argv[1],"r") as filteredfile:
	for eachline in filteredfile:
		temp = eachline.strip().split("\t")
		number1 = temp[7]
		number2 = temp[20]
		elementsvec = temp[24]
		location = temp[6]
		chromosome = temp[31]
		startlocTE = temp[32]
		endlocTE = temp[33]
		filename = temp[44]
		folder = temp[48]
		uniqid = temp[57]
		strand = temp[37]

		if location == 'intron':
			elementvecint = [int(x) for x in elementsvec.split(',')]
			startloc = min(elementvecint) - distbuffer
			endloc = max(elementvecint) + distbuffer
			commandroot = 'samtools view -q 255 -h /ref/tcga/' + folderdic[folder] + '/' + folder + "/*/" + filename + '_gdc_realn_rehead.bam ' + chromosome + ':' + str(startloc) + '-' + str(endloc) + " > " + uniqid + '--' + filename + ".sam" + ' ; samtools view -q 255 /ref/tcga/' + folderdic[folder] + '/' + folder + "/*/" + filename + '_gdc_realn_rehead.bam ' + chromosome + ':' + str(startlocTE) + '-' + str(endlocTE) + ' | cut -f1  | sort | uniq > ' + uniqid + '--' + filename + 'IDs.txt ; cat <(samtools view -H ' + uniqid + '--' + filename + '.sam) <(grep -w -F -f ' + uniqid + '--' + filename + 'IDs.txt ' + uniqid + '--' + filename + '.sam)'
		else:
			strand = temp[37]
			elementvecint = [int(x) for x in elementsvec.split(',')]
			if strand == '+':
				startloc = int(temp[32]) - distbuffer
				endloc = max(elementvecint) + distbuffer
			else:
				startloc = min(elementvecint) - distbuffer
				endloc = int(temp[33]) + distbuffer
			commandroot = 'samtools view -q 255 -h /ref/tcga/' + folderdic[folder] + '/' + folder + "/*/" + filename + '_gdc_realn_rehead.bam ' + chromosome + ':' + str(startloc) + '-' + str(endloc) + " > " + uniqid + '--' + filename + ".sam" + ' ; samtools view -q 255 /ref/tcga/' + folderdic[folder] + '/' + folder + "/*/" + filename + '_gdc_realn_rehead.bam ' + chromosome + ':' + str(startlocTE) + '-' + str(endlocTE) + ' | cut -f1  | sort | uniq > ' + uniqid + '--' + filename + 'IDs.txt ; cat <(samtools view -H ' + uniqid + '--' + filename + '.sam) <(grep -w -F -f ' + uniqid + '--' + filename + 'IDs.txt ' + uniqid + '--' + filename + '.sam)'


		commandfinish = ' | samtools sort -n -m 1G - | bedtools bamtobed -i stdin 2>/dev/null | ~/scripts/rmsk_annotate_bedse_speed.py ' + chromosome + ',' + startlocTE + ',' + endlocTE + ' ' + elementsvec + ' ' + strand + ' ' + number1 + ' ' + number2 + ' > ../filterreadstats/' +  uniqid + '--' + filename + ".stats ; " + "rm " +  uniqid + '--' + filename + ".sam ; rm " +  uniqid + '--' + filename + 'IDs.txt'
	
		commandall = commandroot + commandfinish
		print >> fout, commandall

			




