# Summary

These scripts were the ones used in: <Insert Paper Citation Here>

This pipeline takes assembled RNA-sequencing data (.gtf format) in human or mouse, and then will be able to assemble the data into transcripts, predict the transcripts starting from transposable elements, and calculate expression of the transcripts in comparison to the whole. 

# Requirements

## Software

stringtie >= 1.3.3

samtools >= 1.3.1

python 2.7 (cPickle, pytabix 0.1)

R >= 3.4.1 (ggplot2, bsgenome.hsapiens.ucsc.hg38 (or genome of your choosing)) 

## Reference Files

These are the required reference files for the code to work. The `arguments.txt` file that is within the directory of the scripts or that can be specified by the user will define where these files are located. We have default files ready for hg38, and plan on providing other assemblies as well. 

### Default hg38

We have created a default set of reference files for hg38. For simple use, download the directory, place it in the rnapipeline folder, and extract. Make sure to update the arguments.txt file with the full path of these files (Explained below in Usage).

Download Link hg38: [External Download Link](https://wangftp.wustl.edu/~nshah/rnapipeline_public_link/rnapipelinerefhg38.tar.gz 'Compressed Directory')

### Gencode Dictionary

1. Download Gencode GTF reference file desired here: https://www.gencodegenes.org/

2. Sort the File

```
cat <GENCODE GTF>  | awk '{if($3=="transcript"||$3=="exon"||$3=="start_codon"){print}}' |  awk -F "; " '{print $0"\t"$2}' > <OUTPUT_sorted.gtf>`
```
 
3. Use custom script to create dictionary

```
genecode_to_dic.py <OUTPUT_sorted.gtf>
```
 
This will generate 2 files: (1) `genecode_plus.dic` and (2) `genecode_minus.dic`
 
4. Rename as needed

### Repeatmasker Files

Though this pipeline is optimized for looking at repetitive-element derived transcripts, any bed file for alternative promoter locations can be used. 

1. Download Repeatmasker Annotation. We use the one from the UCSC Table Browser: [External Link](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=693256623_kh0RR0o6vajdA2WLLTA8OeaAPNB6&clade=mammal&org=Human&db=hg38&hgta_group=rep&hgta_track=hg38Patch11&hgta_table=0&hgta_regionType=genome&position=chr12%3A20816734-20825794&hgta_outputType=bed&hgta_outFileName=repeats.bed 'UCSC Table Browser')

Note:
> If you would like to check for overlap with other features such as CpG islands or you need to create a custom file of your features, make sure it has the following format and tab-delimitted:
> 
> chromosome start end label wildcard strand

2. bgzip the file (bgzip comes with samtools)

```
bgzip rmsk.bed > rmsk.bed.gz
```

3. Create tabix index

```
tabix -p bed rmsk.bed.gz
```

4. Both the bgzipped file and the tabix index must be in the same directory

### Intron Annotations (Optional)

A useful feature that can be used is that the pipeline will annotate the first intron of each transcript based on a reference file. This can help in deciding whether a candidate has been previously annotated as an alternative transcript. 

Note:
> If you have already created the gencode dictionary by yourself in the 'Gencode Dictionary' step , then you have already done step 1 and 2 and you can skip to step 3

1. Download Gencode GTF reference file desired here: https://www.gencodegenes.org/

2. Sort the File

```
cat <GENCODE GTF>  | awk '{if($3=="transcript"||$3=="exon"||$3=="start_codon"){print}}' |  awk -F "; " '{print $0"\t"$2}' > <OUTPUT_sorted.gtf>`
```
 
3. Use custom script to create intron annotations

```
genecode_introns.py <OUTPUT_sorted.gtf>
```
 
This will generate 2 files: (1) `<OUTPUT_sorted.gtf>_introns_plus` and (2) `<OUTPUT_sorted.gtf>_introns_minus`

4. Sort the intron annotations bed file

```
sort -k1,1 -k2,2n -k3,3n <OUTPUT_sorted.gtf>_introns_plus > <OUTPUT_sorted.gtf>_introns_plus_sorted
sort -k1,1 -k2,2n -k3,3n <OUTPUT_sorted.gtf>_introns_minus > <OUTPUT_sorted.gtf>_introns_minus_sorted
```

5. bgzip the file (bgzip comes with samtools)

```
bgzip <OUTPUT_sorted.gtf>_introns_plus_sorted > <OUTPUT_sorted.gtf>_introns_plus_sorted.gz
bgzip <OUTPUT_sorted.gtf>_introns_minus_sorted > <OUTPUT_sorted.gtf>_introns_minus_sorted.gz
```

6. Create tabix index

```
tabix -p bed <OUTPUT_sorted.gtf>_introns_plus_sorted.gz
tabix -p bed <OUTPUT_sorted.gtf>_introns_minus_sorted.gz
```

7. Both the bgzipped file and the tabix index need to be in the same folder

### Gene Filter List (Optional)

For large sets of analysis, it might be computationally advantageous to limit the analysis to only a small set of genes for subsequent analysis. Create a file with a genesymbol per line. Example:

```
TP63
GAPDH
SYT1
PKIB
...
```

## Input Files

Our pipeline is optimized to work with `stringtie` generated gtf files. We have not tested it with other assembly software such as `cufflinks`. We plan on adding support for other input file styles in future releases. 

With raw FASTQ files, the following steps need to be taken:

1. QC: Using software such as [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ 'fastqc')

2. Adapter Trimming

3. Alignment: We have used either [STAR](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf 'star') or [HISAT2](https://ccb.jhu.edu/software/hisat2/manual.shtml 'hisat2'). Make sure that the -XS tag will be outputted since the splice-aware assembly software such as `stringtie` has the information it needs. 

Note:

> Future steps require sorted and indexed bam files. Thus, outputting sorted BAM files will be beneficial. 

4. Filtering: Due to problems with mapping to repetitive elements, we usually will filter to only include uniquelly mapped reads. 

5. Assembly: [STRINGTIE](https://ccb.jhu.edu/software/stringtie/ 'stringtie') Default parameters will work, but for discovering low coverage transcripts we sometimes set the following flags (-c 1 -m 100). Note, future steps will require being able to map the gtf file to the bam file, so the easiest way to do this would be to just change the extension of the files from `*.bam` to `*.gtf` (test1.bam to test1.gtf).

Example of Command for steps 4 & 5:

```
samtools view -q 255 -h EGAR00001163970_636782_1_1.rnaseq.fastqAligned.sortedByCoord.out.bam | stringtie - -o EGAR00001163970_636782_1_1.rnaseq.fastqAligned.sortedByCoord.out.gtf -m 100 -c 1
```

# Usage

Once all the gtf files are generated by stringtie for your experiments, the following steps can be used to annotate for alternative promoters from TEs. We have designed this pipeline to be flexible, and thus we provide explanations for each intermediate file that can be used for more custom usage of our code. 

## 1. Setup arguments.txt

The locations of the reference files need to be specified in a tab delimitted file. There should be no extra lines or headers. Here is an example:

```
rmsk	/bar/nshah/reference/rmskhg38.bed6.gz
gencodeplusdic	/bar/nshah/reference/genecode_plus_hg38.dic
gencodeminusdic	/bar/nshah/reference/genecode_minus_hg38.dic
focusgenes	/bar/nshah/reference/oncogenes_augmented.txt
plusintron	/bar/nshah/reference/gencode.v25.annotation.sorted.gtf_introns_plus_sorted.gz
minusintron	/bar/nshah/reference/gencode.v25.annotation.sorted.gtf_introns_minus_sorted.gz
```

### Required Arguments

rmsk: tabix formatted bed6 files with repeatmasker or other file that user wants to check start location for
gencodeplusdic: Dictionary of all gencode elements including introns and exons for the plus (+) strand
gencodeminusdic: Dictionary of all gencode elements including introns and exons for the minus (-) strand

Note:
> These are all the arguments that are needed. The following arguments.txt file would work:

```
rmsk	/bar/nshah/reference/rmskhg38.bed6.gz
gencodeplusdic	/bar/nshah/reference/genecode_plus_hg38.dic
gencodeminusdic	/bar/nshah/reference/genecode_minus_hg38.dic
```

### Optional Arguments

focusgenes: The program has two outputs (1) on a focus set of genes (2) with all genes. This file lists the genes that the user wants to filter for originally (Gene Filter List)
plusintron: Tabix file of all the plus strand introns (Intron Annotations)
minusintron: Tabix file of all the minus strand introns (Intron Annotations)

Note:
> If you do not want to use these options, remove them from the file. There should be no extra lines in the file or it will not work. 

## 2. Run annotation on each GTF file

Run the Program

```
rmskhg38_annotate_gtf_update_test_tpm.py <gtffile> <argumentfile.txt*>
```
*This is optional. If it is not included then the program will default to the argumentfile.txt within the rnapipeline folder.

Description of Output File(s):

(1) \<gtffile\>_annotated_test_all
 
(2) \<gtffile\>_annotated_filtered_test_all
 
(3) \<gtffile\>_annotated_test*
 
(4) \<gtffile\>_annotated_filtered_test*
 
 * Files (3) and (4) will only be produced if the user would like to focus the analysis on a set of genes such as a list of oncogenes and tumor suppresor genes

Excel Document with Description: [External Link](https://wangftp.wustl.edu/~nshah/rnapipeline_public_link/Transcript%20Annotation%20Description.xlsx 'GTF Annotation Output')

## 3. Aggregate annotations across all files in study

## 4. Calculate Read Information

## 5. Filter Candidates and Merge with Reference GTF

## 6. Calculate Transcript-Level Expression with Stringtie
