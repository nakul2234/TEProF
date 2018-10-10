# Summary

These scripts were the ones used in: <Insert Paper Citation Here>

This pipeline takes assembled RNA-sequencing data (.gtf format) in human or mouse, and then will be able to assemble the data into transcripts, predict the transcripts starting from transposable elements, and calculate expression of the transcripts in comparison to the whole. 

# Requirements

## Software

stringtie >=1,3,3
samtools >=1.3.1
python 2.7 (cPickle, pytabix 0.1)
R >=3.4.1 (ggplot2, bsgenome.hsapiens.ucsc.hg38 (or genome of your choosing)) 

## Reference Files

These are the required reference files for the code to work. The `arguments.txt` file that is within the directory of the scripts or that can be specified by the user will define where these files are located. We have default files ready for hg38, and plan on providing other assemblies as well. 

### Default hg38

We have created a default set of reference files for hg38. For simple use, download the directory, place it in the rnapipeline folder, and extract. The arguments.txt file by default has these filenames and paths set. 

Download Link hg38: 

### Gencode Dictionary

1. Download Gencode GTF reference file desired here: https://www.gencodegenes.org/
2. Sort the File

`cat <GENCODE GTF>  | awk '{if($3=="transcript"||$3=="exon"||$3=="start_codon"){print}}' |  awk -F "; " '{print $0"\t"$2}' > <OUTPUT_sorted.gtf>`
 
3. Use custom script to create dictionary

`genecode_to_dic.py <OUTPUT_sorted.gtf>`
 
This will generate 2 files: (1) `genecode_plus.dic` and (2) `genecode_minus.dic`
 
4. Rename as needed

### Repeatmasker Files

### Intron Annotations (Optional)

### Gene Filter List (Optional)

## Input Files

# Usage
