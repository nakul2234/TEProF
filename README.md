# Summary

These scripts were the ones used in: <Insert Paper Citation Here>

This pipeline takes assembled RNA-sequencing data (.gtf format) in human or mouse, and then will be able to assemble the data into transcripts, predict the transcripts starting from transposable elements, and calculate expression of the transcripts in comparison to the whole. 

# Requirements

## (1) Software

stringtie >= 1.3.3

samtools >= 1.3.1

python 2.7 (cPickle, pytabix 0.1)

R >= 3.4.1 (ggplot2, bsgenome.hsapiens.ucsc.hg38 (or genome of your choosing), Xmisc) 

## (2) Reference Files

These are the required reference files for the code to work. The `arguments.txt` file that is within the directory of the scripts or that can be specified by the user will define where these files are located. We have default files ready for hg38, and plan on providing other assemblies as well. 

### (A) Default hg38

We have created a default set of reference files for hg38. For simple use, download the directory, place it in the rnapipeline folder, and extract. Make sure to update the arguments.txt file with the full path of these files (Explained below in Usage).

Download Link hg38: [External Download Link](https://wangftp.wustl.edu/~nshah/rnapipeline_public_link/rnapipelinerefhg38.tar.gz 'Compressed Directory')

Note:
>If you use these default files then you can skip B-E. However, be sure to look at E if you want to use a list of genes different then the default list that is based on the publication. 

### (B) Gencode Dictionary

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

### (C) Repeatmasker Files

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

### (D) Intron Annotations (Optional)

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

### (E) Gene Filter List (Optional)

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
rmsk /bar/nshah/reference/rmskhg38.bed6.gz
rmskannotationfile /bar/nshah/reference/repeatmasker_description_uniq.lst
gencodeplusdic /bar/nshah/reference/genecode_plus_hg38.dic
gencodeminusdic /bar/nshah/reference/genecode_minus_hg38.dic
focusgenes /bar/nshah/reference/oncogenes_augmented.txt
plusintron /bar/nshah/reference/gencode.v25.annotation.sorted.gtf_introns_plus_sorted.gz
minusintron /bar/nshah/reference/gencode.v25.annotation.sorted.gtf_introns_minus_sorted.gz
```

### Required Arguments

**rmsk:** tabix formatted bed6 files with repeatmasker or other file that user wants to check start location for

**gencodeplusdic:** Dictionary of all gencode elements including introns and exons for the plus (+) strand

**gencodeminusdic:** Dictionary of all gencode elements including introns and exons for the minus (-) strand


Note:
> These are all the arguments that are needed. The following arguments.txt file would work:

```
rmsk	/bar/nshah/reference/rmskhg38.bed6.gz
rmskannotationfile /bar/nshah/reference/repeatmasker_description_uniq.lst
gencodeplusdic	/bar/nshah/reference/genecode_plus_hg38.dic
gencodeminusdic	/bar/nshah/reference/genecode_minus_hg38.dic
```

### Optional Arguments

**focusgenes:** The program has two outputs (1) on a focus set of genes (2) with all genes. This file lists the genes that the user wants to filter for originally (Gene Filter List)

**plusintron:** Tabix file of all the plus strand introns (Intron Annotations)

**minusintron:** Tabix file of all the minus strand introns (Intron Annotations)

Note:
> If you do not want to use these options, remove them from the file. There should be no extra lines in the file or it will not work. 

## 2. Run annotation on each GTF file

Run the Program

```
rmskhg38_annotate_gtf_update_test_tpm.py <gtffile> <argumentfile.txt>*
```
>*This is optional. If it is not included then the program will default to the argumentfile.txt within the rnapipeline folder.

Output File(s):

(1) \<gtffile\>_annotated_test_all
 
(2) \<gtffile\>_annotated_filtered_test_all
 
(3) \<gtffile\>_annotated_test*
 
(4) \<gtffile\>_annotated_filtered_test*
 
 >*Files (3) and (4) will only be produced if the user would like to focus the analysis on a set of genes such as a list of oncogenes and tumor suppresor genes

Excel Document with Description of Columns: [External Link](https://wangftp.wustl.edu/~nshah/rnapipeline_public_link/Transcript%20Annotation%20Description.xlsx 'GTF Annotation Output')

## 3. Process annotations files to get rough estimate of relative expression of transcript versus all other transcripts of gene

Run the program on the dataset desired. If you would like to stick to the selection of genes then run this on all the \<gtffile\>_annotated_filtered_test files. For all genes, run it on the \<gtffile\>_annotated_filtered_test_all files. 

```
annotationtpmprocess.py \<gtffile\>_annotated_filtered_test_all
```

Output File(s):

(1) \<gtffile\>_annotated_filtered_test_(all)_c

>Three additional columns are made: (1) Maximum coverage of gene of interest (2) Maximum tpm of gene of interest (3) Total of gene of interest.
>Note: For those transcripts that are from a TE but do not splice into a gene, these stats will be compared across all the transcripts like this. Thus the fraction will be very low.  

## 4. Aggregate annotation samples across samples

To aggregate all the data across samples, we have an R script that will be able to aggregate the statistics and summarize an output.

Arguments are positional, and thus making sure they are in correct order is important.

This script can be run standalone with specifying the parameters, or you can do this step in RStudio to add some more custom filters if needed. The columns in the final file must be the same, however for the subsequent steps to work. 

Options with Defaults:

-e \<ext_treatment\> (default: filtered): The label in the treatment file names that will identify them. If there is not a treatment versus untreated experimental design, the default will call everything as treatment since all the c files should have filter in their name. 

-l \<exon 1 length max\> (default: 2588): The maximum length of exon 1. We are using the 99th percentile of gencode v25 transcripts

-s \<exon skipping max\> (default: 2): Based on genomic contamination, assembly can create transcripts that have intron retention. These could be real, but oftentimes these look like noise. A maximum level of 2 is recommended.

-n \<sample total min\> (default: 1): Minimum number of samples that must have a candidate for it to be considered. With low number of samples 1 is recommended. For larger studies this can be increased.

-t \<treatment total min\> (default: 0): In case the user would like to only consider candidates that are present in a certain number of treatment samples.

-x \<treatment exclusive\> (default: no): In case the user would like a maximum for the untreated samples. This could be if the user wants to only consider treatment-exclusive samples. In that case, the user should put yes instead of no. 

-k \<keep none\> (default: no): There can be transcripts from TEs that have splicing, but do not end up splicing into a gene. By default these are removed. This can be changed to yes if the user would like to look at these.

-f \<filter for TEs\> (default: yes): Repeatmasker files include many repeats besides TE. This could be an interesting benchmark to compare to, but by default they are removed. The user can specify no if they would like to keep these. 

-a \<argument file\> (default: \<directory of pipeline\>/arguments.txt): The arguments.txt file location. By default, the arguments.txt file that is in the rnapipeline directory will be used. If another is desired for use then the fullpath of the file can be specified. 

```
aggregateProcessedAnnotation.R <options>
```

Output File(s):

(1) filter_combined_candidates.tsv: A file with every TE-gene transcript. This file is used for calculating read information in subsequent steps.
 
(2) initial_candidate_list.tsv: A summary of filter_combined_candidates.tsv for each unique candidate. Also lists the treatment and untreated files that the candidate is present in.
 
(3) Step4.RData: Workspace file with data loaded from R session. Subsequent steps load this to save time.

Note:
>This step will give ALL assembled transcripts. Assembly is very noisy, and thus the output of this step will have a high false positive rate.
>Filtering based on tpm or fraction of tpm can help.
>It is reccomended that the user looks at initial_candidate_list.tsv and assures that the candidates are in the format that they desire. 


## 5. Calculate Read Information

It is recommended that candidates are confirmed using read information. We have found this to greatly increase the specificity of the pipeline since assembly can make many errors especially without a reference.

We have developed a fast parser using samtools and bedtools to be able to get read information. There are separate programs for either single and or paired end reads. The only argument needed is the path of the bam files.

Note:
>Bam files should be sorted by position and indexed for this to work.
>Bam file naems should be the same as the gtf files that were the oriignal input for the pipeline except the extention is different.
>It is recommended that paired-end reads are used. Single-end performance is not as robust.

(A) Make the directory where the stats will be stored
```
mkdir filterreadstats
```

(B) Create commands to calculate read information

```
commandsmax_speed.py filter_combined_candidates.tsv <bam file location*>
```
Note:
>*the bam file location should be a full path and should end in a forward slash (e.g. /scratch/nakul/epikate/)

Output File(s):

(1) filterreadcommands.txt: A list of all the commands needed to calculate read information. These commands use a combination of samtools, bedtools, and custom scripts. For every combination of candidate and bamfile present in teh dataset, there will be a command.

(C) Run the commands
It is reccomended that [parallel_GNU](https://www.gnu.org/software/parallel/) or a similar method be used to run the list of commands in parallel. 

```
parallel_GNU -j 4 < filterreadcommands.txt
```
Note:
>(-j) will tell parallel_GNU the number of processes
>The alias for parallel_GNU might just parallel or something else on your server configuration

(D) Combine all the read information files

```
find ./filterreadstats -name "*.stats" -type f -maxdepth 1 -print0 | xargs -0 -n128 -P16 grep e > resultgrep_filterreadstatsdone.txt
cat resultgrep_filterreadstatsdone.txt | sed 's/\:/\t/g' > filter_read_stats.txt
```

This will create the file **filter_read_stats.txt** that has all the read information that is used in subsequent analysis.

## 6. Filter Candidates based on read information



## 7. Merge with Reference GTF

## 8. Annotate Merged GTF and Calculate Transcript-Level Expression with Stringtie

## 9. Aggregate Stringtie Information

This can also be done with the package [Ballgown](https://github.com/alyssafrazee/ballgown 'Ballgown Github') for more advanced statistics on transcript level quantification. In the case of simply wanting to find existence of TE-transcripts we have a helper script to aggregate relevant data. 
