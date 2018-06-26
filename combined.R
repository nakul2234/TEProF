setwd('/tavern/hlee/PC/WTS/2_stringtie')

#HJoo Data script
files_annotated <- list.files(pattern = "filtered_test_all$")
firsttime = TRUE
i = 0
for (file_name in files_annotated){
  
  if (file.size(file_name) == 0){
    next
  }
  filtered_table = read.delim(file_name, sep = "\t", stringsAsFactors = FALSE, header = FALSE)
  
  columnlabels <- c("transcriptname", "type1", "gene1", "chr1", "start1", "end1", "exonintron1", "number1" ,"startcodon1" , "transcriptstart1", "transcriptend1", "elements1", "id1", "splicing", "type2", "gene2", "chr2", "start2", "end2", "exonintron2", "number2" ,"startcodon2", "transcriptstart2", "transcriptend2", "elements2", "id2", "chromtrans", "starttrans", "endtrans", "transcoord", "intronanno", "chrTE", "startTE", "endTE", "subfamTE", "numTE", "strandTE", "strand","covtrans", "covexon")
  
  colnames(filtered_table) <- columnlabels
  
  filtered_table$file_label <- rep(strsplit(file_name, split = "_")[[1]][1], length(filtered_table$strand))
  
  file_name2 <- gsub("filtered_","",file_name)
  
  annotated_table = read.delim(file_name2, sep = "\t", stringsAsFactors = FALSE, header = FALSE)
  colnames(annotated_table) <- columnlabels
  
  maxtranscovfun <- function(tablein){
    return(max(tablein$covtrans))
  }
  maxexoncovfun <- function(tablein){
    return(max(tablein$covexon))
  }
  
  filtered_table$maxtranscov <- apply(filtered_table[,c('gene2','chr2')], 1, function(x) maxtranscovfun(annotated_table[annotated_table$gene2 == x[1],]))
  filtered_table$maxexoncov <- apply(filtered_table[,c('gene2','chr2')], 1, function(x) maxexoncovfun(annotated_table[annotated_table$gene2 == x[1],]))
  
  if (grepl("N",file_name)){
    filtered_table$tumor_count <- rep(0, nrow(filtered_table))
    filtered_table$normal_count <- rep(1, nrow(filtered_table))
  } else {
    filtered_table$tumor_count <- rep(1, nrow(filtered_table))
    filtered_table$normal_count <- rep(0, nrow(filtered_table))
  }
  
  
  if (firsttime == TRUE){
    combined_table = filtered_table
    firsttime = FALSE
  } else {
    combined_table <- rbind(combined_table,filtered_table,stringsAsFactors = FALSE)
  }
  
  print(i)
  i = i + 1
}

saveRDS(filter_combined_table, 'filtertablefinal06202018.rds')
#This is a filter to make sure that candidate transcripts represent a certain percentage of the overall transcripts for a gene
pfilter <- .25
filter_combined_table = combined_table[combined_table$covtrans > pfilter*combined_table$maxtranscov,]
#Remove those that start in exons so its just intronic and intergenic
TEreftable <- read.delim("/bar/nshah/programs/Pyscript/TE.lst", skip=1, header=FALSE, stringsAsFactors = FALSE)

indexsubfams <- match(filter_combined_table$subfamTE, TEreftable$V1)
TEclass_v1 <- TEreftable$V2[indexsubfams]
TEfamily_v1 <- TEreftable$V3[indexsubfams]

filter_combined_table$classTE <- TEclass_v1
filter_combined_table$familyTE <- TEfamily_v1
filter_combined_table <- filter_combined_table[filter_combined_table$exonintron1 != "exon",]

filter_combined_table$uniqidfile <- paste(filter_combined_table$subfamTE,filter_combined_table$startTE,filter_combined_table$gene1,filter_combined_table$exonintron1, filter_combined_table$number1, filter_combined_table$gene2, filter_combined_table$exonintron2, filter_combined_table$number2,filter_combined_table$transcriptstart2, filter_combined_table$file_label, sep = "_")
freqmultipletable <- as.data.frame(table(filter_combined_table$uniqidfile))
i <- sapply(freqmultipletable, is.factor)
freqmultipletable[i] <- lapply(freqmultipletable[i], as.character)
freqmultipletable <- freqmultipletable[freqmultipletable$Freq >1,]

#Make sure that if there are multiple isoforms per candidate they are not counted multiple times
exon1length <- c()
intron1start <- c()
intron1end <- c()
for (transcoords in filter_combined_table$transcoord){
  elements <- strsplit(transcoords, ",")[[1]]
  strandtrans <- elements[1]
  if (strandtrans == "+"){
    sortedvec <- sort(as.numeric(elements[3:length(elements)]))
    exon1length <- c(exon1length, sortedvec[2] - sortedvec[1] + 1)
    intron1start <- c(intron1start, sortedvec[2] + 1)
    intron1end <- c(intron1end, sortedvec[3] - 1)
  } else {
    sortedvec <- sort(as.numeric(elements[3:length(elements)]), decreasing = TRUE)
    exon1length <- c(exon1length, sortedvec[1] - sortedvec[2] + 1)
    intron1start <- c(intron1start, sortedvec[3] + 1)
    intron1end <- c(intron1end, sortedvec[2] - 1)
  }
}

filter_combined_table$exon1length <- exon1length
filter_combined_table$intron1start <- intron1start
filter_combined_table$intron1end <- intron1end

#Make sure that if there are multiple isoforms per candidate they are not counted multiple times
for (multipleid in freqmultipletable$Var1){
  locations_with <- which(filter_combined_table$uniqidfile %in% multipleid)
  filter_combined_table$normal_count[locations_with[2:length(locations_with)]] <- 0
  filter_combined_table$tumor_count[locations_with[2:length(locations_with)]] <- 0
  filter_combined_table$exon1length[locations_with[2:length(locations_with)]] <- 0
}

counts_table = filter_combined_table[,c('subfamTE', 'classTE', 'familyTE', 'chrTE', 'startTE', 'endTE', 'gene1', 'gene2', 'exonintron1', 'number1', 'exonintron2', 'number2', 'transcriptstart2', 'strand','exon1length', 'tumor_count', 'normal_count')]
resexonlength <- aggregate(. ~ subfamTE + classTE + familyTE + chrTE + startTE + endTE + gene1 + gene2 + exonintron1 + number1 + exonintron2 + number2 + transcriptstart2 + strand, counts_table, FUN = sum)

i <- sapply(resexonlength, is.factor)
resexonlength[i] <- lapply(resexonlength[i], as.character)

#resexonlength$exonlengthaverage <- as.numeric(resexonlength$exon1length)/(as.numeric(resexonlength$normalcount) + as.numeric(resexonlength$tumorcount))

exonskipped <- c()
for (i in 1:length(resexonlength$number1)){
  if (resexonlength$gene1[i] != resexonlength$gene2[i]){
    exonskipped <- c(exonskipped, -1)
  } else {
    exonskipped <- c(exonskipped, as.numeric(resexonlength$number2[i])-as.numeric(resexonlength$number1[i]))
  }
}

resexonlength$exonskipped <- exonskipped
resexonlength$ExonlengthAverage <- resexonlength$exon1length/(resexonlength$tumor_count + resexonlength$normal_count)

#Filters we applied in out paper and get rid of noise
# 1. At least 4 tumors have it
# 2. No normal samples have it
# 3. No more than 2 exons skipped
# 4. Not longer than the 99th percentile of exon1s in gencode

resexonlength <- resexonlength[resexonlength$tumor_count > 3,]
resexonlength <- resexonlength[resexonlength$normal_count == 0,]
resexonlength <- resexonlength[resexonlength$exonskipped <3,]
resexonlength <- resexonlength[resexonlength$ExonlengthAverage < 2588,]

#Filter down the data further so that 
resexonlength$uniqid <- paste(resexonlength$subfamTE,resexonlength$startTE,resexonlength$gene1,resexonlength$exonintron1, resexonlength$number1, resexonlength$gene2, resexonlength$exonintron2, resexonlength$number2,resexonlength$transcriptstart2,sep = "_")
filter_combined_table$uniqid <- paste(filter_combined_table$subfamTE,filter_combined_table$startTE,filter_combined_table$gene1,filter_combined_table$exonintron1, filter_combined_table$number1, filter_combined_table$gene2, filter_combined_table$exonintron2, filter_combined_table$number2,filter_combined_table$transcriptstart2,sep = "_")
filter_combined_table <- filter_combined_table[filter_combined_table$uniqid %in% resexonlength$uniqid,]

filter_combined_table$file_label <- gsub( ".gtf", "", as.character(filter_combined_table$file_label))
write.table(filter_combined_table, "filter_combined_candidates.tsv", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE) 

#PAUSE. After this step I need to run some code with filter_combined_candidates.tsv to extract the bam reads

# Load the stats
readstats <- read.delim('filter_combined_candidates_stats.tsv', sep = "\t", stringsAsFactors = FALSE, header = FALSE)
colnames(readstats) <- c('uniqid', 'file_name', 'read', 'startread', 'endread', 'filetype')
readstats$uniqidfile <- paste0(readstats$uniqid, '_', readstats$file_name)

filter_combined_table$uniqidfile <- gsub( ".gtf", "", as.character(filter_combined_table$uniqidfile))

filter_combined_table <- filter_combined_table[!duplicated(filter_combined_table$uniqidfile),]
readstats <- readstats[!duplicated(readstats$uniqidfile),]

filter_combined_table_stat <- merge(filter_combined_table,readstats, by = c('uniqidfile'), all.x = TRUE, all.y = FALSE)
filter_combined_table_stat$uniqidfile <- NULL
filter_combined_table_stat$uniqid.y <- NULL
names(filter_combined_table_stat)[names(filter_combined_table_stat) == 'uniqid.x'] <- 'uniqid'

stataggregate <- aggregate(filter_combined_table_stat[,c("read",'startread', 'endread')], by=list(filter_combined_table_stat[,c("uniqid")]), FUN=sum, na.rm=TRUE)
names(stataggregate)[names(stataggregate) == 'Group.1'] <- 'uniqid'

resexonlengthstat <- merge(resexonlength, stataggregate, by=c('uniqid'), all.x = TRUE, all.y = FALSE)
res2stat <- resexonlengthstat

startreadper <- c()
endreadper <- c()
seper <- c()
distanceTE <- c()
for (idcandidate in res2stat$uniqid){
  testtable <- filter_combined_table_stat[filter_combined_table_stat$uniqid == idcandidate, ]
  totalsamples <- nrow(testtable)
  totalsamplesstart <- sum(testtable$startread > 0)
  totalsamplesend <- sum(testtable$endread > 0)
  sesamples <- sum(testtable$filetype == 'se')
  startreadper <- c(startreadper, totalsamplesstart/totalsamples)
  endreadper <- c(endreadper, totalsamplesend/totalsamples)
  seper <- c(seper, sesamples/totalsamples)
  
  row1 <- testtable[1,]
  strand1 <- row1$strand
  if (strand1 == '+'){
    distanceTE <- c(distanceTE, row1$transcriptstart2 - row1$endTE)
  } else {
    distanceTE <- c(distanceTE, row1$startTE - row1$transcriptend2)
  }
}

res2stat$endreadper <- endreadper
res2stat$startreadper <- startreadper
res2stat$seper <- seper
res2stat$distanceTE <- distanceTE

res2statfil <- res2stat[res2stat$read > 0, ]
res2statfil <- res2statfil[res2statfil$endreadper <= .1, ]
res2statfil <- res2statfil[res2statfil$gene1 == 'None' | res2statfil$startread >0,]
res2statfil <- res2statfil[res2statfil$seper < 1.0,]
res2statfil <- res2statfil[res2statfil$distanceTE >2500 | res2statfil$distanceTE < 0,]

#Fantom5 Annotation
df <- data.frame(chr=res2statfil$chrTE,start=res2statfil$startTE,end=res2statfil$endTE)

library(rtracklayer)
library(GenomicRanges)

chain <- import.chain("~/reference/hg38ToHg19.over.chain")
chrTE19 <- c()
startTE19 <- c()
endTE19 <- c()
for (i in 1:length(df$chr)){
  
  rowcandidate <- df[i,]
  gr <- makeGRangesFromDataFrame(rowcandidate, TRUE)
  hg19coordinate <- liftOver(gr, chain)
  result <- as.data.frame(hg19coordinate, stringsAsFactors = FALSE)
  if (length(result$start)!=0){
    chrTE19 <- c(chrTE19, as.character(result$seqnames[1]))
    startTE19 <- c(startTE19, result$start[1])
    endTE19 <- c(endTE19, result$end[1])
  } else {
    chrTE19 <- c(chrTE19, "None")
    startTE19 <- c(startTE19, "None")
    endTE19 <- c(endTE19, "None")
  }
}

# Before doing the rest of this code you should run liftover on the coordinates (liftOver TEdf.bed ~/reference/hg38ToHg19.over.chain TEdf_hg19.bed unMapped)
#dfreplace <- read.table("TEdf_hg19.bed",stringsAsFactors = FALSE, header = FALSE)

res2statfil$hg19chrTE <- chrTE19
res2statfil$hg19startTE <- startTE19
res2statfil$hg19endTE <- endTE19

#readRDS
fantom5 <- readRDS("/bar/nshah/reference/fantom5.RDS")

fantom5anno <- function(chromosome, start, end, strand, fantom5){
  if (chromosome!="None"){
    searchspace <- fantom5[fantom5$strand == strand & fantom5$chromosome == chromosome, ]
    anyoverlap <- apply(searchspace[,c('start','end')], 1, function(x) as.numeric(max(0,min(end,x[2])-max(start,x[1])))/(as.numeric(x[2]-x[1])))
    candidatepeaks <- searchspace[anyoverlap>.5,]
    candidatepeaks$chromosome <- NULL
    candidatepeaks$start <- NULL
    candidatepeaks$end <- NULL
    candidatepeaks$strand <- NULL
    
    tpmlevels <- data.frame(sumtpms = colSums(candidatepeaks), labels = colnames(candidatepeaks), stringsAsFactors = FALSE)
    tpmlevels <- tpmlevels[tpmlevels$sumtpms > .01,]
    tpmlevels <- tpmlevels[order(-tpmlevels$sumtpms),]
    
    return(c(length(tpmlevels$sumtpms),paste(as.character(tpmlevels$sumtpms),collapse=","),paste(tpmlevels$labels,collapse=",")))
  } else {
    return(c("None","None","None"))
  }
}

resultsof <- apply(res2statfil[,c('hg19chrTE','hg19startTE','hg19endTE','strand')], 1, function(x) fantom5anno(x[1],as.numeric(x[2]),as.numeric(x[3]),x[4],fantom5))
resultsof <- data.frame(t(resultsof))
colnames(resultsof) <- c("Number_Samples_with_Peak", "Tpm_Samples", "Name_Samples")

res2final <- cbind(res2statfil,resultsof)

library('xlsx')
write.xlsx(res2final, '062018_prostate-cancer-filtered-candidates.xlsx', row.names = FALSE)

#GFF3 Creation
chrxamples <- c()
soucexamples <- c()
featurexamples <- c()
startxamples <- c()
endxamples <- c()
scorexamples <- c()
strandxamples <- c()
framexamples <- c()
groupxamples <- c()
for(i in 1:nrow(res2final)) {
  row1 <- res2final[i,]
  theexamples <- filter_combined_table_stat[filter_combined_table_stat$uniqid == row1$uniqid,]
  uniqid <- paste(row1$subfamTE,row1$startTE,row1$gene1,row1$exonintron1, row1$number1, row1$gene2, row1$exonintron2, row1$number2, row1$exonintron1start, row1$intron1end,sep = "_")
  theexamples <- theexamples[order(-theexamples$covtrans),]
  examplestructure <- theexamples[1,]
  structurestring <- examplestructure$transcoord[1]
  chrexample <- strsplit(structurestring, ",")[[1]][2]
  strandexample <- strsplit(structurestring, ",")[[1]][1]
  elementsvec <- as.numeric(tail(strsplit(structurestring, ",")[[1]],-2))
  elementsvec <- sort(elementsvec)
  numexons <- length(elementsvec)/2
  chrxamples <- c(chrxamples, chrexample)
  soucexamples <- c(soucexamples, "curated")
  featurexamples <- c(featurexamples, "mRNA")
  startxamples <- c(startxamples, min(elementsvec))
  endxamples <- c(endxamples, max(elementsvec))
  scorexamples <- c(scorexamples, ".")
  strandxamples <- c(strandxamples, strandexample)
  framexamples <- c(framexamples, ".")
  groupxamples <- c(groupxamples, paste0('ID=',uniqid))
  for (j in 1:numexons){
    chrxamples <- c(chrxamples, chrexample)
    soucexamples <- c(soucexamples, "curated")
    featurexamples <- c(featurexamples, "exon")
    startxamples <- c(startxamples, elementsvec[j*2-1])
    endxamples <- c(endxamples, elementsvec[j*2])
    scorexamples <- c(scorexamples, ".")
    strandxamples <- c(strandxamples, strandexample)
    framexamples <- c(framexamples, ".")
    groupxamples <- c(groupxamples, paste0('Parent=',uniqid))
  }
}

gffexamples <- data.frame(chrxamples, soucexamples, featurexamples, startxamples, endxamples, scorexamples, strandxamples, framexamples, groupxamples)
write.table(gffexamples, "062018_pipelineexamples.gff3", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

#Transcript structure selection
# Currently, in the gff3 file I do this using 

for(i in 1:nrow(res2final)) {
  row1 <- res2final[i,]
  theexamples <- filter_combined_table_stat[filter_combined_table_stat$uniqid == row1$uniqid,]
  uniqid <- row1$uniqid
  theexamples <- theexamples[order(-theexamples$covtrans),]
  examplestructure <- theexamples[1,]
  structurestring <- examplestructure$transcoord[1]
  chrexample <- strsplit(structurestring, ",")[[1]][2]
  strandexample <- strsplit(structurestring, ",")[[1]][1]
  elementsvec <- as.numeric(tail(strsplit(structurestring, ",")[[1]],-2))
  elementsvec <- sort(elementsvec)
  numexons <- length(elementsvec)/2
  chrxamples <- c(chrxamples, chrexample)
  soucexamples <- c(soucexamples, "curated")
  featurexamples <- c(featurexamples, "mRNA")
  startxamples <- c(startxamples, min(elementsvec))
  endxamples <- c(endxamples, max(elementsvec))
  scorexamples <- c(scorexamples, ".")
  strandxamples <- c(strandxamples, strandexample)
  framexamples <- c(framexamples, ".")
  groupxamples <- c(groupxamples, paste0('ID=',uniqid))
  for (j in 1:numexons){
    chrxamples <- c(chrxamples, chrexample)
    soucexamples <- c(soucexamples, "curated")
    featurexamples <- c(featurexamples, "exon")
    startxamples <- c(startxamples, elementsvec[j*2-1])
    endxamples <- c(endxamples, elementsvec[j*2])
    scorexamples <- c(scorexamples, ".")
    strandxamples <- c(strandxamples, strandexample)
    framexamples <- c(framexamples, ".")
    groupxamples <- c(groupxamples, paste0('Parent=',uniqid))
  }
}


