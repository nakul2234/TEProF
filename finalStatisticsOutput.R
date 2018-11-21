#!/usr/bin/env Rscript
library('Xmisc')
library('reshape2')

load('Step10.RData')

#Get the directory of the rnapipeline scripts
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

#Default parameters
minGeneFPKM <- 1
minFracFPKM <- .1
minIntronRead <- 1
parseTreatment <- ''
argumentFile <- paste0(script.basename,'/arguments.txt')

#User argument parsing code
parser <- ArgumentParser$new()

parser$add_argument('-e',type='character',help='A character of parsing the treatment samples')
parser$add_argument('-i',type='character',help='Minimum reads that have to support the splice junction into the gene')
parser$add_argument('-g',type='character',help='Minimum total FPKM of gene that candidate splices into')
parser$add_argument('-f',type='character',help='Minimum fraction of total FPKM of gene of candidate')
parser$add_argument('-a',type='character',help='argument.txt file being used for this run')
argsparse <- parser$get_args()

#Argument parse to change defaults
if (!identical(argsparse$e, character(0))){
  parseTreatment <- argsparse$e
}

if (!identical(argsparse$i, character(0))){
  minIntronRead <- as.numeric(argsparse$i)
}

if (!identical(argsparse$g, character(0))){
  minGeneFPKM <- as.numeric(argsparse$g)
}

if (!identical(argsparse$f, character(0))){
  minFracFPKM <- as.numeric(argsparse$f)
}

if (!identical(argsparse$a, character(0))){
  argumentFile <- argsparse$a
}

print("Arguments")
print(paste0("Treatment Parse: ", parseTreatment))
print(paste0("Minimum Intron Read for Presence: ", minIntronRead))
print(paste0("Minimum Gene FPKM: ", minGeneFPKM))
print(paste0("Minimum Fraction of Gene FPKM: ", minFracFPKM))
print(paste0("ArgumentFile: ", argumentFile))

fracexpressiontable <- read.delim("table_frac_tot_cand",header=FALSE, stringsAsFactors = FALSE)
fracexpressiontable <- as.data.frame(t(fracexpressiontable))
fracexpressiontable <- data.frame(lapply(fracexpressiontable, as.character), stringsAsFactors=FALSE)
colnames(fracexpressiontable) <- as.character(fracexpressiontable[1,])
fracexpressiontable <- fracexpressiontable[c(-1),]
fracexpressiontable[-1] <- lapply(fracexpressiontable[-1], as.numeric)

totexpressiontable <- read.delim("table_tot_cand",header=FALSE, stringsAsFactors = FALSE)
totexpressiontable <- as.data.frame(t(totexpressiontable))
totexpressiontable<- data.frame(lapply(totexpressiontable, as.character), stringsAsFactors=FALSE)
colnames(totexpressiontable) <- as.character(totexpressiontable[1,])
totexpressiontable <- totexpressiontable[c(-1),]
totexpressiontable[-1] <- lapply(totexpressiontable[-1], as.numeric)

fracexpressiontable.m <- melt(fracexpressiontable, id.vars = c('TranscriptID'))
totexpressiontable.m <- melt(totexpressiontable, id.vars = c('TranscriptID'))

if (!(identical(fracexpressiontable.m[,c(1)], totexpressiontable.m[,c(1)]))){
  stop(call = TRUE)
} #Proof that the order of the samples is exactly the same and thus instead of an expensive merge I can do a simple rbind. If this is not the case then something has gone wrong with table creation

colnames(fracexpressiontable.m)[3] <- "fractotal"
fracexpressiontable.m$totexpression <- totexpressiontable.m$value

annotatedcufftranscripts <- annotatedcufftranscripts[order(annotatedcufftranscripts$uniqid),]
annotatedcufftranscripts$duplicated <- duplicated(annotatedcufftranscripts$uniqid)

oldtranscripts <- c()
newtranscripts <- c()

#Get a list of the transcripts

currenttranscript = "None"
for (i in 1:nrow(annotatedcufftranscripts)){
  uniqidrow = annotatedcufftranscripts$uniqid[i]
  transcriptrow = annotatedcufftranscripts$transcriptname[i]
  duplicatedsample = annotatedcufftranscripts$duplicated[i]
  if (duplicatedsample == TRUE){
    oldtranscripts <- c(oldtranscripts, transcriptrow )
    newtranscripts <- c(newtranscripts, currenttranscript)
    
  } else {
    currenttranscript = transcriptrow
  }
}

fracexpressiontable.m$variable <- as.character(fracexpressiontable.m$variable)
j = 1
for (oldtranscript in oldtranscripts){
  fracexpressiontable.m$variable[fracexpressiontable.m$variable==oldtranscript] <- newtranscripts[j]
  j=j+1
}

#Remove duplicate unique ids that could be tripping up the program. Witin the same sample they will be added together since they are pretty much the same just a downstream exon could be causing the difference.
fracexpressiontable.m <- aggregate(. ~ TranscriptID + variable + totexpression,fracexpressiontable.m,FUN = sum)

#Incorporate intron annotation information
intronexpressiontable <- read.delim("table_i_all",header=FALSE, stringsAsFactors = FALSE)
intronexpressiontable <- as.data.frame(t(intronexpressiontable))
intronexpressiontable <- data.frame(lapply(intronexpressiontable, as.character), stringsAsFactors=FALSE)

intronlabels <- paste(intronexpressiontable[1,2:ncol(intronexpressiontable)], intronexpressiontable[2,2:ncol(intronexpressiontable)], intronexpressiontable[3,2:ncol(intronexpressiontable)], intronexpressiontable[4,2:ncol(intronexpressiontable)], sep = "_")
intronlabels <- c("TranscriptID",intronlabels)

colnames(intronexpressiontable) <- intronlabels

intronexpressiontable <- intronexpressiontable[c(-1,-2,-3,-4),]

intronexpressiontable[-1] <- lapply(intronexpressiontable[-1], as.numeric)

intronexpressiontable.m <- melt(intronexpressiontable, id.vars = c('TranscriptID'))

annotatedcufftranscripts$intronname <- paste(annotatedcufftranscripts$chrTE,annotatedcufftranscripts$strand,annotatedcufftranscripts$intronjunstart,annotatedcufftranscripts$intronjunend, sep = "_")

indexfileintron1 <- match(fracexpressiontable.m$variable, annotatedcufftranscripts$transcriptname)

fracexpressiontable.m$intronlabel <- annotatedcufftranscripts$intronname[indexfileintron1]

fracexpressiontable.m$fileintronlabel <- paste(fracexpressiontable.m$TranscriptID, fracexpressiontable.m$intronlabel, sep = '_')
intronexpressiontable.m$fileintronlabel <- paste(intronexpressiontable.m$TranscriptID, intronexpressiontable.m$variable, sep = '_')

indexfileintron <- match(fracexpressiontable.m$fileintronlabel, intronexpressiontable.m$fileintronlabel)

fracexpressiontable.m$intronread <- intronexpressiontable.m$value[indexfileintron]

fracexpressiontable.m.exp.t <- fracexpressiontable.m[fracexpressiontable.m$fractotal >= minFracFPKM & fracexpressiontable.m$totexpression >= minGeneFPKM & fracexpressiontable.m$intronread >= minIntronRead & grepl(parseTreatment,fracexpressiontable.m$TranscriptID), ]
fracexpressiontable.m.exp.n <- fracexpressiontable.m[fracexpressiontable.m$fractotal >= minFracFPKM & fracexpressiontable.m$totexpression >= minGeneFPKM & fracexpressiontable.m$intronread >= minIntronRead & !grepl(parseTreatment,fracexpressiontable.m$TranscriptID), ]

annotatedcufftranscripts$tumor_count <- apply(annotatedcufftranscripts[,c('transcriptname','gene2')],1,function(x) sum(fracexpressiontable.m.exp.t$variable==x[1]))

annotatedcufftranscripts$normal_count <- apply(annotatedcufftranscripts[,c('transcriptname','gene2')],1,function(x) sum(fracexpressiontable.m.exp.n$variable==x[1]))

annotatedcufftranscripts$tumor_fracmean <- apply(annotatedcufftranscripts[,c('transcriptname','gene2')],1,function(x) mean(fracexpressiontable.m[fracexpressiontable.m$variable==x[1] & grepl(parseTreatment,fracexpressiontable.m$TranscriptID), c('fractotal')]))

annotatedcufftranscripts$normal_fracmean <- apply(annotatedcufftranscripts[,c('transcriptname','gene2')],1,function(x) mean(fracexpressiontable.m[fracexpressiontable.m$variable==x[1] & !grepl(parseTreatment,fracexpressiontable.m$TranscriptID), c('fractotal')]))

annotatedcufftranscripts$tumor_genefpkm_mean <- apply(annotatedcufftranscripts[,c('transcriptname','gene2')],1,function(x) mean(fracexpressiontable.m[fracexpressiontable.m$variable==x[1] & grepl(parseTreatment,fracexpressiontable.m$TranscriptID), c('totexpression')]))

annotatedcufftranscripts$normal_genefpkm_mean <- apply(annotatedcufftranscripts[,c('transcriptname','gene2')],1,function(x) mean(fracexpressiontable.m[fracexpressiontable.m$variable==x[1] & !grepl(parseTreatment,fracexpressiontable.m$TranscriptID), c('totexpression')]))

annotatedcufftranscripts$tumor_intronjuncount_mean <- apply(annotatedcufftranscripts[,c('transcriptname','gene2')],1,function(x) mean(fracexpressiontable.m[fracexpressiontable.m$variable==x[1] & grepl(parseTreatment,fracexpressiontable.m$TranscriptID), c('intronread')]))

annotatedcufftranscripts$normal_intronjuncount_mean <- apply(annotatedcufftranscripts[,c('transcriptname','gene2')],1,function(x) mean(fracexpressiontable.m[fracexpressiontable.m$variable==x[1] & !grepl(parseTreatment,fracexpressiontable.m$TranscriptID), c('intronread')]))

annotatedcufftranscripts <- annotatedcufftranscripts[annotatedcufftranscripts$duplicated == FALSE,]

numberTumorSamples <- length(unique(fracexpressiontable.m$TranscriptID[grepl(parseTreatment, fracexpressiontable.m$TranscriptID)]))
numberNormalSamples <- length(unique(fracexpressiontable.m$TranscriptID[!grepl(parseTreatment, fracexpressiontable.m$TranscriptID)]))

annotatedcufftranscripts$tumor_enrichment <- ((annotatedcufftranscripts$tumor_count)/numberTumorSamples)/((annotatedcufftranscripts$normal_count)/numberNormalSamples)

TEreftable <- read.delim(repeatAnnotationFile, skip=1, header=FALSE, stringsAsFactors = FALSE)

indexsubfams <- match(annotatedcufftranscripts$subfamTE, TEreftable$V1)
TEclass_v1 <- TEreftable$V2[indexsubfams]
TEfamily_v1 <- TEreftable$V3[indexsubfams]

annotatedcufftranscripts$classTE <- TEclass_v1
annotatedcufftranscripts$familyTE <- TEfamily_v1


#Final output tables for analysis

##The Candidates that are present in at least one sample

exportTable2 <- annotatedcufftranscripts[,c("subfamTE", "classTE", "familyTE","chrTE", "startTE", "endTE", "exonintron1", "number1", "gene2", "exonintron2", "number2","strand","tumor_count", "normal_count", "tumor_fracmean", "normal_fracmean", "tumor_genefpkm_mean", "normal_genefpkm_mean", "tumor_intronjuncount_mean", "normal_intronjuncount_mean", "tumor_enrichment")]

exportTable2$LocationTE <- paste0(exportTable2$exonintron1,'_',exportTable2$number1)
exportTable2$LocationTE[exportTable2$LocationTE == "None_None"] <- "Intergenic"
exportTable2$SpliceTarget <- paste0(exportTable2$exonintron2,'_',exportTable2$number2)

exportTable2 <- exportTable2[,c("subfamTE", "classTE", "familyTE","chrTE", "startTE", "endTE", "LocationTE", "gene2", "SpliceTarget","strand","tumor_count", "normal_count", "tumor_fracmean", "normal_fracmean", "tumor_genefpkm_mean", "normal_genefpkm_mean", "tumor_intronjuncount_mean", "normal_intronjuncount_mean", "tumor_enrichment")]

colnames(exportTable2) <- c('Subfam','Class',	'Family',	'Chr TE',	'Start TE',	'End TE',	'Location TE', 'Gene', 'Splice Target',	'Strand',	'Treatment Total',	'Normal Total', 'Fraction of Total Expression (Treatment)','Fraction of Total Expression (Normal)', 'Mean Gene Expression (Treatment)','Mean Gene Expression (Normal)','Mean Intron Read Count (Treatment)','Mean Intron Read Count (Normal)', 'Treatment Sample Enrichment')

library(xlsx)
exportTable2 <- exportTable2[order(-exportTable2$`Treatment Total`),]
write.xlsx(exportTable2, "All TE-derived Alternative Isoforms Statistics.xlsx", row.names = FALSE)

##Expression, fraction expression, and intron coverage information across all of the files.
exportAllStats <- fracexpressiontable.m[,c('TranscriptID', "variable", "totexpression", "fractotal", "intronread")]
colnames(exportAllStats) <- c('File', 'Transcript_Name', 'Total Gene Expression (FPKM)', 'Fraction of Total Gene Expression (FPKM)', 'Intron Read Count')
write.table(exportAllStats, "allCandidateStatistics.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

##Final image file 
save.image('Step11_FINAL.RData')

