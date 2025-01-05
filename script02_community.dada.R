rm(list=ls()) #cleans up working space

library(microbiome)
library(dada2)
library(phyloseq)
library(DECIPHER)
library(plyr)
library(tidyverse)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq

print('Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq')
path<-("/PATH/TO/SEQUENCE/DATA")
list.files(path)

fnFs <- sort(list.files(path, pattern=glob2rx("*_R1.fastq"), full.names = TRUE))
fnRs <- sort(list.files(path, pattern=glob2rx("*_R2.fastq"), full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
print('Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq and create dataframe with metadata')
sample.names <- substring(sapply(strsplit(basename(fnFs), "_"), `[`, 1),15,20)


# Place filtered files in filtered/ subdirectory
print('Place filtered files in filtered/ subdirectory')
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220,190),orient.fwd ="TACG",
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

# Learn the error rates of each base
print('Learn the error rates of each base')
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Remove replications
print('Remove replications')
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
print('Name the derep-class objects by the sample names')
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Sample inference
print('Sample inference')
dadaFs <- dada(derepFs, err=errF, multithread=TRUE,pool=T)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE,pool=T)

dadaFs[[1]]
dadaRs[[1]]

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
print('Inspect the merger data.frame from the first sample')
head(mergers[[1]])

# Construct sequence table
print('Construct sequence table')
seqtab <- makeSequenceTable(mergers)
dim(seqtab)


# Inspect distribution of sequence lengths
print('Construct sequence table')
table(nchar(getSequences(seqtab)))

# Remove chimeras
print('remove chimeras')
seqtab.nochim <- removeBimeraDenovo(seqtab, method="pooled", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
print(paste (sum(seqtab.nochim)/sum(seqtab),': fraction of recovered reads after removal of chimeras'))


# Track the number of reads that made it through so far
print('Track the number of reads that made it through so far')
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track,file="/PATH/TO/OUT.DIRECTORY/dada.Number_of_reads_per_sample.tab", sep = "\t")

# Assign taxonomy.
print('Assign taxonomy')
taxa.gt <- assignTaxonomy(seqtab.nochim, "GTDB_bac-arc_ssu_r86.fa.gz",  multithread=TRUE) 

# Inspect taxonomic assignments
print('Inspect taxonomic assignments')
taxa.print <- taxa.gt # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
taxa.print.names<-as.data.frame(taxa.print)


# Extract of sequences
print('Extract of sequences')
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs
names(seqs)<- paste0("SV_", seq(ntaxa(ps0)), "_", taxa.print.names$Order)
seqs.tab<-as.data.frame(seqs)


# Simple data.frame construction from the information encoded in the filenames
print('Simple data.frame construction from the information encoded in the filenames')
seqtab.nochim2<-seqtab.nochim
colnames(seqtab.nochim2)<-rownames(seqs.tab)
rownames(taxa.gt)<-rownames(seqs.tab)
samples.out2 <- rownames(seqtab.nochim)
samdf <- data.frame(sample.ID=sample.names)
rownames(samdf) <- samples.out2

# Create phyloseq object with the data with the new sample names.
print('phyloseq object with the data with the new sample names')
ps <- phyloseq(otu_table(seqtab.nochim2, taxa_are_rows=FALSE),
               sample_data(samdf),
               tax_table(taxa.gt))


# Extract count table (raw count data)
print('extract count table (raw count data)')
OTU1 = as(otu_table(ps), "matrix") #temporary file
# transpose if necessary
if(taxa_are_rows(ps)){OTU1 <- t(OTU1)}
# Coerce to data.frame
OTUdf.ps = as.data.frame(OTU1)
OTUdf.ps <- OTUdf.ps[with(OTUdf.ps, order(row.names(OTUdf.ps))), ]


# Keep samples that are relevant for publication
print('Keep samples that are relevant for publication')
ps <-prune_samples(sample_names(ps)[c(1:13,17:25)],ps)  #ps: phyloseq object with raw output from dada2

# Extract count table (raw count data after sample removal)
print('Extract count table (raw count data)')
OTU1 = as(otu_table(ps), "matrix") #temporary file
if(taxa_are_rows(ps)){OTU1 <- t(OTU1)} # transpose if necessary
OTUdf.ps = as.data.frame(OTU1) # Coerce to data.frame

# Rarefy data
print('Rarefy data')
OTUdf.rar.ps <-rrarefy(OTUdf.ps, min (rowSums(OTUdf.ps)))

# Metadata
print('Create dataframe for metadata')
source <- c(rep('C',9),rep('CS',4),rep('S',9))
DOM <- c(rep('org',3), rep('SW-DOM',3), rep ('S-DOM',3), rep('SW-DOM',2), rep('S-DOM',2),rep('SW-DOM',3), rep('org',3), rep ('S-DOM',3))
STR <- data.frame(source, DOM)
rownames(STR) <- sample_names(ps)
STR$treat <- paste (STR$source, STR$DOM, sep='.')
STR

#  Create phyloseq object with rarefied relative count data
print('Create phyloseq object with rarefied relative count data')
ps.rar <- phyloseq(otu_table(OTUdf.rar.ps, taxa_are_rows=FALSE),
                   sample_data(STR),
                   tax_table(taxa.gt))
ps.rar

# Create relative count data from rarefied dataset
print('Create relative count data from rarefied dataset')
norm_rel = function(x) (x / sum(x)) #function for renormalisation after subsampling
ps.rel= transform_sample_counts(ps.rar, norm_rel ) #transform to relative count data
sample_sums(ps.rel)

# Create counttable for PICRUSt2 input
print('Create counttable for PICRUSt2 input')
rOTU <- as.data.frame(t(otu_table(ps.rel))) %>%
  rownames_to_column(var="ASV")

# Data export
print('Data export')
# Export count table (rarefied, relative, rows: ASVs) >> PICRUSt input
write.table(rOTU,"/PATH/TO/OUT.DIRECTORY/dada-cryo.rcounts.ASV.2.tab", row.names=FALSE,quote=FALSE, sep='\t')
# Export ASV fasta file >> PICRUSt input
write.table(seqs, "/PATH/TO/OUT.DIRECTORY/dada-cryo.seqs.nochim.fasta",quote=FALSE, col.names=FALSE)
# Export Rdata
save.image("/PATH/TO/OUT.DIRECTORY/dada-coal.img")



