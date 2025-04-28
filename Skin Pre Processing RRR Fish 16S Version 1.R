
### Script N. 1 - Fish skin microbiome data pre-processing ###

## Most recent edit: October 23, 2023

#### Code for reading and cleaning 16S sequences in R ####

## Base code adapted from dada2 pipeline tutorial by Matthieu Leray (ML)
## Edited by author; September/October 2022 -- McGill & STRI

## Code re-run January 2023 after fixing cutadapt error in trimming step

#### Read-Me ####
# This code takes trimmed sequencing reads from Cutadapt (done in Terminal), filters, trims, and dereplicates reads.
# Then, ASVs are inferred (code for 3 possible methods; used pseudo-pooling), paired reads are merged, 
# and a sequence table is made. Chimeras are removed, taxonomy is assigned using the SILVA dataset (v 138.1)
# and reads with taxonomy assigned is saved as a phyloseq object. Contaminants are identified and cleaned,
# then control samples are removed from the dataset. We remove reads assigned to chloroplasts, mitochondria, and eukaryotes.
# The final cleaned dataset is saved as a new phyloseq object (ROHR_06_obj_clean).

#### RRR Fish Microbiome - TEP Skin Swabs ####

## Project summary ##
# Fish were collected from 2 sites in the Tropical Eastern Pacific: Las Perlas & Coiba
# Sampling during both upwelling & non-upwelling seasons (collections done 2021-22)
# 10 target species
# Swabs = 9 per species
# Control swabs (swabs exposed to air in field) collected at each site and season
# Swabs were stored in liquid nitrogen for transport from field to lab, then at -80ºC
# DNA was extracted using ZymoBIOMICS 96 well DNA kit - 4 plates total
# DNA was amplified following a modified version of the Earth Microbiome Project
# 16S protocol using phased primers (1 set per plate)
# Negative PCR controls and extraction controls included in each plate
# Index PCR performed to attach unique barcodes, then all four plates (384 samples)
# pooled for sequencing on the Illumina MiSeq sequencing platform of the Smithsonian
# Tropical Research Institute's (STRI) Naos facilities in Panama (courtesy of Marta Vargas)
# Sequenced libraries were demultiplexed at STRI by Marta Vargas using the MiSeq Reporter Software
# These libraries were trimmed with Cutadapt to remove primers
# before being read into R to start the DADA2 pipeline (following tutorial v 1.16) - see code below
# Silva version 138.1 (Mar 10, 2021 update) used for taxonomic assignments


#### Cutadapt notes (for primer removal) ####

# using template provided by Matt, created excel spreadsheet with sample metadata
# excel spreadsheet = "ROHR_06_Swab_22"
# included primer IDs and sequences
# in separate tab = cutadapt code
# open Terminal (applications > utilities > Terminal), set directory to folder where
# decompressed sequence files are held (in this case = cd ~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_06_decompressed)
# activate cutadapt = conda activate cutadaptenv
# copy cutadapt code from excel and run
# move trimmed sequences to new folder

###downloading new packages###
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("phyloseq")

#BiocManager::install("decontam")

#### Loading required libraries ####

library(dada2)
library(ggplot2)
library(ff)
library(phyloseq)
library(gridExtra)
library(tidyr) #cleaning data & plotting
library(dplyr) # data manipulation
library(decontam) #identify contaminant sequences
library(ShortRead) #for optional read length viz
library(microbiome) #for looking at final phyloseq object

#### Reading in data and visualizing read quality ####

#check that we are in the correct working directory with trimmed fastq files



##Creating filepaths to data
path <- "~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_06_trimmed"
#look at files to check if the path works
list.files(path)


##File preparation
#extracting Forward (fnFs) and Reverse (fnRs) reads from files
fnFs <- sort(list.files(path, pattern = "_R1_001.trimmed.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.trimmed.fastq", full.names = TRUE))

#extract sample names from sequence reads (RRR#s)
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


#plotting quality profiles

#plot quality profile of first 4 forward reads
plotQualityProfile(fnFs[1:4])

#plot quality profile of first 4 reverse reads
plotQualityProfile(fnRs[1:4])

## [optional] ##

#visualize the length of the sequence reads (see that majority of reads 230 nt)

fn <- "~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_06_trimmed/7805_R1_001.trimmed.fastq"
trimmed <- "~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_06_trimmed/7805_R1_001.trimmed_2.fastq"
srq <- readFastq(fn)
seqlen.tab <- table(width(srq)) # Table of how many RSVs have a given length in nts
seqlen.tab
seqlens <- as.integer(names(seqlen.tab))
counts <- as.integer(seqlen.tab)
plot(x=seqlens, y=counts)
plot(x=seqlens, y=log10(counts))

#### Filtering & trimming reads ####


#placing filtered files in a new filtered subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


#filtering and trimming, here truncation at 220 (Fwd) and 180 (Rev) bp, 
#2 expected errors max (N discarded automatically)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220,180),
                     maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

#check how filtering and trimming worked
head(out)

##### Learning error rates ####
errF <- learnErrors(filtFs, multithread = TRUE)
    #100873740 total bases in 458517 reads from 27 samples will be used for learning the error rates.
errR <- learnErrors(filtRs, multithread = TRUE)
    #101115900 total bases in 561755 reads from 35 samples will be used for learning the error rates.

#plotting errors
plotErrors(errF, nominalQ=TRUE)
  
plotErrors(errR, nominalQ=TRUE)


##### Dereplicating reads ####
# (this is not in the current dada2 pipeline (1.16))

sam.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
derepFs <- derepFastq(filtFs)
names(derepFs) <- sam.names
derepRs <- derepFastq(filtRs)
names(derepRs) <- sam.names

#### Inferring Sequence Variants (ASVs) ####
# standard sample-by-sample inference (tested this on 10/19/22)

#dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
#dadaFs[[1]]
#dadaRs <- dada(derepRs, err = errR, multithread = TRUE)
#dadaRs[[1]]

## pooled sample inference - to capture rare ASVs (longest run time)

#dadaFs_pool <- dada(derepFs, err = errF, pool = TRUE, multithread = TRUE)
#dadaFs_pool[[1]]
#dadaRs_pool <- dada(derepRs, err = errR, pool = TRUE, multithread = TRUE)
#dadaRs_pool[[1]]

## pseudo-pooled sample inference - to capture rare ASVs (intermediate run time)
#test this method 1/31/23
##work with these data from here forward

dadaFs_pseudo <- dada(derepFs, err = errF, pool = "pseudo", multithread = TRUE)
dadaFs_pseudo[[1]]
#449 sequence variants were inferred from 1809 input unique sequences.

dadaRs_pseudo <- dada(derepRs, err = errR, pool = "pseudo", multithread = TRUE)
dadaRs_pseudo[[1]]
#395 sequence variants were inferred from 1920 input unique sequences.

#### Merging paired ends (forward & reverse reads) ####
mergers <- mergePairs(dadaFs_pseudo, derepFs, dadaRs_pseudo, derepRs)
#look at merger from first sample
head(mergers[[1]])

ROHR_06 <- makeSequenceTable(mergers)
dim(ROHR_06) #384 33406 (384 samples, 33406 columns - corresponding to unique sequences?)

#look at distribution of sequence lengths
table(nchar(getSequences(ROHR_06)))
    ##wide range of sequence lengths = from 220  to 388 (most clustered around 253)

#exporting files to use in the next part of the workflow
saveRDS(ROHR_06, "~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_06_trimmed/ROHR_06.rds")

#read RDS file (if not re-running code from start)

ROHR_06 <- readRDS("~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_06_trimmed/ROHR_06.rds")

#### Identifying and removing chimeras ####
ROHR_06.nochim <- removeBimeraDenovo(ROHR_06, 
                                    method="pooled", 
                                    multithread=TRUE, verbose = TRUE)
  #Identified 1794 bimeras out of 32776 input sequences (10/19/22) - w/ sequence-by-sequence ASV inference & consensus method
  #Identified 7533 bimeras out of 34203 input sequences w/ pseudo-pooled ASV inference
dim(ROHR_06.nochim)
  #check proportion of merged sequence reads removed
sum(ROHR_06.nochim)/sum(ROHR_06)
  #0.8627225

# tracking changes through each step
getN <- function(x) sum(getUniques(x))
track <-    cbind(out, sapply(dadaFs_pseudo, getN), sapply(dadaRs_pseudo, getN),
                  sapply(mergers, getN), rowSums(ROHR_06.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR",
                     "merged", "nonchim")
rownames(track) <- sam.names
head(track)

write.table(track, "~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_06_trimmed/ROHR_06_read_changes_pseudo.txt", sep = "\t", quote = FALSE,
            col.names=NA)
saveRDS(ROHR_06.nochim, "~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_06_trimmed/ROHR_06.nochim_pseudo.rds")

# Read RDS file of ROHR_06.nochim (if not running code from start)
ROHR_06.nochim <- readRDS("~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_06_trimmed/ROHR_06.nochim_pseudo.rds")

#### Assign taxonomy ####
# Reference datasets formatted for DADA2 can be found here: https://benjjneb.github.io/dada2/training.html
set.seed(119)
ROHR_06.nochim.tax <- assignTaxonomy(ROHR_06.nochim, 
                                    "~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_06_trimmed/silva_nr99_v138.1_train_set.fa.gz",
                                    multithread=TRUE)

# Add species-level assignments based on exact matching of ASVs ##not done (error: vector memory exhausted when run on local computer)
# ROHR_06.nochim.tax <- addSpecies(ROHR_06.nochim.tax, "~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_06_trimmed/silva_species_assignment_v138.1.fa.gz")

# saveRDS(ROHR_06.nochim.tax, "~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_06_trimmed/ROHR_06.nochim.tax.rds")

# Read Phyloseq object
ROHR_06.nochim.tax = readRDS("~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_06_trimmed/ROHR_06.nochim.tax.rds")
ROHR_06.nochim = readRDS("~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_06_trimmed/ROHR_06.nochim.rds")


##############################################################################################

#### Inspect taxonomic assignments ####
taxa.print <- ROHR_06.nochim.tax # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

tax_silva_tax  = tax_table(ROHR_06.nochim.tax)
write.csv(tax_silva_tax, file='~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_06_trimmed/ROHR_06.nochim.tax.csv')

# read in .csv file (if not starting from top)
tax_silva_tax <- read.csv('~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_06_trimmed/ROHR_06.nochim.tax.csv')


tax_silva_otu  = otu_table(ROHR_06.nochim, taxa_are_rows = FALSE)
write.csv(tax_silva_otu, file='~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_06_trimmed/ROHR_06.nochim.tax_otu.csv')

# read in .csv file (if not starting from top)
tax_silva_otu <- read.csv('~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_06_trimmed/ROHR_06.nochim.tax_otu.csv')


#### Creating phyloseq object ####
## Opening and extracting sample data from a .csv file
# creating path to the .csv 
ROHR_06_sam <- "~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_06_trimmed/ROHR_6_metadata.csv"
#reading the sample data sheet
ROHR_06_metadata <- read.csv(ROHR_06_sam, header = TRUE, sep = ",", row.names = 1)

## Creating a unique phyloseq object with sample-by-sequence feature table, 
# the sample metadata and the sequence taxonomies 
ROHR_06_obj <- phyloseq(otu_table(tax_silva_otu, taxa_are_rows = FALSE), 
                       sample_data(ROHR_06_metadata), tax_table(tax_silva_tax)) 

# removing ASVs that are not present in any sample (if any)
ROHR_06_obj <- prune_taxa(taxa_sums(ROHR_06_obj) > 0, ROHR_06_obj)

# Save Phyloseq object
saveRDS(ROHR_06_obj, "~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_06_trimmed/ROHR_06_obj.rds")

# Read Phyloseq object (if not starting from top)
ROHR_06_obj = readRDS("~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_06_trimmed/ROHR_06_obj.rds")


#### Identify & clean contaminants - Prevalence ####

# use control samples ("control_sample" metadata column = control swabs, extraction, and PCR controls) to ID contaminants
sample_data(ROHR_06_obj)$is.neg <- sample_data(ROHR_06_obj)$control_sample == "Y"
contamdf.prev <- isContaminant(ROHR_06_obj, method="prevalence", neg="is.neg")
contaminants <- table(contamdf.prev$contaminant)
#check how many contaminants were ID'd
contaminants #52 true (cont), 25530 false (not cont)

# write .csv with these contaminants
write.csv(contamdf.prev, file='~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_06_trimmed/ROHR_06_contamdf.prev.csv')


# Make phyloseq object of presence-absence in negative controls and true samples
ROHR_06_obj.pa <- transform_sample_counts(ROHR_06_obj, function(abund) 1*(abund>0))
ROHR_06_obj.pa.neg <- prune_samples(sample_data(ROHR_06_obj.pa)$control_sample == "Y", ROHR_06_obj.pa)
ROHR_06_obj.pa.pos <- prune_samples(sample_data(ROHR_06_obj.pa)$control_sample == "N", ROHR_06_obj.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ROHR_06_obj.pa.pos), pa.neg=taxa_sums(ROHR_06_obj.pa.neg),
                    contaminant=contamdf.prev$contaminant)
# Graph to visualize contaminants & target sequences
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")


# Create a new phyloseq object without contaminants identified in contamdf.prev
# (This list was created manually from the table - could be coded to run directly)
badTaxa = c("TACGGAGAGGGCTAGCGTTATTCGGAATTATTGGGCGTAAAGAGCGCGTAGGCTGATTAGTAAGTTAAGAGTGAAATCCCAGAGCTCAACTTTGGAATTGCTTTTAAAACTGCTAATCTAGAGATTGAAAGAGGATAGAGGAATTCCTAGTGTAGAGGTGAAATTCGTAAATATTAGGAGGAACACCAGTGGCGAAGGCGTCTATCTGGTTCAAATCTGACGCTGATGCGCGAAGGCGTGGGGAGCAAACAGG,
            TACAGAGGGTGCGAGCGTTAATCGGATTTACTGGGCGTAAAGCGTGCGTAGGCGGCTTTTTAAGTCGGATGTGAAATCCCCGAGCTTAACTTGGGAATTGCATTCGATACTGGGAAGCTAGAGTATGGGAGAGGATGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGATGGCGAAGGCAGCCATCTGGCCTAATACTGACGCTGAGGTACGAAAGCATGGGGAGCAAACAGG,
            TACAGAGGGTGCGAGCGTTAATCGGATTTACTGGGCGTAAAGCGTGCGTAGGCGGCTAATTGAGTCGGATGTGAAATCCCCGAGCTTAACTTGGGAATTGCATTCGATACTGGTTAGCTAGAGTGTGGGAGAGGATGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGATGGCGAAGGCAGCCATCTGGCCTAACACTGACGCTGAGGTACGAAAGCATGGGGAGCAAACAGG,
            TACGTAGGGGGCGAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGCGGTCCCTTAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGTCATTGGAAACTGGGGGACTTGAGTGCAGAAGAGGAGAGCGGAATTCCACGTGTAGCGGTGAAATGCGTAGAGATGTGGAGGAACACCAGTGGCGAAGGCGGCTCTCTGGTCTGTAACTGACGCTGAGGCGCGAAAGCGTGGGGAGCAAACAGG,
            TACGAAGGGTGCAAGCGTTACTCGGAATTACTGGGCGTAAAGCGTGCGTAGGTGGTTTGTTAAGTCTGCTGTGAAAGCCCTGGGCTCAACCTGGGAATTGCAGTGGATACTGGCTGACTAGAATGTGGCAGAGGGTAGCGGAATTCCTGGTGTAGCAGTGAAATGCGTAGAGATCAGGAGGAACATCCGTGGCGAAGGCGGCTACCTGGGCCAACATTGACACTGAGGCACGAAAGCGTGGGGAGCAAACAGG,
            TACGTAGGGTCCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGTTGTGCAAGACCGATGTGAAATCCCCGGGCTTAACCTGGGAATTGCATTGGTGACTGCACGGCTAGAGTGTGTCAGAGGGGGGTAGAATTCCACGTGTAGCAGTGAAATGCGTAGAGATGTGGAGGAATACCGATGGCGAAGGCAGCCCCCTGGGATAACACTGACGCTCATGCACGAAAGCGTGGGGAGCAAACAGG,
            TACCGAGGGTGCAAGCTTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATCCAAAACTGGCAAGCTAGAGTACGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAAACAGG,
            GACAGAGGATGCCAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTCAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTACCAAGCTGGAGTACGGTAGGGGCAGAGGGAATTTCCGGTGGAGCGGTGAAATGCGTAGAGATCGGAAAGAACACCAACGGCGAAAGCACTCTGCTGGGCCGACACTGACACTGAGAGACGAAAGCTAGGGGAGCAAATGGG,
            TACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTGATTTAAGTCTGATGTGAAAGCCCCCAGCTCAACTGGGGAGGGTCATTGGAAACTGGATCACTTGAGTGCAGAAGAGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCAGTGGCGAAGGCGGCTCTCTGGTCTGTAACTGACGCTGAGGCTCGAAAGCGTGGGTAGCAAACAGG,
            TACGAAGGGGGCTAGCGTTGTTCGGATTTACTGGGCGTAAAGCGCACGTAGGCGGACTTTTAAGTCAGGGGTGAAATCCCGGGGCTCAACCCCGGAACTGCCTTTGATACTGGAAGTCTTGAGTATGGTAGAGGTGAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAGGAACACCAGTGGCGAAGGCGGCTCACTGGACCATTACTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGG,
            TACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGTTTGCTAAGACCGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTGGTGACTGGCAGGCTAGAGTATGGCAGAGGGGGGTAGAATTCCACGTGTAGCAGTGAAATGCGTAGAGATGTGGAGGAATACCGATGGCGAAGGCAGCCCCCTGGGCCAATACTGACGCTCATGCACGAAAGCGTGGGGAGCAAACAGG,
            TACGTAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGTTATGCAAGACAGAGGTGAAATCCCCGGGCTCAACCTGGGAACTGCCTTTGTGACTGCATGGCTAGAGTACGGTAGAGGGGGATGGAATTCCGCGTGTAGCAGTGAAATGCGTAGATATGCGGAGGAACACCGATGGCGAAGGCAATCCCCTGGACCTGTACTGACGCTCATGCACGAAAGCGTGGGGAGCAAACAGG,
            TACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGCGGTCCTTTAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGTCATTGGAAACTGGGGGACTTGAGTGCAGAAGAGGAGAGCGGAATTCCACGTGTAGCGGTGAAATGCGTAGAGATGTGGAGGAACACCAGTGGCGAAGGCGGCTCTCTGGTCTGTAACTGACGCTGAGGCGCGAAAGCGTGGGGAGCAAACAGG,
            GTAAGACTTAAGGATTTATTTTTTAATAAAAAGCAAAAAGCGTGTTAAGGATTTTTTAAAAAAAAAAATAAATAGAATTTTTTTCGTAATTGTAATATGTTAAAATGAAAAAAAGAATTTTTTATATGAAGATAATTTATTTTTTTTTTCTTAAATACGAAGGTTTGGGGAGCAAATAGGATTAGAAACCCCTGTAGTCCAGTAGATCGGAAGAGCACACGTCT,
            TACGAAGGGGGCTAGCTTTGTTCGGATTTACTGGGCGTAAAGCGCACGTAGGCGGACTATTAAGTCAGGGGTGAAATCCCGGGGCTCAACCCCGGAACTGCCTTTGATACTGGTAGTCTTGAGTTCGAGAGAGGTGAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAGGAACACCAGTGGCGAAGGCGGCTCACTGGCTCGATACTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGG,
            TACGTAGGGGGCGAGCGTTGTCCGGAATGATTGGGCGTAAAGCGCGCGCAGGCGGTCCTTTAAGTCTGATGTAAAAGCCCGCGGCTTAACCGCGGAAGGTCATTGGAAACTGGGGGACTTGAGGCTAGGAGAGGGAAGTGGAATTCCTGGTGTAGCGGTGAAATGCGTAGAGATCAGGAGGAATACCGATGGCGAAAGCAACTTCCTGGCTTAGAACTGACGCTGAGGCGCGAAAGCGTGGGGAGCAAACAGG,
            GTAAGACTTAAGGATTTATTTTTTAATAAAAAGCAAAAAGCGTGTTAAGGATTTTTTAAAAAAAAAAATAAATAGAATTTTTTTCGTAATTGTAATATGTTAAAATGAAAAAAAGAATTTTTTATATGAAGATAATTTATTTTTTTTTTCTTAAATACGAAGGTTTGGGGAGCAAATAGGATTAGATACCCTTGTAGTCCAGTAGATCGGAAGAGCACACGTCT,
            TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGG,
            TACGGAGGGTGCAAGCGTTGTCCGGATTTATTGGGTTTAAAGGGTGCGTAGGTGGTTTATTAAGTCAGTGGTGAAAGGCTACAGCTCAACTGTAGGAGTGCCATTGATACTGGTAGACTTGAGTTTGGTCGAGGTAGGCGGAATTTATGGTGTAGCGGTGAAATGCATAGATACCATAAAGAACACCGATAGTGAAGACAGCTTACTAGGCCTGAACTGACACTGAGGCACGAAAGCGTGGGGAGCGAACAGG,
            TACGTAAAAGACTAGTGTTAGTCATCTTTATTAGGTTTAAAGGGTACCTAGACGGTAAATTAAACTCTAAATGAGTACTTTTTTACTAGAGTTTTATAAGAGAAGGAAGAATTTCTGGAGTAGTGATAGAATATTTTAATACCAGAAGGACTGATAACGGCGAAGGCGTCCTTCTATGTAAAAACTGACGTTGAGGGACGAAGGCTTGGGTAGCAATAAGG,
            TACGGAGGGGGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGATCGGAAAGTTGGGGGTGAAATCCCGGGGCTCAACCTCGGAACTGCCTTCAAAACTACTGGTCTGGAGTTCGAGAGAGGTGAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAGGAACACCAGTGGCGAAGGCGGCTCACTGGCTCGATACTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGG,
            GACGTAGGGGGCGAGCGTTGTCCGGAGTTACTGGGCGTAAAGGGCGCGTAGGCGGTTTAGCAAGTCAGGTGTAAAAGGCCACGGCTCAACCGTGGATATGCATCTGAAACTGCTGAGCTAGAGGGCAGGAGAGGGGAGTGGAATTCCCGGTGTAGCGGTGAAATGCGTAGATATCGGGAGGAATACCAGTGGCGAAGGCGACTCCCTGGACTGGCCCTGACGCTGAGGCGCGAGAGCGTGGGGAGCAAACAGG,
            TACGGAGGGTGCAAGCGTTGTCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGATCTTTAAGTCAGTGGTGAAAGCCTGCAGCTTAACTGTAGAACTGCCATTGATACTGGAGGTCTTGAGTACACTAGAGGTAGGCGGAATTTATGGTGTAGCGGTGAAATGCATAGATACCATAAAGAACACCTATAGCGTAGGCAGCTTACTGGAGTGTAACTGACGCTGAGGCACGAAAGCATGGGGAGCGAACAGG,
            TACGTAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGTTGTGCAAGACAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTTGTGACTGCACGGCTGGAGTGCGGCAGAGGGGGATGGAATTCCGCGTGTAGCAGTGAAATGCGTAGATATGCGGAGGAACACCGATGGCGAAGGCAATCCCCTGGGCCTGCACTGACGCTCATGAACGAAAGCGTGGGGAGCAAACAGG,
            TACGGAGGGAGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGCTTTGTAAGTTAGAGGTGAAAGCCCGGGGCTCAACCCCGGAACTGCCTTTAAGACTGCATCGCTTGAATCCAGGAGAGGTGAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAAGAACACCAGTGGCGAAGGCGGCTCACTGGACTGGTATTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGG,
            TACGAAGGGGGCTAGCTTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGATATTTAAGTCAGGGGTGAAATCCCAGAGCTCAACTCTGGAACTGCCTTTGATACTGGGTATCTTGAGTATGGAAGAGGTAAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAGGAACACCAGTGGCGAAGGCGGCTTACTGGTCCATTACTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGG,
            TACGAAGGGTGCAAGCGTTACTCGGAATTACTGGGCGTAAAGCGTGCGTAGGTGGTGGCTTAAGTCCGTTGTGAAAGCCCTGGGCTCAACCTGGGAATTGCAGTGGATACTGGGTCACTAGAGTGTGGTAGAGGGTGGCGGAATTCCCGGTGTAGCAGTGAAATGCGTAGAGATCGGGAGGAACACCCGTGGCGAAGGCGGCCACCTGGGCCAACACTGACACTGAGGCACGAAAGCGTGGGGAGCAAACAGG,
            TACGAAGGGGGCTAGCGTTGCTCGGAATCACTGGGCGTAAAGGGCGCGTAGGCGGCTATTTAAGTCGGGGGTGAAAGCCTGTGGCTCAACCACAGAATTGCCTTCGATACTGGATAGCTTGAGTATGGTAGAGGTTGGTGGAACTGCGAGTGTAGAGGTGAAATTCGTAGATATTCGCAAGAACACCGGTGGCGAAGGCGGCCAACTGGACCATTACTGACGCTGAGGCGCGAAAGCGTGGGGAGCAAACAGG,
            TACGGAGGGCGCGAGCTTTACCCGGATTCACTGGGCGTAAAGGGCGTGTAGGCGGCCTGGGGCGTCCCATGTGAAAGACCACGGCTCAACCGTGGGGGAGCGTGGGATACGCTCAGGCTAGACGGTGGGAGAGGGTGGTGGAATTCCCGGAGTAGCGGTGAAATGCGCAGATACCGGGAGGAACGCCGATGGCGAAGGCAGCCACCTGGTCCACCCGTGACGCTGAGGAGCGAAAGCGTGGGGAGCAAACCGG,
            TACGGAGGGTGCGAGCGTTGTCCGGATTTATTGGGTTTAAAGGGTGCGTAGGTGGCTGATTAAGTCAGTGGTGAAAGCCCATCGCTTAACGATGGAAATGCCATTGATACTGGTTAGCTTGAGTATCGTTGAAGTAGGCGGAATTTATGGTGTAGCGGTGAAATGCATAGATACCATAAAGAACACCGATAGCGTAGGCAGCTTACTAAGCGATAACTGACACTGAGGCACGAAAGCATGGGGAGCGAACAGG,
            TACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGTGATTTAAGTCTGATGTGAAAGCCCCCAGCTCAACTGGGGAGGGTCATTGGAAACTGGATCACTTGAGTGCAGAAGAGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCAGTGGCGAAGGCGGCTCTCTGGTCTGTAACTGACGCTGAGGCTCGAAAGCGTGGGTAGCAAACAGG,
            TACGGAGGGTGCAAGCGTTACCCGGAATCACTGGGCGTAAAGGGCGTGTAGGCGGTTATTTAAGTCTGGTTTTAAAGACCGGGGCTCAACCCCGGGAGTGGACTGGATACTGGATGACTTGACCTCTGGAGAGGGAACTGGAATTCCTGGTGTAGCGGTGGAATGCGTAGATACCAGGAGGAACACCAATGGCGAAGGCAAGTTCCTGGACAGAAGGTGACGCTGAGGCGCGAAAGTGTGGGGAGCGAACCGG,
            TACGAAGGGTGCAAGCTTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTCAGCAAGTTGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCCAAAACTACTGAGCTAGAGTACGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAAACAGG,
            TACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGTTGTGTAAGACAGGTGTGAAATCCCCGGGCTTAACCTGGGAACTGCGCTTGTGACTGCACGGCTAGAGTATGGCAGAGGGGGGTGGAATTCCACGTGTAGCAGTGAAATGCGTAGAGATGTGGAGGAACACCGATGGCGAAGGCAGCCCCCTGGGCCAATACTGACGCTCATGAACGAAAGCGTGGGGAGCAAACAGG,
            TACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTTCTTAAGTCTGATGTGAAAGCCCCCGGCTCAACCGGGGAGGGTCATTGGAAACTGGGAGACTTGAGTGCAGAAGAGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCAGTGGCGAAGGCGGCTCTCTGGTCTGCAACTGACGCTGAGGATCGAAAGCGTGGGGAGCAAACAGG,
            TACGTAGGGGGCAAGCGTTATCCGGATTCATTGGGCGTAAAGCGCTCGTAGGCGGTCTGTTAGGTCGGGAGTTAAATCCGGGGGCTCAACCCCCGTTCGCTCCCGATACCGGCAGACTTGAGTTTGGTAGGGGAAGGTGGAATTCCTAGTGTAGCGGTGGAATGCGCAGATATTAGGAAGAACACCAGTGGCGAAGGCGGCCTTCTGGGCCATAACTGACGCTGAGGAGCGAAAGCTAGGGGAGCAAACAGG,
            TACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGTTTTTTAAGTCAGATGTGAAAGCCCCGGGCTCAACCTGGGAATTGCATTTGAAACTGGGAAACTAGAGTGTGTGAGAGGGGGGTAGAATTCCAAGTGTAGCGGTGAAATGCGTAGAGATTTGGAGGAATACCAGTGGCGAAGGCGGCCCCCTGGCACAACACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGG,
            TACGGAGGGGGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGCTTTGTAAGTTAGAGGTGAAAGCCTGGAGCTCAACTCCAGAACTGCCTTTAAGACTGCATCGCTTGAATCCAGGAGAGGTGAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAAGAACACCAGTGGCGAAGGCGGCTCACTGGACTGGTATTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGG,
            TACGGAGGATGCAAGCGTTATTCGGAATTATTGGGCGTAAAGCGTCTGTAGGTGGTTTTTTAAGTCTACTGTTAAATATTAAGGCTTAACCTTAAAAAAGCGGTATGAAACTAAAAAACTTGAGTTTAGTAGAGGTAGAGGGAATTCTCGGTGTAGTGGTGAAATGCGTAGAGATCGAGAAGAACACCGGTAGCGAAAGCGCTCTACTGGGCTAAAACTGACACTGAGAGACGAAAGCTAGGGGAGCAAATAGG,
            TACGGAGGGGGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGACCGGAAAGTTGGGGGTGAAATCCCGGGGCTCAACCTCGGAACTGCCTTCAAAACTATCGGTCTTGAGTTCGAGAGAGGTGAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAGGAACACCAGTGGCGAAGGCGGCTCACTGGCTCGATACTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGG,
            TACGTAGGTGACAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTCTTTTAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGTCATTGGAAACTGGAAGACTTGAATGCAGAAGAGGAAAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGTAATTGACGCTGAGGATCGAAAGCGTGGGGAGCGAACAGG,
            TACGTAGGGTGCGAGCGTTGTCCGGAATTACTGGGCGTAAAGAGCTCGTAGGTGGTTTGTCGCGTCGTCTGTGAAATTCCGGGGCTTAACTCCGGGCGTGCAGGCGATACGGGCATAACTTGAGTGCTGTAGGGGAGACTGGAATTCCTGGTGTAGCGGTGAAATGCGCAGATATCAGGAGGAACACCGATGGCGAAGGCAGGTCTCTGGGCAGTTACTGACGCTGAGGAGCGAAAGCATGGGTAGCGAACAGG,
            TACGTAGGGTGCGAGCGTTGTCCGGAATTATTGGGCGTAAAGAGCTTGTAGGCGGTTTGTCGCGTCTGCTGTGAAAATCCGGGGCTCAACCCCGGACTTGCAGTGGGTACGGGCAGACTAGAGTGTGGTAGGGGAGACTGGAATTCCTGGTGTAGCGGTGAAATGCGCAGATATCAGGAGGAACACCGATGGCGAAGGCAGGTCTCTGGGCCACTACTGACGCTGAGAAGCGAAAGCATGGGGAGCGAACAGG,
            TACAGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGCGGTTTGTTAAGTCGGATGTGAAATCCCCGGGCTCAACCTGGGAATTGCATTCGATACTGGCAAGCTAGAGTATGGGAGAGGAAGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGATGGCGAAGGCAGCCTTCTGGCCTAATACTGACGCTGAGGTGCGAAAGCATGGGGAGCAAACAGG,
            TACGGAGGGTGCAAGCATTAATCGGATTTACTGGGTGTAAAGAGCGCGTAGGCGGGTAGGCAAGTCAGATGTGAAATTCCGGAGCTCAACTCCGGAGCTGCATTTGAAACTACTTATCTTGAGGGTGGACGGAGAAAACGGAATTCCACGTGTAGCGGTGAAATGCGTAGATATGTGGAAGAACACCGGTGGCGAAGGCGGTTTTCTAGTTCATTCCTGACGCTGAGGCGCGAGAGCAAGGGGAGCAAACAGG,
            TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTCTGTCAAGTCGGATGTGAAATCCCCGGGCTCAACCCGGGAACTGCATTCGAAACTGGCAGACTAGAGTCTTGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACAAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGG,
            TACGAAGGGGGCTAGCGTTGCTCGGAATCACTGGGCGTAAAGGGCGCGTAGGCGGCCGATTAAGTCGGGGGTGAAAGCCTGTGGCTCAACCACAGAATTGCCTTCGATACTGGTTGGCTTGAGACCGGAAGAGGACAGCGGAACTGCGAGTGTAGAGGTGAAATTCGTAGATATTCGCAAGAACACCAGTGGCGAAGGCGGCTGTCTGGTCCGGTTCTGACGCTGAGGCGCGAAAGCGTGGGGAGCAAACAGG,
            TACGCAGGTGGCGAGCGTTGCCCGGAATTACTGGGCGTAAAGGGTGCGTAGGCGGCCGGACAAGTCATAGGTTAAAGCCCGGAGCTCAACTCCGGAAAGGCCTATGATACTGTCTGGCTTGAGGGCCGGAGAGGCTGGCGGAATTCCCGGTGTAGGGGTGAAATCCGTAGATATCGGGAGGAACACCGGTGGGGAAGCCGGCCAGCTGGACGGTTCCTGACGCTGAGGAACGAAAGCGTGGGGAGCAAACCGG,
            TACGTAGGGGGCGAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGTAGGCGGCCGGTTAAGTCCTGTGTGAAAGACCACGGCTCAACCGTGGGGGTGCACGGGAAACTGGCCGGCTTGAGTGCAGGAGAGGGGAGCGGAATTCCCGGTGTAGCGGTGGAATGCGTAGAGATCGGGAGGAACACCAGTGGCGAAGGCGGCTCTCTGGCCTGTAACTGACGCTGAGGCGCGAAAGCGTGGGGAGCGAACAGG,
            TACGTAGGGTGCGAGCGTTGTCCGGAATTATTGGGCGTAAAGGGCTCGTAGGCGGTCTGTCACGTCGGGAGTGAAAACTCGGGGCTCAACCCCGAGCCTGCTTCCGATACGGGCAGACTAGAGGTATGCAGGGGAGAACGGAATTCCTGGTGTAGCGGTGAAATGCGCAGATATCAGGAGGAACACCGGTGGCGAAGGCGGTTCTCTGGGCATTACCTGACGCTGAGGAGCGAAAGTGTGGGGAGCGAACAGG,
            TACGGAGGGGGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGACGGAGAAGTCAGAGGTGAAATCCCAGGGCTCAACCTTGGAACTGCCTTTGAAACTTTCTGTCTTGAGGTCGAGAGAGGTGAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAGGAACACCAGTGGCGAAGGCGGCTCACTGGCTCGATACTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGG,
            TACGAAGGGGGCTAGCGTTGCTCGGAATCACTGGGCGTAAAGGGTGCGTAGGCGGGTCTTTAAGTCAGGGGTGAAATCCTGGAGCTCAACTCCAGAACTGCCTTTGATACTGAGGATCTTGAGTCCGGAAGAGGTGAGTGGAACTGCGAGTGTAGAGGTGAAATTCGTAGATATTCGCAAGAACACCAGTGGCGAAGGCGGCTCACTGGTCCGGTACTGACGCTGAGGCACGAAAGCGTGGGGAGCAAACAGG")

goodTaxa <- setdiff(taxa_names(ROHR_06_obj), badTaxa)
ROHR_06_obj <- prune_taxa(goodTaxa, ROHR_06_obj)

#### Remove all control samples (swab, EXT, PCR) ####

ROHR_06_obj = subset_samples(ROHR_06_obj, sample_names(ROHR_06_obj) != "8733-CC-NU")
ROHR_06_obj = subset_samples(ROHR_06_obj, sample_names(ROHR_06_obj) != "15106-CC-NU")
ROHR_06_obj = subset_samples(ROHR_06_obj, sample_names(ROHR_06_obj) != "8740-CC-NU")
ROHR_06_obj = subset_samples(ROHR_06_obj, sample_names(ROHR_06_obj) != "15107-CC-NU")
ROHR_06_obj = subset_samples(ROHR_06_obj, sample_names(ROHR_06_obj) != "15110-CC-NU")
ROHR_06_obj = subset_samples(ROHR_06_obj, sample_names(ROHR_06_obj) != "8736-CC-NU")
ROHR_06_obj = subset_samples(ROHR_06_obj, sample_names(ROHR_06_obj) != "8743-CC-NU")
ROHR_06_obj = subset_samples(ROHR_06_obj, sample_names(ROHR_06_obj) != "8734-CC-NU")
ROHR_06_obj = subset_samples(ROHR_06_obj, sample_names(ROHR_06_obj) != "8744-CC-NU")
ROHR_06_obj = subset_samples(ROHR_06_obj, sample_names(ROHR_06_obj) != "CLP-NU-1")
ROHR_06_obj = subset_samples(ROHR_06_obj, sample_names(ROHR_06_obj) != "CC-UP-1")
ROHR_06_obj = subset_samples(ROHR_06_obj, sample_names(ROHR_06_obj) != "CLP-UP-1")
ROHR_06_obj = subset_samples(ROHR_06_obj, sample_names(ROHR_06_obj) != "CLP-NU-2")
ROHR_06_obj = subset_samples(ROHR_06_obj, sample_names(ROHR_06_obj) != "CLP-UP-2")
ROHR_06_obj = subset_samples(ROHR_06_obj, sample_names(ROHR_06_obj) != "CLP-NU-3")
ROHR_06_obj = subset_samples(ROHR_06_obj, sample_names(ROHR_06_obj) != "CLP-UP-3")
ROHR_06_obj = subset_samples(ROHR_06_obj, sample_names(ROHR_06_obj) != "CLP-UP-4")
ROHR_06_obj = subset_samples(ROHR_06_obj, sample_names(ROHR_06_obj) != "CLP-UP-5")
ROHR_06_obj = subset_samples(ROHR_06_obj, sample_names(ROHR_06_obj) != "CLP-UP-6")
ROHR_06_obj = subset_samples(ROHR_06_obj, sample_names(ROHR_06_obj) != "pcr-control-1")
ROHR_06_obj = subset_samples(ROHR_06_obj, sample_names(ROHR_06_obj) != "pcr-control-2")
ROHR_06_obj = subset_samples(ROHR_06_obj, sample_names(ROHR_06_obj) != "pcr-control-3")
ROHR_06_obj = subset_samples(ROHR_06_obj, sample_names(ROHR_06_obj) != "ext-control-1")
ROHR_06_obj = subset_samples(ROHR_06_obj, sample_names(ROHR_06_obj) != "pcr-control-4")
ROHR_06_obj = subset_samples(ROHR_06_obj, sample_names(ROHR_06_obj) != "ext-control-2")

#check that all controls were removed (should not have any CC/CLP/pcr-control/ext-control)
sample_names(ROHR_06_obj)

#### Look at distribution of library size ####
df <- as.data.frame(sample_data(ROHR_06_obj)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ROHR_06_obj)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color = control_sample)) + geom_point()

#### Remove outliers, chloroplasts, mitochondria, eukaryotes... ####
## Remove outlier samples
## (In this case, did not remove any at this step)
ROHR_06_obj = subset_samples(ROHR_06_obj, sample_names(ROHR_06_obj) != "")


## Removing ASVs of chloroplasts
ROHR_06_obj <- subset_taxa(ROHR_06_obj, Order != "Chloroplast" | is.na(Order))
ROHR_06_obj

## Removing ASVs of mitochondria
ROHR_06_obj <- subset_taxa(ROHR_06_obj, Family != "Mitochondria" | is.na(Family))
ROHR_06_obj

## Removing ASVs assigned to Eukaryotes and unassigned (NA) at the kingdom level
ROHR_06_obj <- subset_taxa(ROHR_06_obj, Kingdom != "Eukaryota")
ROHR_06_obj

# Remove empty rows
ROHR_06_obj <- prune_taxa(taxa_sums(ROHR_06_obj) > 0, ROHR_06_obj)
ROHR_06_obj

# Check what taxonomic classifications are in the phyloseq object (optional)
rank_names(ROHR_06_obj)
  # have "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"  

# List out sample variables (taken from metadata file - includes RRR, ocean, region, site...)
sample_variables(ROHR_06_obj)

#### Save Phyloseq object ####
saveRDS(ROHR_06_obj, "~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_06_trimmed/ROHR_06_obj_clean.rds")

##############################################################################################

#### Read Phyloseq object ####
ROHR_06_obj_clean = readRDS("/Volumes/LLL_one/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_06_trimmed/ROHR_06_obj_clean.rds")

##############################################################################################

# Summary of this original phyloseq object (microbiome package)
summarize_phyloseq(ROHR_06_obj_clean)

# min number of reads = 1
# max number of reads = 48663
# total number of reads = 4418389
# avg. number of reads = 12307.5
# median number of reads = 7967
# Sparsity = 0.987635948358637
# Any OTU sum to 1 or less? YES
# Number of singletons = 29
# Percent of OTUs that are singletons \n  (i.e. exactly one read detected across all samples) = 0.123720136518771
# number of sample variables = 43

# number of taxa in OTU table
ntaxa(ROHR_06_obj_clean) #23440
nsamples(ROHR_06_obj_clean) #359

# Check how many samples of each fish species were collected in each gulf & season

## Extract the sample metadata
sample_data <- as.data.frame(sample_data(ROHR_06_obj_clean))

# Create a contingency table of species by region and season
species_region_season_counts <- table(sample_data$species, sample_data$region, sample_data$season)

# Convert the table to a data frame for easier manipulation
species_region_season_df <- as.data.frame(species_region_season_counts)

# Rename columns for clarity
colnames(species_region_season_df) <- c("Species", "Region", "Season", "Count")

# Pivot data to have `Season` as columns and add row totals for each Species and Region combination
species_region_season_wide <- species_region_season_df %>%
  pivot_wider(names_from = Season, values_from = Count, values_fill = 0) %>%
  mutate(`n by Region` = rowSums(select(., -Species, -Region)))

# Calculate the summed total for each species across all regions
species_totals <- species_region_season_wide %>%
  group_by(Species) %>%
  summarize(`Total by Species` = sum(`n by Region`))

# Add the species totals to the main table
species_region_season_wide <- species_region_season_wide %>%
  left_join(species_totals, by = "Species")

# Mark "Total by Species" as `NA` except in the last region row for each species
species_region_season_wide <- species_region_season_wide %>%
  group_by(Species) %>%
  mutate(`Total by Species` = ifelse(row_number() == n(), `Total by Species`, NA)) %>%
  ungroup()

# Calculate column totals without duplicating the "n by Region" column in the bottom row
column_totals <- colSums(select(species_region_season_wide, -Species, -Region), na.rm = TRUE)
column_totals <- c(Species = "Total", Region = "", as.list(column_totals))
column_totals["n by Region"] <- NA  # Remove the duplicate in the "n by Region" column for the bottom row

# Bind column totals as the last row
species_region_season_table <- rbind(species_region_season_wide, column_totals)

# Arrange rows to group by Species, and display the final table
species_region_season_table <- species_region_season_table %>%
  arrange(Species, Region)

print(species_region_season_table)

# Save table as csv file

# Specify correct file path
file_path <- "/Volumes/LLL_one/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_06_trimmed/species_region_season_table.csv"

# Save the table as a CSV file
write.csv(species_region_season_table, file = file_path, row.names = FALSE)

# Message to confirm saving & correct location
cat("Table saved as:", file_path)