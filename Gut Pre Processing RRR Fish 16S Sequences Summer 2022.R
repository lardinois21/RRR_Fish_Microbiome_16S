#### Draft code for reading and cleaning 16S sequences ####
## Base code adapted from dada2 pipeline tutorial by Matthieu Leray, edited by Laura Lardinois
# November 2022 -- McGill & STRI

#### Current project: RRR Fish Microbiome - TEP Gut Samples ####

## Project summary ##
# Fish were collected from 2 sites in the Tropical Eastern Pacific: Las Perlas & Coiba
# Sampling during both upwelling & non-upwelling seasons (collections done 2021-22)
# 10 target species, 8 of which are represented in this dataset
# Swabs = 9 per species
# Gut = 8 per species(most matched to skin)
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
# These libraries were sent to me (Laura Lardinois) and trimmed with Cutadapt to remove primers
# before being read into R to start the DADA2 pipeline (following tutorial v 1.16) - see code below
# Silva version 138.1 (mar 10, 2021 update) used for taxonomic assignments


#### Cutadapt notes (for primer removal) ####

# using template provided by Matt, created excel spreadsheet with sample metadata
# make sure that sample file names match the .fasq files specified in the cutadapt command line
# excel spreadsheet = "ROHR_06_Swab_22" --> includes both swab and gut files for now
# included primer IDs and sequences
# in separate tab = cutadapt code
# open Terminal (applications > utilities > Terminal), set directory to folder where
# decompressed sequence files are held (in this case = cd /Users/lllardinois/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_07_decompressed)
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
library(dplyr)
library(decontam) #identify contaminant sequences
library(ShortRead) #for optional read length viz

#### Reading in data and visualizing read quality ####

#check that we are in the correct working directory with trimmed fastq files



##Creating filepaths to data
path <- "~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_07_trimmed"
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

#### [optional] ####

#visualize the length of the sequence reads (see that majority of reads 230 nt)

fn <- "~/Documents/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_06_trimmed/7805_R1_001.trimmed.fastq"
trimmed <- "~/Documents/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_06_trimmed/7805_R1_001.trimmed_2.fastq"
srq <- readFastq(fn)
seqlen.tab <- table(width(srq)) # Table of how many RSVs have a given length in nts
seqlen.tab
seqlens <- as.integer(names(seqlen.tab))
counts <- as.integer(seqlen.tab)
plot(x=seqlens, y=counts)
plot(x=seqlens, y=log10(counts))

#### filtering & trimming ####


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

##### learning error rates ####
errF <- learnErrors(filtFs, multithread = TRUE)
    #101575760 total bases in 461708 reads from 23 samples will be used for learning the error rates.
errR <- learnErrors(filtRs, multithread = TRUE)
    #101007000 total bases in 561150 reads from 28 samples will be used for learning the error rates.

#plotting errors
plotErrors(errF, nominalQ=TRUE)
  
plotErrors(errR, nominalQ=TRUE)


##### Dereplicating reads (this is not in the current dada2 pipeline (1.16)) ####
sam.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
derepFs <- derepFastq(filtFs)
names(derepFs) <- sam.names
derepRs <- derepFastq(filtRs)
names(derepRs) <- sam.names

#### Inferring Sequence Variants ####
# standard sample-by-sample inference (tested this on 11/22/22)

dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaFs[[1]]
#dada-class: object describing DADA2 denoising results
#64 sequence variants were inferred from 7268 input unique sequences.
#Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

dadaRs <- dada(derepRs, err = errR, multithread = TRUE)
dadaRs[[1]]
#dada-class: object describing DADA2 denoising results
#49 sequence variants were inferred from 4849 input unique sequences.
#Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

## pooled sample inference - to capture rare ASVs (longest run time)

dadaFs_pool <- dada(derepFs, err = errF, pool = TRUE, multithread = TRUE)
dadaFs_pool[[1]]
dadaRs_pool <- dada(derepRs, err = errR, pool = TRUE, multithread = TRUE)
dadaRs_pool[[1]]

## pseudo-pooled sample inference - to capture rare ASVs (intermediate run time)
#used this method for swab & gut samples (11/22/22) - started running at 11:08 am - 
##work with these data from here forward

dadaFs_pseudo <- dada(derepFs, err = errF, pool = "pseudo", multithread = TRUE)
dadaFs_pseudo[[1]]
#dada-class: object describing DADA2 denoising results
#72 sequence variants were inferred from 7268 input unique sequences.
#Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

dadaRs_pseudo <- dada(derepRs, err = errR, pool = "pseudo", multithread = TRUE)
dadaRs_pseudo[[1]]
#dada-class: object describing DADA2 denoising results
#60 sequence variants were inferred from 4849 input unique sequences.
#Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16


##Merging paired ends (forward & reverse reads)
mergers <- mergePairs(dadaFs_pseudo, derepFs, dadaRs_pseudo, derepRs)
#look at merger from first sample
head(mergers[[1]])

ROHR_07 <- makeSequenceTable(mergers)
dim(ROHR_07) # 384 21485 (384 samples, 21485 columns - corresponding to unique sequences?)

#look at distribution of sequence lengths
table(nchar(getSequences(ROHR_07)))
    ##wide range of sequence lengths = from  220 to 387 (most clustered around 253)
#remove non-target sequences??

#exporting files to use in the next part of the workflow
saveRDS(ROHR_07, "~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_07_trimmed/ROHR_07.rds")

#read RDS file (if not re-running code from start)

ROHR_07 <- readRDS("~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_07_trimmed/ROHR_07.rds")

#### Identifying and removing chimeras ####
ROHR_07.nochim <- removeBimeraDenovo(ROHR_07, 
                                    method="pooled", 
                                    multithread=TRUE, verbose = TRUE)
  #Identified 6559 bimeras out of 21485 input sequences.

dim(ROHR_07.nochim) #384, 14926
  #check proportion of merged sequence reads removed
sum(ROHR_07.nochim)/sum(ROHR_07)
  # 0.8995282 --> kept ~90% of reads

#tracking changes through each step
getN <- function(x) sum(getUniques(x))
track <-    cbind(out, sapply(dadaFs_pseudo, getN), sapply(dadaRs_pseudo, getN),
                  sapply(mergers, getN), rowSums(ROHR_07.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR",
                     "merged", "nonchim")
rownames(track) <- sam.names
head(track)

write.table(track, "~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_07_trimmed/ROHR_07_read_changes_pseudo.txt", sep = "\t", quote = FALSE,
            col.names=NA)
saveRDS(ROHR_07.nochim, "~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_07_trimmed/ROHR_07.nochim_pseudo.rds")

#read RDS file of ROHR_07.nochim (if not running code from start)
ROHR_07.nochim <- readRDS("~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_07_trimmed/ROHR_07.nochim_pseudo.rds")

#### Assign taxonomy ####
#reference datasets formatted for DADA2 can be found here: https://benjjneb.github.io/dada2/training.html
set.seed(119)
ROHR_07.nochim.tax <- assignTaxonomy(ROHR_07.nochim, 
                                    "~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_07_trimmed/silva_nr99_v138.1_train_set.fa.gz",
                                    multithread=TRUE)

#add species-level assignments based on exact matching of ASVs ##not done yet (error: vector memory exhausted when run on local computer)
ROHR_07.nochim.tax <- addSpecies(ROHR_07.nochim.tax, "~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_07_trimmed/silva_species_assignment_v138.1.fa.gz")

saveRDS(ROHR_07.nochim.tax, "~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_07_trimmed/ROHR_07.nochim.tax.rds")

# Read Phyloseq object
ROHR_07.nochim.tax = readRDS("~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_07_trimmed/ROHR_07.nochim.tax.rds")
ROHR_07.nochim = readRDS("~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_07_trimmed/ROHR_07.nochim.rds")


##############################################################################################

#### Inspect taxonomic assignments ####
taxa.print <- ROHR_07.nochim.tax # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

tax_silva_tax  = tax_table(ROHR_07.nochim.tax)
write.csv(tax_silva_tax, file='~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_07_trimmed/ROHR_07.nochim.tax.csv')

#read in .csv file (if not starting from top)
tax_silva_tax <- read.csv('~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_07_trimmed/ROHR_07.nochim.tax.csv')


tax_silva_otu  = otu_table(ROHR_07.nochim, taxa_are_rows = FALSE)
write.csv(tax_silva_otu, file='~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_07_trimmed/ROHR_07.nochim.tax_otu.csv')

#read in .csv file (if not starting from top)
tax_silva_otu <- read.csv('~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_07_trimmed/ROHR_07.nochim.tax_otu.csv')


##Creating phyloseq object
##Opening and extracting sample data from a .csv file
#creating path to the .csv 
ROHR_07_sam <- "~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_07_trimmed/ROHR_7_metadata.csv"
#reading the sample data sheet
ROHR_07_metadata <- read.csv(ROHR_07_sam, header = TRUE, sep = ",", row.names = 1)

##Creating a unique phyloseq object with sample-by-sequence feature table, 
#the sample metadata and the sequence taxonomies 
ROHR_07_obj <- phyloseq(otu_table(tax_silva_otu, taxa_are_rows = FALSE), 
                       sample_data(ROHR_07_metadata), tax_table(tax_silva_tax)) 

#removing ASVs that are not present in any sample (if any)
ROHR_07_obj <- prune_taxa(taxa_sums(ROHR_07_obj) > 0, ROHR_07_obj)

#fix phyloesq object metadata so all samples are either 'Gut', 'Swab' or 'Control'
unique(ROHR_07_metadata$sample) #currently have an erroneous Stomach in the metadata

ROHR_07_metadata$sample[ROHR_07_metadata$sample == 'Stomach'] <- 'Gut' ##changes wrong label to Gut
unique(ROHR_07_metadata$sample) #now only have either Control, Gut, or Swab (caribbean)

new_sample_data <- sample_data(ROHR_07_metadata)

sample_data(ROHR_07_obj) <- new_sample_data

# Save Phyloseq object
saveRDS(ROHR_07_obj, "~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_07_trimmed/ROHR_07_obj.rds")

# Read Phyloseq object (if not starting from top)
ROHR_07_obj = readRDS("~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_07_trimmed/ROHR_07_obj.rds")


#### Identify & Clean Contaminants - Prevalence ####
sample_data(ROHR_07_obj)$is.neg <- sample_data(ROHR_07_obj)$sample == "Control"
contamdf.prev <- isContaminant(ROHR_07_obj, method="prevalence", neg="is.neg")
contaminants <- table(contamdf.prev$contaminant)
#check how many contaminants ID
contaminants # 27 true (cont), 14899 false (not cont)

write.csv(contamdf.prev, file='~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_07_trimmed/ROHR_07_contamdf.prev.csv')


# Make phyloseq object of presence-absence in negative controls and true samples
ROHR_07_obj.pa <- transform_sample_counts(ROHR_07_obj, function(abund) 1*(abund>0))
ROHR_07_obj.pa.neg <- prune_samples(sample_data(ROHR_07_obj.pa)$sample == "Control", ROHR_07_obj.pa)
ROHR_07_obj.pa.pos <- prune_samples(sample_data(ROHR_07_obj.pa)$sample != "Control", ROHR_07_obj.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ROHR_07_obj.pa.pos), pa.neg=taxa_sums(ROHR_07_obj.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")


# Create a new phyloseq object without contaminants identified in contamdf.prev
#surely better way to do this, but for now - view df, sort by "contaminants", select reads marked as "TRUE"
View(contamdf.prev)

badTaxa = c("TACGTAGGGTCCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGTTGTGCAAGACCGATGTGAAATCCCCGGGCTTAACCTGGGAATTGCATTGGTGACTGCACGGCTAGAGTGTGTCAGAGGGGGGTAGAATTCCACGTGTAGCAGTGAAATGCGTAGAGATGTGGAGGAATACCGATGGCGAAGGCAGCCCCCTGGGATAACACTGACGCTCATGCACGAAAGCGTGGGGAGCAAACAGG,
            TACGTAGGGGGCGAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGCGGTCTCTTAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGTCATTGGAAACTGGGGGACTTGAGGGCAGGAGAGGAGAGCGGAATTCCACGTGTAGCGGTGAAATGCGTAGAGATGTGGAGGAACACCAGTGGCGAAGGCGGCTCTCTGGCCTGCACCTGACGCTGAGGCGCGAAAGCGTGGGGAGCAAACAGG,
            TACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGGGAGCGCAGGTGGTTTCTTAAGTCTGATGTGAAAGCCCACGGCTTAACCGTGGAGGGTCATTGGAAACTGGGAAACTTGAGTACAGAAGAGGAATGTGGAACTCCATGTGTAGCGGTGGAATGCGTAGATATATGGAAGAACACCAGTGGCGAAGGCGACATTCTGGTCTGTTACTGACACTGAGGCTCGAAAGCGTGGGGAGCAAACAGG,
            TACGTAAGTCCCAAGCGTTACCCGGAATTATTGGGTGTAAAGGGTCCGTAGGCGGACTGGTAAGTCAGGTATGAAAGACCGGAGCTCAACTCCGAGTTTGTGCTTGAAACTGCAAGTCTTGAATCAGGGAGAGGTTAGCGGAATTCTAAGTGTAGGGGTGCAATCCGTAGATACTTAGAAGAACACCAAAAGCGAAGGCAGCTAACTGGAACTCGATTGACGCTGAGGGACGAAAGCGTGGGGAGCGAAGAGG,
            AGCAGCCGCGGTAATACGAGAGACTCAAGTTGACAGCCATCGGCGTAAAGTGTGGTCAAGATGAATAAAAACTAGAGCTGAATGCCCTCAAAGCTGTTATACGCCTTCGAAGGTAAGAAACCCAACTACGAAAGTGGCTCTATTATATCTGATTCCACGAAAGCTAGGGAACAAACTGGGATTAGATACCCCCGTAGTCCAGTAGATCGGAAGAGCACACGTCTGAACTCCAGT,
            TACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGCGGTCCTTTAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGTCATTGGAAACTGGAGGACTTGAGTACAGAAGAGGAAAGCGGAATTCCACGTGTAGCGGTGAAATGCGTAGAGATGTGGAGGAACACCAGTGGCGAAGGCGGCTTTCTGGTCTGTAACTGACGCTGAGGCGCGAAAGCGTGGGGAGCAAACAGG,
            TACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGAGCGTAGGTGGCTTGATAAGTCAGATGTGAAATCCCCGGGCTTAACCTGGGAACTGCATCTGATACTGTTAAGCTAGAGTAGGTGAGAGGGAAGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGATGGCGAAGGCAGCTTCCTGGCATCATACTGACACTGAGGCTCGAAAGCGTGGGTAGCAAACAGG,
            TACAGAGGGTGCGAGCGTTAATCGGATTTACTGGGCGTAAAGCGTGCGTAGGCGGTCATTTAAGTCAAATGTGAAATCCCCGAGCTTAACTTGGGAATTGCATTCGATACTGGATGGCTAGAGTATGGGAGAGGATGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGATGGCGAAGGCAGCCATCTGGCCTAATACTGACGCTGAGGTACGAAAGCATGGGGAGCAAACAGG,
            TACGGAGGGAGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGCTTTTCAAGTCAGGGGTGAAATCCCGGGGCTCAACCCCGGAACTGCCCTTGAAACTGGATGGCTAGAATACTGGAGAGGTGAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAAGAACACCAGTGGCGAAGGCGACTCACTGGACAGTTATTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGG,
            TACGTAGGTCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGTGGTTTATTAAGTCTGGTGTAAAAGGCAGTGGCTCAACCATTGTATGCATTGGAAACTGGTAGACTTGAGTGCAGGAGAGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCGGTGGCGAAAGCGGCTCTCTGGCCTGTAACTGACACTGAGGCTCGAAAGCGTGGGGAGCAAACAGG,
            GACGGGGGATGCAAGTGTTATCCGGAATCACTGGGCGTAAAGCGTCTGTAGGTGGTTTAGTAAGTCAACTGTTAAATCTTGAGGCTCAACCTCAAAATCGCAGTCGAAACTATTAGACTTGAGTATAGTAGGGGTAAAGGGAATTTCCAGTGGAGCGGTGAAATGCGTAGAGATTGGAAAGAACACCGATGGCGAAGGCACTTTACTGGGCTATTACTGACACTCAGAGACGAAAGCTAGGGTATCAAATGGG,
            TACGAAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGATAAGTTGGATGTGAAAGCCCCGGGCTCAACCTGGGAATTGCATCCAAAACTGTCTGACTAGAGTATGGCAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGGCTAATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAAACAGG,
            TACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGCGGCCTTTTAAGTCTGATGTGAAAGCCCCCGGCTCAACCGGGGAGGGCCATTGGAAACTGGAAGGCTTGAGTACAGAAGAGAAGAGTGGAATTCCACGTGTAGCGGTGAAATGCGTAGAGATGTGGAGGAACACCAGTGGCGAAGGCGACTCTTTGGTCTGTAACTGACGCTGAGGCGCGAAAGCGTGGGGAGCAAACAGG,
            TACGTAGGTGGCGAGCGTTATCCGGAATTATTGGGCGTAAAGGGAGCGTAGGCGGATACTTAAGTGGGATGTGAAAGACTTGGGCTCAACCCGAGGGCTGCATTCCAAACTGAGTATCTAGAGTGCAGGAGAGGAGAGTGGAATTCCTAGTGTAGCGGTGAAATGCGTAGAGATTAGGAAGAACACCAGTGGCGAAGGCGACTCTCTGGACTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAGG,
            TACGCAGGTGGCGAGCGTTGCCCGGAATTACTGGGCGTAAAGGGTGCGTAGGCGGCCGGACAAGTCATAGGTTAAAGCCCGGAGCTCAACTCCGGAAAGGCCTATGATACTGTCTGGCTTGAGGGCCGGAGAGGCTGGCGGAATTCCCGGTGTAGGGGTGAAATCCGTAGATATCGGGAGGAACACCGGTGGGGAAGCCGGCCAGCTGGACGGTTCCTGACGCTGAGGCACGAAAGCGTGGGGAGCAAACCGG,
            TACAGAGGGTGCGAGCGTTAATCGGATTTACTGGGCGTAAAGCGTGCGTAGGCGGCTGATTAAGTCGGATGTGAAATCCCTGAGCTTAACTTAGGAATTGCATTCGATACTGGTCAGCTAGAGTATGGGAGAGGATGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGATGGCGAAGGCAGCCATCTGGCCTAATACTGACGCTGAGGTACGAAAGCATGGGGAGCAAACAGG,
            TACGTAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGTTTTGTAAGTCTGATGTGAAATCCCCGGGCTCAACCTGGGAATTGCATTGGAGACTGCAAGGCTAGAATCTGGCAGAGGGGGGTAGAATTCCACGTGTAGCAGTGAAATGCGTAGATATGTGGAGGAACACCGATGGCGAAGGCAGCCCCCTGGGTCAAGATTGACGCTCATGCACGAAAGCGTGGGGAGCAAACAGG,
            TACGTAGGGGGCGAGCGTTGTCCGGATTCACTGGGCGTAAAGCGCGCGCAGACGGCTTTCTGCGTCTGGGGTGAAAACCCGGGGCTCAACCCCGGGAGTGCCTTGGATACGGGAGAGCTTGAGGGTCGGAGAGGCAAGGGGAATTCCACGTGTAGCGGTGAAATGCGTAGAGATGTGGAAGAACACCTGTGGCGAAAGCGCCTTGCTGGCCGATCCCTGACGTTGAGGCGCGAAAGCGTGGGGAGCAAACAGG,
            TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTCTGTTAAGTCAGATGTGAAATCCCCGGGCTTAACCTGGGAACTGCATTTGAAACTGGCAGGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGG,
            TACGGAGGGCGCGAGCGTTACCCGGATTTACTGGGCGTAAAGGGCGTGTAGGCGGCTTGGGGCGTCCCATGTGAAAGACCACGGCTCAACCGTGGGGGAGCGTGGGATACGCTCAGGCTAGACGGCGGGAGGGGGTGGTGGAATTCCCGGAGTAGCGGTGAAATGCGCAGATACCGGGAGGAACGCCGATGGCGAAGGCAGCCACCTGGCTCGTTCGTGACGCTGAGGCGCGAAAGCGTGGGGAGCAAACCGG,
            TACGAAGGGGGCTAGCGTTGCTCGGAATTACTGGGCGTAAAGGGCGCGTAGGCGGATCGTTAAGTCAGAGGTGAAATCCCAGGGCTCAACCCTGGAACTGCCTTTGATACTGGCGATCTTGAGTATGAGAGAGGTATGTGGAACTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAAGAACACCAGTGGCGAAGGCGACATACTGGCTCATTACTGACGCTGAGGCGCGAAAGCGTGGGGAGCAAACAGG,
            TACGTAGGGTGCGAGCGTTGTCCGGAATTACTGGGCGTAAAGAGCTCGTAGGTGGTTTGTCGCGTTGTTCGTGAAAACTCACAGCTTAACTGTGGGCGTGCGGGCGATACGGGCAGACTAGAGTACTGCAGGGGAGACTGGAATTCCTGGTGTAGCGGTGGAATGCGCAGATATCAGGAGGAACACCGGTGGCGAAGGCGGGTCTCTGGGCAGTAACTGACGCTGAGGAGCGAAAGCGTGGGGAGCGAACAGG,
            TACGGAGGGGGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGACTGGAAAGTTGGGGGTGAAATCCCGGGGCTCAACCTCGGAACTGCCTCCAAAACTATCAGTCTGGAGTTCGAGAGAGGTGAGTGGAATACCGAGTGTAGAGGTGAAATTCGTAGATATTCGGTGGAACACCAGTGGCGAAGGCGGCTCACTGGCTCGATACTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGG,
            TACGTAGGGTGCGAGCGTTGTCCGGATTTATTGGGCGTAAAGGGCTCGTAGGTGGTTGATCGCGTCGGAAGTGTAATCTTGGGGCTTAACCCTGAGCGTGCTTTCGATACGGGTTGACTTGAGGAAGGTAGGGGAGAATGGAATTCCTGGTGGAGCGGTGGAATGCGCAGATATCAGGAGGAACACCAGTGGCGAAGGCGGTTCTCTGGGCCTTTCCTGACGCTGAGGAGCGAAAGCGTGGGGAGCGAACAGG,
            TACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCTTAAAGAGTTCGTAGGTGGTTGAAAAAGTTGGTGGTGAAATCCCAGAGCTTAACTCTGGAACTGCCATCAAAACTTTTCAGCTAGAGTTTGATAGAGGAAAGCAGAATTTCTAGTGTAGAGGTGAAATTCGTAGATATTAGAAAGAATACCAATTGCGAAGGCAGCTTTCTGGATCAATACTGACACTGAGGAACGAAAGCATGGGTAGCGAAGAGG,
            GACGGAGGATGCAAGTGTTATCCGGAATCACTGGGCGTAGAGCGTCTGTAGGTGGTCAAATAAGTCAACTGTTAAATCTTGAGGCTCAACCTCAAAACCGCAGTCGAAACTATTAGACTAGAGTATAGTAGGGGTAAAGGGAATTTCCAGTGGAGCGGTGAAATGCGTAGAGATTGGAAGGAACACCGATGGCGAAGGCACTTTACTGGGCTATTACTAACACTCAGAGACGAAAGCTAGGGTAGCAAATGGG,
            TACGAAGGTCCCGAGCGTTATTCGGAATCACTGGGCGTAAAGGGAGCGTAGGCTGCGCGGTAAGTCAGATGTGAAATCTCAGGGCTCAACCCTGAAACTGCATCCGATACTGCCGTGCTAGAGTAATGGAGGGGTAACTGGAATTCTCGGTGTAGCAGTGAAATGCGTGGATATCGAGAGGAAGACCAATGGCGAAGGCAAGTTACTGGACATTTACTGACGCTGAGGCTCGAAGGCTAGGGTAGCGAAAGGG")

goodTaxa <- setdiff(taxa_names(ROHR_07_obj), badTaxa)
ROHR_07_obj <- prune_taxa(goodTaxa, ROHR_07_obj)


# Now remove all control samples
##is this correct? did I accidentally get rid of a bunch of other sample sequences?
ROHR_07_obj = subset_samples(ROHR_07_obj, sample_data(ROHR_07_obj) != "Control")

#check that all controls were removed (should not have any control-ext, control-pcr)
sample_names(ROHR_07_obj) #still have the control samples that had RRR numbers as names

#remove remaining controls
ROHR_07_obj = subset_samples(ROHR_07_obj, sample_names(ROHR_07_obj) != "16860")
ROHR_07_obj = subset_samples(ROHR_07_obj, sample_names(ROHR_07_obj) != "16835")
ROHR_07_obj = subset_samples(ROHR_07_obj, sample_names(ROHR_07_obj) != "16890")
ROHR_07_obj = subset_samples(ROHR_07_obj, sample_names(ROHR_07_obj) != "16870")

#check that all controls were removed (should not have any control-ext, control-pcr)
sample_names(ROHR_07_obj) #should have 372 total samples, given that 12/384 were controls
#--> all good to go

# Look at distribution of library size
df <- as.data.frame(sample_data(ROHR_07_obj)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ROHR_07_obj)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color = sample)) + geom_point()

# Remove outlier samples ##remove any??
ROHR_07_obj = subset_samples(ROHR_07_obj, sample_names(ROHR_07_obj) != "")


##Removing ASVs of chloroplasts
ROHR_07_obj <- subset_taxa(ROHR_07_obj, Order != "Chloroplast" | is.na(Order))
ROHR_07_obj

#otu_table()   OTU Table:         [ 14555 taxa and 372 samples ]
#sample_data() Sample Data:       [ 372 samples by 44 sample variables ]
#tax_table()   Taxonomy Table:    [ 14555 taxa by 6 taxonomic ranks ]

##Removing ASVs of mitochondria
ROHR_07_obj <- subset_taxa(ROHR_07_obj, Family != "Mitochondria" | is.na(Family))
ROHR_07_obj

#otu_table()   OTU Table:         [ 13438 taxa and 372 samples ]
#sample_data() Sample Data:       [ 372 samples by 44 sample variables ]
#tax_table()   Taxonomy Table:    [ 13438 taxa by 6 taxonomic ranks ]

##Removing ASVs assigned to Eukaryotes and unassigned (NA) at the kingdom level
ROHR_07_obj <- subset_taxa(ROHR_07_obj, Kingdom != "Eukaryota")
ROHR_07_obj

#otu_table()   OTU Table:         [ 13415 taxa and 372 samples ]
#sample_data() Sample Data:       [ 372 samples by 44 sample variables ]
#tax_table()   Taxonomy Table:    [ 13415 taxa by 6 taxonomic ranks ]

ROHR_07_obj <- prune_taxa(taxa_sums(ROHR_07_obj) > 0, ROHR_07_obj)
ROHR_07_obj
#otu_table()   OTU Table:         [ 13389 taxa and 372 samples ]
#sample_data() Sample Data:       [ 372 samples by 44 sample variables ]
#tax_table()   Taxonomy Table:    [ 13389 taxa by 6 taxonomic ranks ]


#check what taxonomic classifications are in the phyloseq object (optional)
rank_names(ROHR_07_obj)
  # have "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"  

#list out sample variables (taken from metadata file - includes RRR, ocean, region, site...) (optional)
sample_variables(ROHR_07_obj)

# Save Phyloseq object
saveRDS(ROHR_07_obj, "~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_07_trimmed/ROHR_07_obj_clean.rds")

##############################################################################################

# Read Phyloseq object
ROHR_07_obj_clean = readRDS("~/OneDrive - McGill University/McGill/Lab_Docs/Fish_Data/Fish_2022/2021-2022_TEP_Fish/16S_Summer_2022/ROHR_07_trimmed/ROHR_07_obj_clean.rds")

##############################################################################################




