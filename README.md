# Cryptic biodiversity: untangling coral reef fishes' microbiomes in changing environments

## TLDR
Repository for fish skin microbiome work (16S), to be included in an upcoming manuscript (stay tuned for pre-print soon!)

This repository holds code for processing, analyzing, and visualizing 16S sequences of fish skin microbiome samples collected from tropical reef fish in the Tropical Eastern Pacific of Panama, as well as associated environmental water samples. There is also code to parse through various types of environmental data, including water temperature and dissolved oxygen, collected in various sites in Panama.

These scripts should work for any kind of 16S metabarcoding data, and can easily be adapted to work with other primers/DNA barcodes (i.e. 18S, 23S), with some modifications at the taxonomic assignment step in sequence pre-processing. The alpha and beta diversity analyses are adapted for microbiome work from classic community ecology metrics, following current practices for microbial ecology. 

Current samples were collected in 2021-22 by members of the Rohr Reef Resilience project at the Smithsonian Tropical Research Institute in Panama. 


## Authors

- [@llardinois](https://www.github.com/lardinois21)


## Project Details

Fish were collected from 2 island archipelagos in the Tropical Eastern Pacific: Las Perlas & Coiba. Sampling was done during both upwelling & non-upwelling seasons (collections done 2021-22). We had 12 target species, 10 of which are represented in this dataset. Skin swabs = 9 per species/season/region. (paired Gut samples = 8 per species, were also collected, but will be part of a separate project, stay tuned...) Control swabs (swabs exposed to air in field) were collected at each site and season. Swabs were stored in liquid nitrogen for transport from field to lab, then at -80ºC. We also collected environmental data (temperature, DO) from each site.

## Project Scripts

#Note:# Start with scripts #1 & #2, which perform basic metabarcoding data cleaning and append a phylogenetic tree, then use any of the other scripts as needed, depending on desired analyses and visualization. The list below is given in order of typical use - environmental data can be run independently.

(1) Skin Pre-Processing RRR Fish 16S Version 1.R
    - sequence pre-processing (following cutadapt step in terminal), quality control & creation of initial phyloseq object with tax table, OTU table, and metadata
    - follows DADA2 pipeline

(2)

(3)

(4)

(5)



