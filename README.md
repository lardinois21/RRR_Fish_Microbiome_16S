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

Note: Start with scripts #1 & #2, which perform basic metabarcoding data cleaning and append a phylogenetic tree, then use any of the other scripts as needed, depending on desired analyses and visualization. The list below is given in order of typical use - environmental data can be run independently. Note that all analyses were run on both unrarefied and rarefied data.

### Fish Microbiome Scripts
(1) Skin Pre-Processing RRR Fish 16S Version 1.R
 - sequence pre-processing (following cutadapt step in terminal), quality control & creation of initial phyloseq object with tax table, OTU table, and metadata
 - follows DADA2 pipeline

(2) Skin Alpha Diversity 16S Version 1.Rmd
 - additional cleaning & visualizing of phyloseq object
 - alpha diversity analyses, including Hill numbers
 - calculating phylogenetic tree for microbial ASVs & adding to new phyloseq
 - rarefying dataset
 - creating various a-div plots
    
(3) Skin Beta Diversity 16S Version 1.Rmd
- beta diversity metrics run on full fish dataset (many host species)
- Bray Curtis, Jaccard, Unifrac, WUnifrac... = different diversity metrics run to assess microbial community composition.
- PERMANOVAs to test significance of different predictors in explaining microbial community variation across species, regions, and seasons
- ordinations (PCoAs, NMDS) to visualize microbiome (dis)similarity across groups

(4) Skin Species by Species Version 1.Rmd
- scripts splits whole dataset into subsets by host species, for both the rarefied and unrarefied datasets
- beta diversity metrics run on each host species dataset [ can skip to script below, which makes this run more quickly by looping ]

(5) Skin Species by Species B div code V1.Rmd
- beta diversity metrics run species-by-species
- runs loops on species-by-species datasets created in script #4 to test various b-div metrics and visualizations, matching those run on the whole dataset in script #3
    
(6) Species by Species Differential Abundance Version 1.Rmd
- differential abundance analyses (DESeq2) run on species with significant differences in microbiome between seasons (from script #5) to look at which microbial taxa are changing and driving these community level differences

(7) Skin Figures Version 1.Rmd
- script used to refine and clean manuscript figures, especially the multi-panel plots (some figures coded directly in scripts above)

### Water Microbiome Scripts
Note: first run water microbiome samples through the same pre-processing script (Fish script #1) and/or merge the datasets (fish + water) prior to pre-processing and run everything together

(1) Water Alpha Diversity 16S Version 1.Rmd
- Same as fish above, but with water samples:
        - additional cleaning & visualizing of phyloseq object
        -  alpha diversity analyses, including Hill numbers
        - calculating phylogenetic tree for microbial ASVs & adding to new phyloseq
        - rarefying dataset
        - creating various a-div plots
        
(2) Water Beta Diversity 16S Version 1.Rmd
- Same as fish above, but with water samples:
        - Bray Curtis, Jaccard, Unifrac, WUnifrac... = different diversity metrics run to assess microbial community composition.
        - PERMANOVAs to test significance of different predictors in explaining microbial community variation across species, regions, and seasons
        - ordinations (PCoAs, NMDS) to visualize microbiome (dis)similarity across groups
        
(3) Water Differential Abundance 16S Version 1.Rmd
- Same as fish above, but with water samples & comparing water to fish:
        - 2 differential abundance analyses (MaAsLin2 and DESeq2) tested on water samples

### Environmental Data Scripts
(1) Environmental Data Temp and DO Version 1.Rmd
- cleans HOBO logger temperature and dissolved oxygen data from the RRR project
- takes averages of publicly available data from the STRI Environmental Monitoring Program and adds these to the RRR data
- creates plots for temperature, DO, and the two together




