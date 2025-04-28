# Cryptic biodiversity: untangling coral reef fishes' microbiomes in changing environments

## TLDR
Repository for fish skin microbiome work (16S), to be included in an upcoming manuscript (in prep as of April 2025)

This repository holds code for processing 16S sequences of fish skin microbiome samples collected from tropical reef fish in the Tropical Eastern Pacific of Panama, as well as associated environmental water samples. 

Current samples were collected in 2021-22 by members of the Rohr Reef Resilience project at the Smithsonian Tropical Research Institute. 

## Authors

- [@llardinois](https://www.github.com/lardinois21)


## Project Details

Fish were collected from 2 sites in the Tropical Eastern Pacific: Las Perlas & Coiba. Sampling was done during both upwelling & non-upwelling seasons (collections done 2021-22). We had 12 target species, 10 of which are represented in this dataset. Skin swabs = 9 per species. (paired Gut samples = 8 per species, were also collected, but will be part of a separate project, stay tuned...) Control swabs (swabs exposed to air in field) were collected at each site and season. Swabs were stored in liquid nitrogen for transport from field to lab, then at -80ºC. We also collected environmental data (temperature, DO) from each site.

## Project Scripts


(1) Skin Pre-Processing RRR Fish 16S Version 1.Rmd 
    - sequence pre-processing (following cutadapt step in terminal & creation of initial phyloseq object with tax table, OTU table, and metadata

