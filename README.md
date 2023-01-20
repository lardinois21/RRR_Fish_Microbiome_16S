# RRR_Fish_Microbiome_16S
 Repository for fish microbiome work (16S)
This repository holds code for processing 16S sequences of fish skin and gut microbiome samples collected from tropical reef fish in Panama, as well as associated environmental water samples. 
Current samples were collected in 2021-22 by members of the Rohr Reef Resilience project at the Smithsonian Tropical Research Institute. 

Project Summary:

Fish were collected from 2 sites in the Tropical Eastern Pacific: Las Perlas & Coiba. Sampling was done during both upwelling & non-upwelling seasons (collections done 2021-22). We had 12 target species, 10 of which are represented in this dataset. Gut samples = 8 per species. Skin swabs = 9 per species. Control swabs (swabs exposed to air in field) were collected at each site and season. Swabs were stored in liquid nitrogen for transport from field to lab, then at -80ºC. 

DNA was extracted using ZymoBIOMICS 96 well DNA kit - 4 plates total (per sample type = skin & gut). DNA was amplified following a modified version of the Earth Microbiome Project 16S protocol using phased primers (1 set per plate). Negative PCR controls and extraction controls were included in each plate. Index PCR performed to attach unique barcodes, then all four plates (384 samples) pooled for sequencing on the Illumina MiSeq sequencing platform of the Smithsonian Tropical Research Institute's (STRI) Naos facilities in Panama (courtesy of Marta Vargas). Sequenced libraries were demultiplexed at STRI by Marta Vargas using the MiSeq Reporter Software. These libraries were sent to me (Laura Lardinois) and trimmed with Cutadapt to remove primers before being read into R to start the DADA2 pipeline (following tutorial v 1.16). Silva version 138.1 (mar 10, 2021 update) used for taxonomic assignments.

See "Pre Processing RRR Fish 16S Sequences Summer 2022" documents for pipeline and processing up to the creation of the initial phyloseq object.