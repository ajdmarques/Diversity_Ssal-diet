# COI metabarcoding provides insights into the high diverse diet of a generalist salamander, _Salamandra salamandra_

Raw data and R scripts are provided to replicate the results of in Marques, Mata, & Velo-Antón (2022) COI Metabarcoding Provides Insights into the Highly Diverse Diet of a Generalist Salamander, Salamandra salamandra (Caudata: Salamandridae) published in Diversity.

Users should download the contents of this repository. Open the R project 'DietAnalysis.Rproj' and run the R scripts in order from 01 to 04.
Input data (./input/_primer_/) includes read counts for finalized operational taxonomic units (OTUs) '_primer_.read-count.csv, taxonomic assignments as determined by BOLDIGGER '_primer_.boldigger.csv', and sample information '_primer_.sample-info.csv'. 

Scripts can be found in './R-scripts/' and in order of operation will, 
  01) merge input files from different fragments into a single table and apply filters, 
  02) calculate and generate output figures of sample rarefaction curves for each fagments and region,
  03) compare prey richness between fragments and regions,
  04) build a PerMANOVA model and test for differences between fragments and regions.

Before proceeding, ensure the following packages have been installed.

install.packages(c('dplyr','tidyverse','iNEXT','vegan','magrittr','ggplot2'))

Please cite as follows:
Marques, A.J.D.; Mata, V.A.; Velo-Antón, G. COI Metabarcoding Provides Insights into the Highly Diverse Diet of a Generalist Salamander, Salamandra salamandra (Caudata: Salamandridae). Diversity 2022, 14, 89. https://doi.org/10.3390/d14020089 



