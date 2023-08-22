# sc2_canada_variants
Angela McLaughlin et al. 

## Description
Scripts needed to reproduce our analysis as described in the manuscript, 'Effectiveness of Canadian travel restrictions in reducing burden of SARS-CoV-2 variants of concern'
pre-print to be made available at [medRxiv](https://www.medrxiv.org).

## Dependencies
Note this pipeline was developed and tested on mac OX Monterey 12.15.1 

### Programs
* R version 4.3.1
* Python 3.8.3
* pangolin v3.1.20 and and pangoLEARN data release March 22, 2022 https://github.com/cov-lineages/pangolin
* FastTree v2.1.11
* IQTREE-2.1.2 (includes LSD2)
* minimap2 

### R packages
ape 5.7-1, Biostrings 2.68.1, phytools 1.9-16, phangorn 2.11.1, coronavirus 0.4.1, tidyverse 2.0.0, plyr 1.8.8, lubridate 1.9.2, stringi 1.7.12, zoo 1.8-12, ochRe 1.0.0, cowplot 1.1.1, ggstance 0.3.6, ggalluvial 0.12.5, ggmosaic 0.3.3, ggtree 3.8.2, ggplotify 0.1.1, ggrepel 0.9.3, and MASS 7.3-60. 

### Python packages/scripts
* viralMSA.py https://github.com/niemasd/ViralMSA
* masking script https://github.com/W-L/ProblematicSites_SARS-CoV2/blob/master/src/mask_alignment_using_vcf.py

### Additional files cloned from github
* nextstrain exclude list https://github.com/nextstrain/ncov/blob/master/defaults/exclude.txt
* problematic sites https://github.com/W-L/ProblematicSites_SARS-CoV2/blob/master/problematic_sites_sarsCov2.vcf

## Usage
The readme files in the individual subfolders have details on running scripts. This repo does not include a reproducible example with fake data, as in the previous repo for the first two waves of SARS-CoV-2, and GISAID data cannot be shared publicly, therefore inputs reference files not included in the repo. Inputs and output would have to be modified to run the script locally.

## Thank you
* GISAID.org
* Originating and submitting laboratories and personnel who contributed to the GISAID database
* Co-authors: Vincent Montoya, Rachel Miller, Canadian COVID-19 Genomics Network (CanCOGeN) Consortium, Michael Worobey, Jeffrey B. Joy
* The many software and package developers
* See manuscript for full citations

## Funding sources
AM was supported by a Canadian Institutes for Health Research (CIHR) Doctoral grant.  MW was supported by the David and Lucile Packard Foundation. JBJ was supported by Genome Canada BCB 287PHY grant, an operating grant from the CIHR Coronavirus Rapid Response Programme number 440371, and a CIHR variant of concern supplement. The British Columbia Centre for Excellence in HIV/AIDS also provided support for AM, VM, RLM, and JBJ.

