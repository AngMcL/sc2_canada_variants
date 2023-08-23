## 02_trees

In the interest of not breaching GISAID data sharing agreement, we have not included the SARS-CoV-2 phylogenies for this analysis. Trees are available upon request to GISAID approved data users. I have, however, included scripts needed to infer, root, remove temporal outliers, and time scale trees.

I have also included a brief side analysis comparing the strict and relaxed molecular clock rates of variants of concern and interest. In this analysis, I also compared clock rates for wild type SC2 from the first two waves (McLaughlin et al, eLife, 2022) when 25, 50, 75, or 100% of Canadian sequences were sampled, with the difference up to 50,000 from global sequences. (Relatively more Canadian sequences means less diversity and generally slower clock). 

In the results folder, I have included summaries of the removal of temporal outliers. The temporal signal including outliers prior to their removal (InclOutliers) and following their removal (ExclOutliers) were compared.

### Scripts
*00_IterateFasttree.sh*
- A shell script to iteratively run fasttree over all alignments in a folder

*01_Root.am.R*
- R script to outgroup root trees on Wuhan-hu-1

*02_RemoveOutliers.Variant.R*
- R script to identify and remove outliers if their residuals exceeded 0.001 in the relationship of divergence over time (excluding parental lineages from the linear fit) or if they had more than 30 mutations on their pendant edge (terminal branch leading to tip)
- Generates summary plots in results folder

*02_RemoveOutliers.VOC-omi.R*
- As above, but for separate omicron lineages

*03_IterateIQtree.2.sh*
- A shell script to iteratively run IQ-Tree2 with LSD2 to timescale trees under a relaxed molecular clock with 0.2 variance and a GTR+I+R3 substitution model

*04_ExtractSlopesTMRCAFromLog.Rmd*
- R markdown script to extract from the IQTree logfile (all available upon request) the slope (relaxed molecular clock rate)and TMRCA of the tree