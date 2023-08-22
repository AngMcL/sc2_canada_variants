## 01_subSampled

### Ins and Outs
- Input: cleaned, partitioned alignments + corresponding metadata
- No GISAID sequences or metadata included here.
- Output: subsampled alignments and metadata

### Scripts:
*00_MergePartitionsSplitCanada.Rmd*
- Merges the partitioned sequence and metadata, then splits into Canada and global sets

*01_SummarizeCasesSeqsByGeoLin.Rmd*
- Estimates variant cases over time by Canadian province or global region
- Summarizes raw and proportional contributions to cases, sequences pre-subsampling

*02_SubsampleByVOC.Rmd*
- Or by Omicron lineage, BA.1, BA.1.1, BA.2
- Based on a total number of target sequences and sequences available,temporally distributes the  number of sequences to sample from each month to mimic the distribution of variant cases in Canada or globally while also up-sampling sparse early months
- Uses the proportional contributions of Canadian provinces or global regions to variant cases to inform sequences' sampling probabilities
- Subsample 10 bootstraps per variant

*02b_SubsampleByVOC_uniform.Rmd*
- Or by Omicron lineage, BA.1, BA.1.1, BA.2
- Based on a total number of target sequences and sequences available, uniformly distributes the number of sequences to sample from each month to be maximally uniform in the confines of what is available (strategy used for first two waves of SARS-CoV-2 in Canada analysis)
- Uses the proportional contributions of Canadian provinces or global regions to variant cases to inform sequences' sampling probabilities

*03_SummarizeSubsamples.Rmd*
- Summarizes raw and proportional contributions to sequences post-subsampling

*03b_SummarizeSubsample_uniform.Rmd*
- Summarizes raw and proportional contributions to sequences post-subsampling with uniform temporal distribution of sequences

*04_ComparePearsons.Rmd*
- compare subsampling strategies in terms of correlation between cases and sequences overall and by month
1) Uniform: similar to previous paper, uniformly distributed sequences over time
2) Temporally distributed: seqs over time mimic variant cases over time
3) Pre-subsample

*05_QueryFirstSample.Rmd*
- Make sure the first Canadian GISAID sample for each variant (in cleaned data) is accurate 

*06_GenerateEpiSets.Rmd*
- generate gisaid episet for all the uniquely sampled sequences
- generate lists of accession IDs sampled in each bootstrap

*07_Fig1CasesImports.Rmd*
- make a nicely formatted fig 1 of cases and importations by variant over time with cartoon phylo

