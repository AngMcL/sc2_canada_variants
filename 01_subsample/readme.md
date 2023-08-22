## 01_subSampled

Input: cleaned, partitioned alignments + corresponding metadata
Output: subsampled alignments and metadata

No GISAID sequences or metadata included here.

Order of scripts:
00_MergePartitionsSplitCanada.Rmd
01_SummarizeCasesSeqsByGeoLin.Rmd
02_SubsampleByVariant.Rmd
03_SummarizeSubsamples.Rmd


Subsampling strategies to compare:
1) Uniform: similar to previous paper, uniformly distributed sequences over time
2) Temporally distributed: seqs over time mimic variant cases over time

x 10 bootstraps each

