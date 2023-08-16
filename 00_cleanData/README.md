# 00_cleanData

* Note that we can not share GISAID sequences or metadata with unauthorized users, therefore have excluded any raw data here

* Before starting:
    * cp the exclude.txt file from nextstrain ncov
	* Raw data into GISAID/00_original
	
* Run scripts/cleanUnalignedGISAID-partition.R
	* Remove non-human and environmental sequences
    * Clean sequence names, remove seqs with incomplete dates from outside Canada
    * output into GISAID/sequences/01_clean/ and GISAID/meta/01_clean/
    * NOTE mods made here to loop across clade partitions for computational feasibility
    
* Aligned using viralMSA (python package by Niema Moshiri) invoking minimap2 into GISAID/sequences/02_aligned

* Run scripts/cleanAlignedGISAID.R
    * remove seqs with gaps, Ns, or ambigs >10%
    * Output GISAID/sequences/03_cleaner/ and GISAID/metadata/03_cleaner/
    
* Seems to be an issue with Biostrings because the output fasta had gaps inserted every 500 or so bases
	* fixed by read.FASTA, write.FASTA the exported filed from biostring
	* for biggest fasta, re-align the data

* Delay masking problematic sites until after subsample (pre-filter job was killed b/c too big)

