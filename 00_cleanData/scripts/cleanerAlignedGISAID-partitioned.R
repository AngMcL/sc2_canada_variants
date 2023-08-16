#!/usr/bin/env Rscript

# objective
# trimmed down version of the Rmd to clean a set of partitioned and aligned fasta files
# usage:
# Rscript scripts/cleanerAlignedGISAID-partitioned.R <input fasta folder> <input meta folder>
# Rscript scripts/cleanerAlignedGISAID-partitioned.R "GISAID/sequences/02_aligned" "GISAID/metadata/01_clean"


## Load libraries
library(ape)
library(stringr)
library(stringi)
library(Biostrings)
library(gtools)
library(dplyr)
library(tidyr)
library(plyr)

## Inputs and outputs here
args<-commandArgs(trailingOnly=TRUE)
seq.fold.in<-args[1] #folder with sequences
meta.fold.in<-args[2] #folder with metadata

#manual run
seq.fold.in<-"GISAID/sequences/02_aligned"
meta.fold.in<-"GISAID/metadata/01_clean"

#list of files input
list.seq.in<-list.files(seq.fold.in,pattern = ".fasta.aln",full.names = T) %>% mixedsort()
list.meta.in<-list.files(meta.fold.in,pattern = "clean.csv",full.names = T) %>% mixedsort()
if(length(list.seq.in)!=length(list.meta.in)){print("meta and seqs don't match");break}
n.part<-length(list.seq.in)

#make output folders
seq.fold.out<-str_replace(seq.fold.in,"02_aligned","03_cleaner")
meta.fold.out<-str_replace(meta.fold.in,"01_clean","03_cleaner")
if(!dir.exists(seq.fold.out)){dir.create(seq.fold.out)}
if(!dir.exists(meta.fold.out)){dir.create(meta.fold.out)}

#output files
list.seq.out<-list.seq.in %>% str_replace_all(seq.fold.in,seq.fold.out) %>% 
  str_replace_all("clean.fasta.aln","cleaner.fasta")
list.meta.out<-list.meta.in %>% str_replace_all(meta.fold.in,meta.fold.out) %>% 
  str_replace_all("clean.csv","cleaner.csv")
list.summary.out<-str_replace_all(list.meta.out, ".csv","_ExclSummary.csv")

#loop through each partition
for (k in 1:n.part){
  print(paste("Reading in partition",k,sep=" "))
  
  # Curate and clean sequences
  ## Read in the GISAID metadata
  meta<-read.csv(file = list.meta.in[k],header = TRUE)
  
  ## Read alignment
  file.seqs<-readDNAStringSet(list.seq.in[k])

  #size of alignment 
  l.og<-length(file.seqs)
  l.can<-length(str_which(meta$Location, "Canada"))
  
  # Quantify the actg,n,- content and exclude low quality sequences
  #summarize the content of each sequence
  seq.summary<-data.frame(alphabetFrequency(file.seqs))
  seq.summary$length<-29903
  colnames(seq.summary)[16:18]<-c("-","+",".")
  seq.summary<-seq.summary %>% mutate(
    AMBIG = M+R+W+S+Y+K+V+H+D+B,
    a=A/length,
    c=C/length,
    t=T/length,
    g=G/length,
    n=N/length,
    ambig=AMBIG/length,
    gap=`-`/length)
  
  #Identify sequences with >10% gaps or >20% n >10% ambig
  bad.gaps<-which(seq.summary$gap>0.10) #removes sequences with <26912.7 b
  bad.n<-which(seq.summary$n>0.20) #had to bump this up to include the BC sequences...
  bad.ambig<-which(seq.summary$ambig>0.10)
  
  #how many? should match table
  l.bad.n<-length(bad.n)
  l.bad.gaps<-length(bad.gaps)
  l.bad.ambig<-length(bad.ambig)
  #l.bad.n; l.bad.gaps; l.bad.ambig
  
  #merge into one vector
  gaps.ns<-unique(c(bad.n,bad.gaps,bad.ambig))
  
  #REMOVE THEM from alignment
  if (length(gaps.ns)>0){
    file.seqs<-file.seqs[-gaps.ns]
  }
  
  #restrict meta to remaining sequence names
  meta<-meta[which(meta$new.names %in% names(file.seqs)),]
  #get rid of reference in alignment
  file.seqs<-file.seqs[which(names(file.seqs) %in% meta$new.names)]
  l.final<-length(file.seqs)
  l.final.can<-length(str_which(meta$Location,"Canada"))
  l.final==nrow(meta)
  
  #print summary statistics of the cleaning run to a csv file
  summary.out<-str_replace_all(list.meta.out[k], ".csv","_gapn-summary.csv")
  out.colz<-c("n.seq.in","n.seq.in.can","n.badn","n.badgaps","n.ambig","n.gapns",
              "n.seq.out","n.seq.out.canada")
  valuez<-c(l.og, #original length
            l.can, #original number from Canada
            l.bad.n, #>20% n
            l.bad.gaps, #>10% gaps
            l.bad.ambig, #>10% ambig
            length(gaps.ns), #all three above
            l.final, #final number of sequences
            l.final.can #final number of Can sequences
  )
  summary.statz<-data.frame(parameter=out.colz,value=valuez)
  write.csv(summary.statz,file=summary.out,row.names = F)
  
  #export META
  write.csv(meta,file=list.meta.out[k],row.names = F)
  
  #export FASTA - note that have to resave using ape for proper format no gaps - issue with biostrings...
  print("Writing output")
  writeXStringSet(file.seqs,file=list.seq.out[k], format="fasta")
  rm(file.seqs)
  fix<-read.FASTA(list.seq.out[k])
  write.FASTA(fix,list.seq.out[k])
  rm(fix)
}
