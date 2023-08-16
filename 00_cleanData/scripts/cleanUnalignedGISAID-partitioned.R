#!/usr/bin/env Rscript

# objective
# trimmed down version of the Rmd to clean a set of partitioned (into clades) and unaligned fasta files

# Rscript <myscript.R> <input_seq.fasta> <input_meta.tsv>
# Rscript scripts/cleanUnalignedGISAID-partitioned.R GISAID/sequences GISAID/metadata 

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
# seq.fold.in<-"GISAID/sequences"
# meta.fold.in<-"GISAID/metadata"

#list of files input
list.seq.in<-list.files(seq.fold.in,pattern = ".fasta",full.names = T) %>% mixedsort()
list.meta.in<-list.files(meta.fold.in,pattern = ".tsv",full.names = T) %>% mixedsort()
if(length(list.seq.in)!=length(list.meta.in)){print("meta and seqs don't match");break}
n.part<-length(list.seq.in)

#make output folders
seq.fold.out<-paste(seq.fold.in,"clean",sep="/")
meta.fold.out<-paste(meta.fold.in,"clean",sep="/")
if(!dir.exists(seq.fold.out)){dir.create(seq.fold.out)}
if(!dir.exists(meta.fold.out)){dir.create(meta.fold.out)}

#output files
list.seq.out<-list.seq.in %>% str_replace_all(seq.fold.in,seq.fold.out) %>% str_replace_all(".fasta","_clean.fasta")
list.meta.out<-list.meta.in %>% str_replace_all(meta.fold.in,meta.fold.out) %>% str_replace_all(".tsv","_clean.csv")
list.summary.out<-str_replace_all(list.meta.out, ".csv","_ExclSummary.csv")

#loop through each partition
for (k in 1:n.part){
  
  print(paste("Reading in partition ",k,sep=""))
  ## Read in the GISAID metadata
  meta<-read.table(file = list.meta.in[k], sep = '\t', header = TRUE,fill = TRUE,quote = "")
  
  #what is the most recent data date?
  data.date<-last(sort(meta$Collection.date))
  
  ## Read alignment
  file.seqs<-readDNAStringSet(list.seq.in[k])
  
  #extract filenames
  file.names<-names(file.seqs)
  
  #size of alignment 
  l.og<-length(file.seqs)
  w.og<-width(file.seqs)[1]
  l.can<-length(str_which(file.names, "Canada"))
  
  # Connect sequences to metadata, make a meta column in same format as sequence names
  #note that in this download, meta doesn't contain 'virus.name' but sorting is the same
  nrow(meta)==length(file.names)
  meta$Virus.name<-file.names
  # length(which(file.names %in% meta$Virus.name))==length(file.names) 
  
  # Remove sequences with duplicate names
  # different epi_isl, but no way to tell which is which is seqs
  dups<-which(duplicated(meta$Virus.name))
  l.dups<-length(dups)
  if (l.dups>0){
    meta<-meta[-dups,]
    file.names<-file.names[-dups]
    file.seqs<-file.seqs[-dups]
  }
  length(unique(meta$Virus.name))==length(meta$Virus.name)
  
  #look for duplicates in the meta also
  dup.met<-which(duplicated(meta$Accession.ID))
  l.dup.met<-length(dup.met)
  if (l.dup.met>0){
    meta<-meta[-dup.met,]
  }
  length(unique(meta$Accession.ID))==nrow(meta)
  
  ## Read in the exclude from nextstrain and remove these seqs
  #sequences that are duplicates, too divergent, missing mutations, excessive reversions, bat and pangolin, wrong date
  #this file should be pulled from the ncov github as described in readme.md
  exc<-read.table("scripts/exclude.txt")
  exc<-exc[,1]
  #length(exc)
  
  #make a temporary version of file.names in same format as exclude
  file.names.excl<-as.vector(sapply(meta$Virus.name, function(x)
    paste0(unlist(strsplit( unlist(strsplit(x,split="\\|"))[[1]], "/")) [c(2,3,4)],collapse="/")
  ))
  
  #remove any sequences and names that match exc
  excs<-which(file.names.excl %in% exc)
  l.excs<-length(excs) 
  if (length(excs)>0){
    meta<-meta[-excs,]
  }
  
  ## Remove additional outliers, non-human host, environmental
  #outliers due to sequencing errors or temporal signal discordance identified elsewhere
  #see http://virological.org/t/phylodynamic-analysis-176-genomes-6-mar-2020/356/10
  outliers<- c("EPI_ISL_406592", "EPI_ISL_406595")
  l.out<-length(outliers)
  # outliers<-paste(outliers,collapse="|")
  
  #join these in a running list of exclude that we will exclude based on meta after linking below
  exclude<-outliers
  
  #### NON-HUMAN ####
  # table(meta$Host)
  
  #rows to remove b/c not from humans
  nonhum<-which(meta$Host!="Human")
  l.nonhum<-length(nonhum)
  
  #add to the exclude
  non.hum.id<-meta$Accession.ID[nonhum]
  exclude<-unique(c(exclude,non.hum.id))
  
  #### ENVIRONMENTAL ####
  env<-str_which(meta$Virus.name, "env")
  env2<-which(meta$Host=="Environment")
  env<-unique(c(env,env2))
  l.env<-length(env)
  env.id<-meta$Accession.ID[env]
  exclude<-unique(c(exclude,env.id))
  
  #make a wide single character separated by | for str_which to read properly
  exclude.wide<-paste(exclude,collapse="|")
  
  #remove these all at the same time
  if(length(exclude)>0){
    meta<-meta[-str_which(meta$Accession.ID,exclude.wide),]
  }
  
  ## Rename with no weird chars
  #name as: GISAID_ID/Date ### STRAIN NOT AVALBL anymore
  meta<-meta %>% unite("new.names", c(Accession.ID, Collection.date), sep = "/", remove=FALSE)
  #get rid of any whitespace
  meta$new.names<-str_replace_all(meta$new.names, fixed(" "), "")
  #get rid of any special charactres
  if(length(str_which(meta$new.names,pattern=paste0(c("é","ê","è","ë","ç","//'"),collapse="//|")))>0){
    meta$new.names<-stri_trans_general(str = meta$new.names, 
                                       id = "Latin-ASCII")}
  #change these "-XX" to blanks
  meta$Collection.date<-str_replace_all(meta$Collection.date,"-XX","")
  # new.meta$Collection.date<-str_replace_all(new.meta$Collection.date,"-XX","")
  
  ## ONLY keep incomplete dates from Canada and infer using LSD
  #occurences with incomplete dates...grr
  weird<-c()
  for (i in 1:nrow(meta)){
    if(length(unlist(strsplit(meta$Collection.date[i],split="")))<10) {weird<-c(weird,i); next} #if less than 10 characters (including -), incomplete
  }
  l.weird<-length(weird)
  
  #how many of these are canadian
  l.weird.can<-length(str_which(meta$Location[weird],"Canada"))
  
  #### IN THIS VERSION keep ONLY CANADIAN incomplete dates ####
  #remove Canadian sequences from the weird vector
  weird.remove<-c()
  for (i in 1:nrow(meta)){
    if(length(unlist(strsplit(meta$Collection.date[i],split="")))<10 & 
       !str_detect(meta$Location[i],"Canada")) 
    {weird.remove<-c(weird.remove,i); next}
  }
  
  #check
  length(weird.remove) == l.weird - l.weird.can
  
  # RUN THIS TO REMOVE sequences and meta with incomplete date from non-Canadian 
  if(length(weird.remove)>0){
    meta<-meta[-weird.remove,]
  }
  
  #reduce the canadian incomplete versions to match meta
  still.here<-which(file.names %in% meta$Virus.name)
  all(meta$Virus.name %in% file.names)
  length(still.here)==nrow(meta)
  
  # length(still.here)==nrow(meta) #check
  file.names<-file.names[still.here]
  file.seqs<-file.seqs[still.here] #seqs are ordered the same as names
  length(file.names)==length(file.seqs)
  nrow(meta)==length(file.names)
  
  ####replace the file.names with the new.names####
  #double check order the same
  all(file.names==meta$Virus.name)
  #replace with new names
  file.names<-meta$new.names
  #should be the same
  length(file.names)==length(unique(file.names))
  
  ## Merge file names and sequences back together to export
  l.final<-length(file.names)
  l.final.can<-length(str_which(meta$new.names, "Canada"))
  # names(file.seqs)<-file.names
  
  # replace file names 
  names(file.seqs)<-file.names 
  
  ## Export alignment and metadata
  print("Writing output")
  writeXStringSet(file.seqs,file=list.seq.out[k], format="fasta")
  write.csv(meta,file=list.meta.out[k])
  
  #print summary statistics of the cleaning run to a csv file
  out.colz<-c("n.seq.in","n.seq.can.in","nextstrain.exclude","n.duplicates","n.tempoutliers",
              "n.nonhuman","n.environmental","n.incompletedate","n.incompletedate.can",
              "n.seq.out","n.seq.can.out")
  valuez<-c(l.og, #original length
            l.can, #original number from Canada
            l.excs, #nextstrain exclude
            l.dups, #duplicate names
            l.out, #temporal outliers
            l.nonhum, #non human host
            l.env, #environmental
            length(weird.remove), #incomplete dates, not Canada
            l.weird.can, #incomplete canadian dates
            l.final, #final number of sequences
            l.final.can #final number of Can sequences
  )
  summary.statz<-data.frame(parameter=out.colz,value=valuez)
  write.csv(summary.statz,file=list.summary.out[k],row.names = F)

}
