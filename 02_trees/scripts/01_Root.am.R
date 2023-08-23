#!/usr/bin/env Rscript
# Root tree and remove temporal outliers

#usage:
# Rscript 01_Root.am.R "ft" "EPI_ISL_402125"
# Rscript scripts/01_Root.am.R "202210_analysis/ft_deltagamma" "EPI_ISL_402125"
#setup libraries
library(stringr)
library(ape)
library(tidyverse)
library(gtools)

#### Inputs and outputs here ####
args<-commandArgs(trailingOnly=TRUE)
treebase<-args[1] #folder with trees
root.accession<-args[2] #accession number of the sequence to root with as outgroup

## Or, change inputs manually
# treebase<-"ft2"
# root.accession<-"EPI_ISL_402125"

#setup treefiles
trees.in<-list.files(treebase,pattern = ".tre",full.names = T)
if(length(trees.in)<1){
    print("trees not found");break
}
#sort these in numeric order
trees.in<-mixedsort (trees.in)

#number trees = number of bootstraps
b<-length(trees.in)

#set up outputs
f.out<-paste(treebase,"_root",sep="")
if(!dir.exists(f.out)){dir.create(f.out)}

all.out<-str_replace_all(trees.in, treebase, f.out) #replace the folder
trees.out<-str_replace_all(all.out, ".tre","_root.tre") #replace the name in folder

# On rooting 
# In the manuscript, we followed recommendation (Rambaut 2020 Nature) to root on earliest lineage A viruses:Wuhan/WH04/2020 (EPI_ISL_406801), sampled on 5 January 2020, 
# however here, we root on Wuhan-Hu-1 (GenBank accession no. MN908947) sampled on 26 December 2019, because this fake dataset doesn't have EPI_ISL_406801
#read in trees into a list
trees<-replicate(n=b, vector)
for (i in 1:b){
  print(paste("reading tree",i))
  trees[[i]]<-read.tree(file=trees.in[i])

  #force into binary, ie resolve polytomies randomly
  if (!is.binary(trees[[i]])){
      trees[[i]]<-multi2di(trees[[i]])
  }

  #ROOT

  #earliest lineage B virus: "Wuhan-Hu-1" (Accession number "EPI_ISL402125" or Genbank MN908947) sampled on 2019-12-26
  linB<-root.accession
  linB.tip<- grep(paste("\\b",linB,"\\b",sep=""), trees[[i]]$tip.label)
  if(length(linB.tip)<1) {
      print("accession root not found")
      break
  }
    if(length(linB.tip)>1) {
      #drop all but the first tip
      trees[[i]]<-drop.tip(trees[[i]],linB.tip[2])
      linB.tip<-linB.tip[1]
  }
  
  
  
  #on earliest lineage A virus:Wuhan/WH04/2020 sampled 2020-01-05,
  ##### According to Rambaut, MRCA of whole outbreak as it shares two muts not found in wuhan-hu-1
  # linA<-"EPI_ISL_406801"
  # linA.tip<-str_which(trees[[i]]$tip.label,linA)
  
  #Alternative: find the MRCA node of lineages A and B
  # AB.mrca<-getMRCA(trees[[i]],tip=c(linB.tip,linA.tip))
  
  #re-root it on outgroup
  #trees[[i]]<-ape::unroot(trees[[i]])
  trees[[i]]<-ape::root(trees[[i]],outgroup=linB.tip,resolve.root=TRUE)
  #is.rooted(trees[[i]])
  
  # export trees
  write.tree(trees[[i]],file=trees.out[i])
}

