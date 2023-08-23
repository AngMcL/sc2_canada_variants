#!/usr/bin/env Rscript

# IN the VOC version of this script, infer a molecular clock rate for VOC only (not parental lineages) to exclude VOC seqs that have more than expected diversity
# remove temporal outliers
# export fasta and tree with no outliers (input for iqtree2/lsd2)
# export a dateFile corresponding to tips for lsd2

#usage:
# Rscript scripts/02_RemoveOutliers.Variant.R "ft_alphadelta_root" "~/../../Volumes/PhD_2/2022_SC2_Can/20220322_gisaid/01_subsampled/bootsamples_fasta_202208/masked" "../01_subSample/bootsamples_meta_202208" "alpha"
# Rscript scripts/02_RemoveOutliers.Variant.R "ft_alphadelta_root" "~/../../Volumes/PhD_2/2022_SC2_Can/20220322_gisaid/01_subsampled/bootsamples_fasta_202208/masked" "../01_subSample/bootsamples_meta_202208" "delta"
# Rscript scripts/02_RemoveOutliers.Variant.R "ft_betagammaiota_root" "~/../../Volumes/PhD_2/2022_SC2_Can/20220322_gisaid/01_subsampled/bootsamples_fasta_202208/masked2" "../01_subSample/bootsamples_meta_202208" "gamma"
# Rscript scripts/02_RemoveOutliers.Variant.R "ft_betagammaiota_root" "~/../../Volumes/PhD_2/2022_SC2_Can/20220322_gisaid/01_subsampled/bootsamples_fasta_202208/masked2" "../01_subSample/bootsamples_meta_202208" "beta"
# Rscript scripts/02_RemoveOutliers.Variant.R "ft_betagammaiota_root" "~/../../Volumes/PhD_2/2022_SC2_Can/20220322_gisaid/01_subsampled/bootsamples_fasta_202208/masked2" "../01_subSample/bootsamples_meta_202208" "iota"
# Rscript scripts/02_RemoveOutliers.Variant.R "ft_etaepsilon_root" "~/../../Volumes/PhD_2/2022_SC2_Can/20220322_gisaid/01_subsampled/bootsamples_fasta_202208/masked2" "../01_subSample/bootsamples_meta_202208" "eta"
# Rscript scripts/02_RemoveOutliers.Variant.R "ft_etaepsilon_root" "~/../../Volumes/PhD_2/2022_SC2_Can/20220322_gisaid/01_subsampled/bootsamples_fasta_202208/masked2" "../01_subSample/bootsamples_meta_202208" "epsilon"
# Rscript scripts/02_RemoveOutliers.Variant.R "ft_kappamu_root" "~/../../Volumes/PhD_2/2022_SC2_Can/20220322_gisaid/01_subsampled/bootsamples_fasta_202208/masked" "../01_subSample/bootsamples_meta_202208" "kappa"
# Rscript scripts/02_RemoveOutliers.Variant.R "ft_kappamu_root" "~/../../Volumes/PhD_2/2022_SC2_Can/20220322_gisaid/01_subsampled/bootsamples_fasta_202208/masked" "../01_subSample/bootsamples_meta_202208" "mu"
# Rscript scripts/02_RemoveOutliers.Variant.R "202210_analysis/ft_other_root" "~/../../Volumes/PhD_2/2022_SC2_Can/20220322_gisaid/01_subsampled/bootsamples_fasta_202210/diffs_202210/diffs_masked" "../01_subSample/bootsamples_meta_202210" "other"
# Rscript scripts/02_RemoveOutliers.Variant.R "BA_analysis/ft_BA11_root" "~/../../Volumes/PhD_2/2022_SC2_Can/20220322_gisaid/01_subsampled/bootsamples_fasta_omi/BA_masked" "../01_subSample/bootsamples_meta_omi" "BA.1.1"
# Rscript scripts/02_RemoveOutliers.Variant.R "BA_analysis/ft_BA11_unif_root" "~/../../Volumes/PhD_2/2022_SC2_Can/20220322_gisaid/01_subsampled/bootsamples_fasta_omi_unif/BA_unif_masked" "../01_subSample/bootsamples_meta_unif" "BA.1.1"
# Rscript scripts/02_RemoveOutliers.Variant.R "BA_analysis/ft_BA1_root" "~/../../Volumes/PhD_2/2022_SC2_Can/20220322_gisaid/01_subsampled/bootsamples_fasta_omi/BA_masked" "../01_subSample/bootsamples_meta_omi" "BA.1"
# Rscript scripts/02_RemoveOutliers.Variant.R "BA_analysis/ft_BA1_unif_root" "~/../../Volumes/PhD_2/2022_SC2_Can/20220322_gisaid/01_subsampled/bootsamples_fasta_omi_unif/BA_unif_masked" "../01_subSample/bootsamples_meta_unif" "BA.1"

# setup libraries
library(stringr)
library(ape)
library(treeio)
library(adephylo)
library(tidyverse)
library(gtools)
library(lubridate)

#### Inputs and outputs ####
args<-commandArgs(trailingOnly=TRUE)
treebase<-args[1] #folder with rooted trees
fastabase<-args[2] #folder with alignments corresponding to trees
metabase<-args[3] #folder with subsampled meta
focusVOC<-args[4]

# ## Or, change inputs manually
treebase<-"BA_analysis/ft_BA2_root"
fastabase<-"~/../../Volumes/PhD_2/2022_SC2_Can/20220322_gisaid/01_subsampled/bootsamples_fasta_omi/BA_masked"
metabase<-"../01_subSample/bootsamples_meta_omi"
focusVOC<- "BA.2"

# treebase<-"BA_analysis/ft_BA2_unif_root"
# fastabase<-"~/../../Volumes/PhD_2/2022_SC2_Can/20220322_gisaid/01_subsampled/bootsamples_fasta_omi_unif/BA_unif_masked"
# metabase<-"../01_subSample/bootsamples_meta_unif"
# focusVOC<- "BA.2"

#setup treefiles, alignment files using folders
trees.in<-list.files(treebase,pattern = ".tre",full.names = T) %>% mixedsort ()
if(length(trees.in)<1 | all(!file.exists(trees.in))){
    print("trees not found");break
}

#restrict to focus VOC
trees.in<-str_subset(trees.in,paste0("/",focusVOC,"_"))

#number of bootstraps/trees
b<-length(trees.in)

#setup tempest data
tempest.in<-str_replace_all(trees.in, ".tre",".tsv")
if(length(tempest.in)<1 | all(!file.exists(tempest.in))){
    print("tempest tsv not found");break
}

#alignments in (to remove outliers from) - only those matching trees (in batches)
aligns.in<-str_replace_all(trees.in, treebase,fastabase)
aligns.in<-str_replace_all(aligns.in, "_root.tre",".fasta")

if(length(aligns.in)<1 | all(!file.exists(aligns.in))){ 
    print("fasta not found");break
}

#metadata in
meta.in<-str_replace_all(trees.in, treebase,metabase)
meta.in<-str_replace_all(meta.in, "_mask_root.tre",".csv")
if(length(meta.in)<1 | all(!file.exists(meta.in))){ 
  print("meta not found");break
}

#set up outputs
f.out<-paste(treebase,"_res",sep="")
if(!dir.exists(f.out)){dir.create(f.out)}

all.out<-str_replace_all(trees.in, treebase, f.out)
trees.out<-str_replace_all(all.out, ".tre","_res.tre")
aligns.out<-str_replace_all(all.out, ".tre","_res.fasta")
dates.out<-str_replace_all(all.out, ".tre","_res_dates.txt")

#figures out
fold.out<-paste(f.out,"results",sep="_")
if(!dir.exists(fold.out)){dir.create(fold.out)}

#summarytable out
summ.table.out<-paste(fold.out, "/",focusVOC,"_ExclSummary.csv",sep="")
summ.table<-data.frame(trees=trees.in)

#### READ IN DATA #####
print("reading in data")

## Read in the trees into a list
trees<-replicate(n=b, vector)
summ.table$ntip.start<-NA

for (i in 1:b){
  trees[[i]]<-read.tree(file=trees.in[i])
  summ.table$ntip.start[i]<-Ntip(trees[[i]])
  
}

## Read in tempest data with residuals here
resid<-replicate(n=b,vector())
for (i in 1:b){
  resid[[i]]<-read.table(file=tempest.in[i] ,sep="\t",header=T)
  #boot name
  boot.nm<-str_replace_all(last(unlist(strsplit(tempest.in[i],"/"))),"_mask_root.tsv","")
  resid[[i]]$boot<-boot.nm
}

## Read in the metadata and extract the complete dates into resid object
meta<-replicate(n=b,vector())
for (i in 1:b){
  meta[[i]]<-read.table(file=meta.in[i] ,sep=",",header=T)
  meta[[i]]$boot<-i
  meta[[i]]$tip<-meta[[i]]$new.names
  #calculate meta decdates
  meta[[i]]$decdate<-decimal_date(as.Date(meta[[i]]$date))
  #change this so we can use later
  meta[[i]]$incomp[meta[[i]]$country!="Canada"]<-0
  #trim the object down
  meta[[i]]<-meta[[i]][,c("tip","var.WHO","Lineage","decdate","incomp","date")]
  colnames(meta[[i]])[ncol(meta[[i]])]<-"fulldate"
}

#check matchy
# length(which(!meta[[i]]$tip %in% resid[[i]]$tip))

#Link to resid by tiplabel
for (i in 1:b){
  resid[[i]]<-left_join(resid[[i]], meta[[i]],by="tip")
  #replace the resid date with the meta decdates
  resid[[i]]$date<-resid[[i]]$decdate
  resid[[i]]<-resid[[i]][,-which(colnames(resid[[i]])=="decdate")]
}


## No longer needed
##REMOVE residual datapoint if date=0 (incomplete dates)
#can't use for analysis, so will carry forward blindly... 
# resid[[i]]$tip[resid[[i]]$date==0] #all the months with no dates...thanks ontario
# for (i in 1:b){
#   if(any(resid[[i]]$date==0)){
#     resid[[i]]<-resid[[i]][-which(resid[[i]]$date==0),]
#   }
# }

#restrict to variants - NEW
#note this also removes Wuhan-hu-1 from exclusion
if(focusVOC != "other"){ #as long as not analyzing Other specifically:
  for (i in 1:b){
    if(any(is.na(resid[[i]]$var.WHO))){ #NA
      resid[[i]]<-resid[[i]][-which(is.na(resid[[i]]$var.WHO)),]
    }
    if(any(resid[[i]]$var.WHO=="Other")){ 
      resid[[i]]<-resid[[i]][-which(resid[[i]]$var.WHO=="Other"),]
    }
  }
}

#for omicron, stratify by subvar for the fit
###make a separate script for this, pulling from the "...VOC.omi.R" script

#### Make a linear fit and calculate residuals manually #####
print("calculating residuals for variant-specific lm fits")
for (i in 1:b){
  m<-lm(resid[[i]]$distance ~ resid[[i]]$date)
  # plot(resid[[i]]$distance ~ resid[[i]]$date)
  #calculate the residuals and replace in the df
  resid[[i]]$residual<-residuals(m)
}

#### PLOT RESIDUALS ####
print("plot rtt and residuals")
##Look at distribution of residuals and clock in each
##  plot residuals by BOOTS
plot.residuals<-function(res){
  #root-to-tip regression
  lm<-summary(lm(res$distance ~ res$date))
  slope<-formatC(coef(lm)[[2]],format="e", digits=2)
  slope.n<-coef(lm)[[2]] #numeric
  int<-coef(lm)[[1]]
  eq<-paste("y = ",slope," x + ",signif(int,digits=3),sep="")
  R2<-signif(lm$adj.r.squared,digits = 3)
  
  #set positions
  max.x<-max(res$date)
  max.y<-max(res$distance)
  
  BOOT<-res$boot[1]
  
  #overlay line on scatter plot, colored by residual
  ggplot(res)+
    geom_point(aes(x=date,y=distance, color=abs(residual)))+
    scale_color_continuous(type = "viridis")+
    labs(x="Year",y="Root-to-tip distance",color="Abs.Value\nResidual",title=BOOT)+
    geom_abline(slope=slope.n,intercept = int,color="black")+
    annotate("text",label=paste("Adj.R^2=",R2,"\n",eq,sep=""),x=max.x,y=max.y*0.9,hjust=1,color="black",size=4)
  ggsave(paste(fold.out,"/TemporalSignalResidInclOutliers_",BOOT,".png",sep=""),height=3,width=4,units="in")
}

#run it over the list and save plots into results
lapply(resid,plot.residuals)

#pull model fits into summary table
summ.table$modelfit.start<-NA
summ.table$R2.start<-NA
#redundant, but easier
for (i in 1:b){
  lm<-summary(lm(data=resid[[i]], distance ~ date))
  slope<-formatC(coef(lm)[[2]],format="e", digits=2)
  slope.n<-coef(lm)[[2]] #numeric
  int<-coef(lm)[[1]]
  eq<-paste("y = ",slope," x + ",signif(int,digits=3),sep="")
  R2<-signif(lm$adj.r.squared,digits = 3)
  summ.table$modelfit.start[i]<-eq
  summ.table$R2.start[i]<-R2
}

##Look at distribution of residuals in each
##  plot density residuals with cutoff}

cutoff<-0.001

plot.noabs.residuals<-function(res){

  BOOT<-res$boot[1]
  
  #calculate mean + 4 stdev for each for residual
  # cutoff<-mean((res$residual))+(4*sd((res$residual))) #too variable
  
  #calculate n cutoff (on negative or positive side)
  n.cutoff<-length(which(abs(res$residual)>cutoff))
  
  max.x<-max((res$residual))
  max.y<-.95
  
  #density plot, colored by residual
  ggplot(res)+
    geom_density(aes(x=(residual), y=..scaled..))+
    labs(x="Residual",y="Scaled density",title=BOOT)+
    geom_vline(xintercept=cutoff,color="red",linetype=2)+
    geom_vline(xintercept=cutoff,color="red",linetype=2)+
    annotate("text",label=paste("cutoff = ",signif(cutoff,digits=3),"\n", "n excluded = ",n.cutoff,sep=""),
             x=max.x,y=max.y*0.8,hjust=1,color="black",size=3)
  ggsave(paste(fold.out,"/densityResiduals_N-outliers_",BOOT,".png",sep=""),height=3,width=4,units="in")
}

lapply(resid,plot.noabs.residuals)

#Extract this into summary table
summ.table$BOOT<-NA
summ.table$cutoff<-NA
summ.table$Ntip.residual.out<-NA

## Remove residuals past a cutoff
#apply cutoff as above
for (i in 1:b){
  outliers<-which(abs(resid[[i]]$residual)>cutoff ) #abs() here cuts off too low and too high
  ##Extract this into summary table
  summ.table$BOOT[i]<-resid[[i]]$boot[1]
  summ.table$cutoff[i]<-cutoff
  summ.table$Ntip.residual.out[i]<-length(outliers)  
  
  #remove the outliers
  if(length(outliers)>0){
    trees[[i]]<-ape::drop.tip(phy=trees[[i]],tip=resid[[i]]$tip[outliers])
    resid[[i]]<-resid[[i]][-outliers,]
  }
}

#check
# Ntip(trees[[2]])==nrow(resid[[2]])
#not the same if excluded incomplete dates (0s) from resid

# note that the fake data was simulated under a strict clock so shouldn't really have temp outliers, unless made more stringent

# plot the residuals again
plot.new.residuals<-function(res){
  
  # res<-res[which(res$var.WHO!=""),]
  anyna<-which(is.na(res$date))
  if(length(anyna)>0){
   res<-res[-anyna,]
   } 

  #root-to-tip regression
  lm<-summary(lm(res$distance ~ res$date))
  res$residual<-residuals(lm(res$distance ~ res$date))
  slope<-formatC(coef(lm)[[2]],format="e", digits=2)
  slope.n<-coef(lm)[[2]] #numeric
  int<-coef(lm)[[1]]
  eq<-paste("y = ",slope," x + ",signif(int,digits=3),sep="")
  R2<-signif(lm$adj.r.squared,digits = 3)
  
  #set positions
  max.x<-max(res$date)
  max.y<-max(res$distance)
  min.y<-min(res$distance)
  
  BOOT<-res$boot[1]
  print(BOOT)
  #overlay line on scatter plot, colored by residual
  ggplot(res)+
    geom_point(aes(x=date,y=distance, color=abs(residual)))+
    scale_color_continuous(type = "viridis")+
    labs(x="Year",y="Root-to-tip distance",color="Abs.Value\nResidual",title=BOOT)+
    geom_abline(slope=slope.n,intercept = int,color="black")+
    annotate("text",label=paste("Adj.R^2=",R2,"\n",eq,sep=""),x=max.x,y=min.y,vjust=0,hjust=1,color="black",size=3)
  ggsave(paste(fold.out,"/TemporalSignalResidExclOutliers_",BOOT,".png",sep=""),height=3,width=4,units="in")
  
}

#this gets overwritten after exlcuding pendant edges
# lapply(resid,plot.new.residuals)


#pull model fits into summary table
summ.table$modelfit.excl.resid.out<-NA
summ.table$R2.excl.resid.out<-NA
#redundant, but easier
for (i in 1:b){
  if(any(is.na(resid[[i]]$date))){
   print("na dates detected")
   wna<-which(is.na(resid[[i]]$date))
   print(length(wna))
   resid[[i]]<-resid[[i]][-wna,]
  }  
  
  lm<-summary(lm(data=resid[[i]][resid[[i]]$var.WHO!="",], distance ~ date))
  slope<-formatC(coef(lm)[[2]],format="e", digits=2)
  slope.n<-coef(lm)[[2]] #numeric
  int<-coef(lm)[[1]]
  eq<-paste("y = ",slope," x + ",signif(int,digits=3),sep="")
  R2<-signif(lm$adj.r.squared,digits = 3)
  summ.table$modelfit.excl.resid.out[i]<-eq
  summ.table$R2.excl.resid.out[i]<-R2
}


#Exclude outliers based on long terminal branch length
print("pendant edge function")
## calculate the pendant edge length (last branch leading to tip) for all tips 
#PE Function (from Jeffrey Joy)
pendant.edge <-function(tree, scale=F){
  
  if(is.rooted(tree)==FALSE)
    warning("A rooted phylogeny is required for meaningful output of this function", call.=FALSE)
  if(scale==TRUE){
    
    #Scale tree to have unit depth (for an ultrametric tree) or scale all branches to unit length (for an additive tree)
    if(is.ultrametric(tree)==TRUE)
      tree<- rescaleTree(tree, 1) else
        tree$edge.length<- tree$edge.length/sum(tree$edge.length)
  }
  
  edge.index<-which(tree$edge[,2] %in% 1:length(tree$tip.label))
  w<- as.data.frame(tree$edge.length[edge.index])
  results<- cbind(tree$tip.label, w)
  
  names(results)<- c("Species", "w")
  results
}

## See if long branches still an issue
##  pendant edge with cutoffs}
#max number of muts expected on one branch, increased to 30 from 12 (prev analysis)
max.muts<-30
pend.cutoff<-max.muts/29903
summ.table$max.muts<-max.muts
summ.table$Ntip.pendant.out<-NA

for (i in 1:b){
  pend<-pendant.edge(trees[[i]])
  
  n.cutoff<-length(which(pend$w>pend.cutoff))
  print(n.cutoff)
  summ.table$Ntip.pendant.out[i]<-n.cutoff
  
  max.x<-max(pend$w)
  max.y<-.95
  BOOT<-resid[[i]]$boot[1]
  
  #density plot of # mutations
  ggplot(pend)+
    geom_density(aes(x=round(w*29903), y=..scaled..))+
    labs(x="Number of mutations on pendant edge",y="Scaled density")+
    geom_vline(xintercept=max.muts,color="red",linetype=2)+
    annotate("text",label=paste("cutoff = ",max.muts,"\n", "n excluded = ",n.cutoff,sep=""),
             x=30,y=0.95,hjust=1,color="black",size=3)
  ggsave(paste(fold.out,"/PendantEdges_N-outliers_",BOOT,".png",sep=""),height=3,width=4,units="in")
}


## Exclude long pendant edges (more than 12 muts)
print("excluding outliers based on pendant edges")
#apply cutoff
for (i in 1:b){
  pend<-pendant.edge(trees[[i]])
  
  outliers<-which(pend$w > pend.cutoff) 
  if(length(outliers)>0){
    trees[[i]]<-ape::drop.tip(phy=trees[[i]],tip=pend$Species[outliers])
    pend2<-pend[-outliers,]
    resid[[i]]<-resid[[i]][-which(resid[[i]]$tip %in% pend$Species[outliers] ),]
  }
}


#check for one boot
# Ntip(trees[[length(trees)]])==nrow(pend2)
# Ntip(trees[[i]])

## Re-plot and save the residuals, excluding seqs with high residuals and long pendant edges (final strict clock)
print("re-plotting residuals")
lapply(resid,plot.new.residuals)

#pull model fits into summary table
print("summarizing in table")
summ.table$modelfit.end<-NA
summ.table$R2.end<-NA
#redundant, but easier
for (i in 1:b){
  lm<-summary(lm(data=resid[[i]][resid[[i]]$var.WHO!="",], distance ~ date))
  print("making lm")  
  slope<-formatC(coef(lm)[[2]],format="e", digits=2)
  slope.n<-coef(lm)[[2]] #numeric
  int<-coef(lm)[[1]]
  eq<-paste("y = ",slope," x + ",signif(int,digits=3),sep="")
  R2<-signif(lm$adj.r.squared,digits = 3)
  summ.table$modelfit.end[i]<-eq
  summ.table$R2.end[i]<-R2
}

#EXPORT SUMMARY TABLE
write.csv(summ.table, summ.table.out, row.names=F)


## exclude these outliers from the alignment
print("excluding outliers from alignments")

align<-replicate(n=b,vector)
for (i in 1:b){
  fast<-read.FASTA(aligns.in[i]) 
  
  #check, should be true
  length(which(names(fast) %in% trees[[i]]$tip.label))==length(trees[[i]]$tip.label)
  # names(fast)[which(!names(fast) %in% trees[[i]]$tip.label)]
  # trees[[i]]$tip.label [which(!trees[[i]]$tip.label %in% names(fast))]
  
  #now can match on name
  #only keep the keepers
  keep<-which(names(fast)%in% trees[[i]]$tip.label)
  align[[i]]<-fast[keep]
}

# test, these should be the same
for (i in 1:b){
  print(length(align[[i]])==Ntip(trees[[i]]))
}

## make sure that naming order is the same in fasta and tree (for lsd2)
for (i in 1:b){
  if(!all(names(align[[i]])==trees[[i]]$tip.label))
    #reorder the fasta to match tree
    align[[i]]<-align[[i]][order(match(names(align[[i]]), trees[[i]]$tip.label))]
    print(all(names(align[[i]])==trees[[i]]$tip.label)) #should all be T after
}

## Extract the dates for LSDin datefile
print("extracting dates")
## export a datefile from this tree
dates.list<-replicate(b,vector())
for (i in 1:b){
  #extract tipnames
  nm<-trees[[i]]$tip.label
  #get dates from meta/resid
  dates.list[[i]]<-data.frame(node_name=nm,
                              date=resid[[i]]$fulldate[match(nm, resid[[i]]$tip)])
  #if any NA, take from the tipname
  nay<-which(is.na(dates.list[[i]]$date))
  if(length(nay)>0){
      for (j in 1:length(nay)){
        #strsplit it from tip
        name<-dates.list[[i]]$node_name[nay[j]]
        dates.list[[i]]$date[nay[j]]<-last(unlist(strsplit(name,split="/")))
      }
    }
}
  
# head(dates.list[[1]])

##export the trees, alignments, datefiles
print("writing the final objects")
for (i in 1:b){
  write.tree(trees[[i]], trees.out[i])
  write.FASTA(align[[i]],aligns.out[i])
  write.table(dates.list[[i]], dates.out[i], row.names=F, sep="\t", quote = F, col.names=F)
}
