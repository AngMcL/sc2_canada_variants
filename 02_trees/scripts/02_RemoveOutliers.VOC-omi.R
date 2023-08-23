#!/usr/bin/env Rscript

# IN the VOC version of this script, infer a molecular clock rate for VOC only (not parental lineages) to exclude VOC seqs that have more than expected diversity
# remove temporal outliers
# export fasta and tree with no outliers (input for iqtree2/lsd2)
# export a dateFile corresponding to tips for lsd2

#usage:
# $Rscript RemoveOutliers.am.R "ft_root" "bootsamples"
# $Rscript scripts/RemoveOutliers.VOC.R "ft-d1_root" "~/../../Volumes/PhD_2/2022_SC2_Can/20220322_gisaid/01_subsampled/bootsamples_fasta_50k/masked/masked_align" "../01_subSample/bootsamples_meta"
# $Rscript scripts/RemoveOutliers.VOC.R "ft-d2_root" "~/../../Volumes/PhD_2/2022_SC2_Can/20220322_gisaid/01_subsampled/bootsamples_fasta_50k/masked/masked_align" "../01_subSample/bootsamples_meta"
# $Rscript scripts/RemoveOutliers.VOC.R "ft_2_root" "~/../../Volumes/PhD_2/2022_SC2_Can/20220322_gisaid/01_subsampled/bootsamples_fasta_50k/masked/masked_align" "../01_subSample/bootsamples_meta"

# setup libraries
library(stringi)
library(stringr)
library(ape)
library(treeio)
library(adephylo)
library(tidyverse)
library(gtools)

#### Inputs and outputs ####
args<-commandArgs(trailingOnly=TRUE)
treebase<-args[1] #folder with rooted trees
fastabase<-args[2] #folder with alignments corresponding to trees
metabase<-args[3]

# ## Or, change inputs manually
treebase<-"ft2_root"#"ft1_root"
fastabase<-"~/../../Volumes/PhD_2/2022_SC2_Can/20220322_gisaid/01_subsampled/bootsamples_fasta/masked/masked_align"
metabase<-"../01_subSample/bootsamples_meta"
focusVOC<-"omicron"

#setup treefiles, alignment files using folders
trees.in<-list.files(treebase,pattern = ".tre",full.names = T) %>% mixedsort ()
if(length(trees.in)<1 | all(!file.exists(trees.in))){
    print("trees not found");break
}

#restrict to delta only (alpha, beta fine)
trees.in<-str_subset(trees.in,focusVOC)

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

## Read in the metadata
meta<-replicate(n=b,vector())
for (i in 1:b){
  meta[[i]]<-read.table(file=meta.in[i] ,sep=",",header=T)
  meta[[i]]$boot<-i
  meta[[i]]$tip<-meta[[i]]$new.names
  meta[[i]]<-meta[[i]][,c("tip","Variant","Lineage")]
}

#check matchy
# length(which(!meta[[i]]$tip %in% resid[[i]]$tip))

#Link to resid by tiplabel
for (i in 1:b){
  resid[[i]]<-left_join(resid[[i]], meta[[i]],by="tip")
}

##REMOVE residual datapoint if date=0 (incomplete dates)
#can't use for analysis, so will carry forward blindly... 
# resid[[i]]$tip[resid[[i]]$date==0] #all the months with no dates...thanks ontario
for (i in 1:b){
  if(any(resid[[i]]$date==0)){
    resid[[i]]<-resid[[i]][-which(resid[[i]]$date==0),]
  }
}


#restrict to variants - NEW
#note this also removes Wuhan-hu-1 from exclusion
for (i in 1:b){
  if(any(is.na(resid[[i]]$Variant))){ #NA
    resid[[i]]<-resid[[i]][-which(is.na(resid[[i]]$Variant)),]
  }
  if(any(resid[[i]]$Variant=="")){ #blanks
    resid[[i]]<-resid[[i]][-which(resid[[i]]$Variant==""),]
  }
}

#restrict to specific VOC
#for delta, remove B.1.617.1 (Kappa) from Delta residuals
if(focusVOC=="delta"){
  for (i in 1:b){
    resid[[i]]<-resid[[i]][-str_which(resid[[i]]$Variant,"Kappa"),]
  }
}
#for gamma, remove P.2 (Zeta), P.3 (Theta)
if(focusVOC=="gamma"){
  for (i in 1:b){
    resid[[i]]<-resid[[i]][-str_which(resid[[i]]$Variant,"Theta|Zeta"),]
  }
}

#for omicron, stratify by subvar for the fit
#NOTE will have to make this more complicated for futre builds with BA.4, BA.5 etc
library(stringi)
if(focusVOC=="omicron"){
  for (i in 1:b){
    resid[[i]]$subvar<-NA
    resid[[i]]$subvar[str_which(resid[[i]]$Lineage, "BA.1")]<-"BA.1" #will also include desc
    resid[[i]]$subvar[str_which(resid[[i]]$Lineage, "BA.2")]<-"BA.2"
    resid[[i]]$subvar[str_which(resid[[i]]$Lineage, "BA.3")]<-"BA.3"
  }
}

#need to change the strucutre of the summary for different fits for each subvar
subvars<-sort(unique(resid[[i]]$subvar))
l.sub<-length(subvars)
summ.table<-data.frame(trees=rep(trees.in,length(unique(resid[[i]]$subvar))), 
                       subvar=c(rep(subvars[1],length(trees.in)),
                                rep(subvars[2],length(trees.in)),
                                rep(subvars[3],length(trees.in))), 
                       #might need to add more here if more subvars
                       ntip.start=NA) #for each subvar in the meta
for (i in 1:b){ # for all boots
  for (j in 1:l.sub){ #and all subvars
    mysub<-subvars[j]
    n.sub<-length(str_which(meta[[i]]$Lineage,mysub))
    summ.table$ntip.start[which(summ.table$trees==trees.in[i] & summ.table$subvar==mysub)]<-n.sub
  }
}

#### Make a linear fit and calculate residuals manually FOR EACH SUBVAR#####
for (i in 1:b){
  for (j in 1:l.sub){
    index.var<-which(resid[[i]]$subvar==subvars[j])
    sub.df<-resid[[i]] [index.var,]
    m<-lm(sub.df$distance ~ sub.df$date)
    # plot(resid[[i]]$distance ~ resid[[i]]$date)
    #calculate the residuals and replace in the df
    resid[[i]]$residual[index.var]<-residuals(m) 
  }
}

#### PLOT RESIDUALS IN EACH SUBVAR ####
print("remove outliers based on residuals")

##Look at distribution of residuals and clock in each
#prepare to pull model fits into summary table
summ.table$modelfit.start<-NA
summ.table$R2.start<-NA

##  plot residuals by BOOTS and SUBvAR
for (i in 1:b){
  for (j in 1:l.sub){
    mysub<-subvars[j]
    index.var<-which(resid[[i]]$subvar==mysub)
    res<-resid[[i]][index.var,]
    BOOT<-paste(res$boot[1],mysub)
    
    #root-to-tip regression
    lm<-summary(lm(res$distance ~ res$date))
    slope<-formatC(coef(lm)[[2]],format="e", digits=2)
    slope.n<-coef(lm)[[2]] #numeric
    int<-coef(lm)[[1]]
    eq<-paste("y = ",slope," x + ",signif(int,digits=3),sep="")
    R2<-signif(lm$adj.r.squared,digits = 3)
    
    #add these to the table
    index.tree.var<-which(summ.table$trees==trees.in[i] & summ.table$subvar==mysub)
    summ.table$modelfit.start[index.tree.var]<-eq
    summ.table$R2.start[index.tree.var]<-R2    
    
    #set positions
    max.x<-max(res$date)
    max.y<-max(res$distance)
    
    #overlay line on scatter plot, colored by residual
    ggplot(res)+
      geom_point(aes(x=date,y=distance, color=abs(residual)))+
      scale_color_continuous(type = "viridis")+
      labs(x="Year",y="Root-to-tip distance",color="Abs.Value\nResidual",title=BOOT)+
      theme_bw()+
      geom_abline(slope=slope.n,intercept = int,color="black")+
      annotate("text",label=paste("Adj.R^2=",R2,"\n",eq,sep=""),
               x=max.x,y=max.y*0.9,hjust=1,color="black",size=4)
    ggsave(paste(fold.out,"/TemporalSignalResidInclOutliers_",BOOT,".png",sep=""),
           height=3,width=4,units="in")
  }
}

##Look at distribution of residuals in each

#Prepare to extract this into summary table
summ.table$BOOT<-NA
summ.table$cutoff<-NA
summ.table$Ntip.residual.out<-NA

##  plot density residuals with cutoff AND 
## Remove residuals past a cutoff
cutoff<-0.001

for (i in 1:b){
  for (j in 1:l.sub){
    mysub<-subvars[j]
    index.var<-which(resid[[i]]$subvar==mysub)
    res<-resid[[i]][index.var,]
    BOOT<-paste(res$boot[1],mysub)
    
    #calculate mean + 3 stdev for each for residual
  
    #calculate n cutoff (on negative or positive side)
    n.cutoff<-length(which(abs(res$residual)>cutoff))
    max.x<-max((res$residual))
    max.y<-.95
    
    #density plot, colored by residual
    ggplot(res)+
      geom_density(aes(x=(residual), y=..scaled..))+
      labs(x="Residual",y="Scaled density",title=BOOT)+
      theme_bw()+
      geom_vline(xintercept=cutoff,color="red",linetype=2)+
      geom_vline(xintercept=cutoff,color="red",linetype=2)+
      geom_vline(xintercept=-cutoff,color="red",linetype=2)+
      annotate("text",label=paste("cutoff = ",signif(cutoff,digits=3),"\n", "n excluded = ",n.cutoff,sep=""),
               x=max.x,y=max.y*0.8,hjust=1,color="black",size=3)
    ggsave(paste(fold.out,"/densityResiduals_N-outliers_",BOOT,".png",sep=""),height=3,width=4,units="in")
    
    #which ones are outliers?
    
    outliers<-which(abs(res$residual)>cutoff ) #abs() here cuts off too low and too high
    ##Extract this into summary table
    index.tree.var<-which(summ.table$trees==trees.in[i] & summ.table$subvar==mysub)
    summ.table$BOOT[index.tree.var]<-BOOT
    summ.table$cutoff[index.tree.var]<-cutoff
    summ.table$Ntip.residual.out[index.tree.var]<-length(outliers)  
    
    #remove the outliers
    if(length(outliers)>0){
      outlier.tips<-res$tip[outliers]
      trees[[i]]<-ape::drop.tip(phy=trees[[i]],tip=outlier.tips)
      resid[[i]]<-resid[[i]][-which(resid[[i]]$tip %in% outlier.tips),]
    }
  } #end of loop over subvars
} #end of boot loop


#check
# Ntip(trees[[1]])==nrow(resid[[1]])
#not the same if excluded incomplete dates (0s) from resid

#pull model fits excluding residuals into summary table
#dont really need to plot these
summ.table$modelfit.excl.resid.out<-NA
summ.table$R2.excl.resid.out<-NA

for (i in 1:b){
  for (j in 1:l.sub){
    mysub<-subvars[j]
    index.var<-which(resid[[i]]$subvar==mysub)
    res<-resid[[i]][index.var,]
      
    lm<-summary(lm(data=res, distance ~ date))
    slope<-formatC(coef(lm)[[2]],format="e", digits=2)
    slope.n<-coef(lm)[[2]] #numeric
    int<-coef(lm)[[1]]
    eq<-paste("y = ",slope," x + ",signif(int,digits=3),sep="")
    R2<-signif(lm$adj.r.squared,digits = 3)
    #add to the table
    index.tree.var<-which(summ.table$trees==trees.in[i] & summ.table$subvar==mysub)
    summ.table$modelfit.excl.resid.out[index.tree.var]<-eq
    summ.table$R2.excl.resid.out[index.tree.var]<-R2
  }
}


#Exclude outliers based on long terminal branch length
print("pendant edge exclusion")
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

##  calculate pendant edges and apply cutoff
#max number of muts expected on one branch, increased to 30 from 12 (prev analysis)
max.muts<-30
pend.cutoff<-max.muts/29903
summ.table$max.muts<-max.muts
summ.table$Ntip.pendant.out<-NA

#####
### COME BACK TO HERE ####
#####

print("excluding outliers based on pendant edges")
for (i in 1:b){
  pendants<-pendant.edge(trees[[i]])
  
  for (j in 1:l.sub){
    mysub<-subvars[j]
    index.var<-which(resid[[i]]$subvar==mysub)
    res<-resid[[i]][index.var,]
    BOOT<-paste(res$boot[1],mysub)

    pend.var<-pendants[match(res$tip, pendants$Species),]
    outliers<-which(pend.var$w > pend.cutoff) 
    n.cutoff<-length(outliers)
    
    #add to the summary table
    index.tree.var<-which(summ.table$trees==trees.in[i] & summ.table$subvar==mysub)
    summ.table$Ntip.pendant.out[index.tree.var]<-n.cutoff
    
    max.x<-max(pend.var$w)
    max.y<-.95

    #density plot of # mutations on pend edges for this subvar
    ggplot(pend.var)+
      geom_density(aes(x=round(w*29903), y=..scaled..))+
      labs(x="Number of mutations on pendant edge",y="Scaled density")+
      theme_bw()+
      geom_vline(xintercept=max.muts,color="red",linetype=2)+
      annotate("text",label=paste("cutoff = ",max.muts,"\n", "n excluded = ",n.cutoff,sep=""),
               x=30,y=0.95,hjust=1,color="black",size=3)
    ggsave(paste(fold.out,"/PendantEdges_N-outliers_",BOOT,".png",sep=""),height=3,width=4,units="in")

    #remove outliers
    if(n.cutoff>0){
      trees[[i]]<-ape::drop.tip(phy=trees[[i]],tip=pend.var$Species[outliers])
      #if any of these tips are in the residual df...remove them 
      if(length(which(res$tip %in% pend.var$Species[outliers] ))>0){
        outlier.tips<-res$tip [which(res$tip %in% pend.var$Species[outliers] )]
        resid[[i]]<-resid[[i]][-which(resid[[i]]$tip %in% outlier.tips),]
      }
    }
    
    
    } # end of loop across subvars
}


# plot the residuals and temporal linear fit excluding outliers


#pull model fits into summary table
summ.table$modelfit.end<-NA
summ.table$R2.end<-NA

for (i in 1:b){
  for (j in 1:l.sub){
    mysub<-subvars[j]
    index.var<-which(resid[[i]]$subvar==mysub)
    res<-resid[[i]][index.var,]
    BOOT<-paste(res$boot[1],mysub)
    
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

    #overlay line on scatter plot, colored by residual
    ggplot(res)+
      geom_point(aes(x=date,y=distance, color=abs(residual)))+
      scale_color_continuous(type = "viridis")+
      labs(x="Year",y="Root-to-tip distance",color="Abs.Value\nResidual",title=BOOT)+
      theme_bw()+
      geom_abline(slope=slope.n,intercept = int,color="black")+
      annotate("text",label=paste("Adj.R^2=",R2,"\n",eq,sep=""),x=max.x,y=min.y,vjust=0,hjust=1,color="black",size=3)
    ggsave(paste(fold.out,"/TemporalSignalResidExclOutliers_",BOOT,".png",sep=""),height=3,width=4,units="in")
    
    #Add to summary table
    index.tree.var<-which(summ.table$trees==trees.in[i] & summ.table$subvar==mysub)
    summ.table$modelfit.end[index.tree.var]<-eq
    summ.table$R2.end[index.tree.var]<-R2
  }
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
  #parse the dates
  dates.list[[i]]<-data.frame(node_name=nm,date=NA)
  for (j in 1:length(nm)){
    dates.list[[i]]$date[j]<-last(unlist(strsplit(nm[j],split="/")))
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
