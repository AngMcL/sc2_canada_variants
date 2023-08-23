#!/usr/bin/env Rscript
# Ancestral state reconstruction of viral geography
# last update 2022-10-21

#usage:
# cd 03_ancestral
# Rscript scripts/AncestralReconstruction.Variant.R "../02_trees/202208_analysis/delta_lsd_out" "../01_subsample/bootsamples_meta_202208"
# Rscript scripts/AncestralReconstruction.Variant.R "../02_trees/BA_analysis/BA1_lsd_out" "../01_subsample/bootsamples_meta_omi"
# Rscript AncestralReconstruction.Variant.R "ft_BA2_root_res/BA2_lsd_out" "BA2_meta"
# Rscript AncestralReconstruction.Variant.R "ft_BA2_unif_root_res/BA2_unif_lsd_out" "BA2_unif_meta"

# use the subsampled data/trees that have had outliers removed based on tempest analysis
# conduct ancestral reconstruction on the tree
# Then, identify introductions (can nodes preceded by non-can nodes)
# split each intro into a sublineage
# parse the lineage from each and then add a .1, .2, .3
# for each sublineage intro, quantify number of desc, whether provincial mixing, from where
# export the parent state of tips 
# export the full set of internal node states and just for canada transitions


#library setup
library(ape)
library(stringr)
library(dplyr)
library(tidyr)
library(phytools)
library(phangorn)
library(forcats)
library(gtools)

#### Setup inputs and outputs ####
args<-commandArgs(trailingOnly=TRUE)
trees.in.f<-args[1] 
meta.in.f<-args[2]

#manual input
#trees.in.f<-"../02_trees/ft2.2_lsd_out"
#meta.in.f<-"../01_subsample/bootsamples"
#trees.in.f<-"../02_trees/202208_analysis/delta_lsd_out"
#meta.in.f<-"../01_subsample/bootsamples_meta_202208"

#tree inputs
trees.in<-list.files(path=trees.in.f,pattern=".timetree.nex",full.names = T)
no.want<-str_which(trees.in, ".nex.binary.nex|nex.binary.tre")
if (length(no.want)>0){trees.in<-trees.in[-no.want]}
trees.in<-mixedsort (trees.in)
b<-length(trees.in) #number of bootstraps
BOOTS<-1:b

if(length(trees.in)<1){
 print("trees not found");break
}

#subsampled meta input corresponding to each bootstrap tree
#note that this requires the same naming base in meta
meta.in<-str_replace_all(trees.in, trees.in.f,meta.in.f) #changes folder
meta.in<-str_replace_all(meta.in, ".timetree.nex",".csv") #changes csv
meta.in<-str_replace_all(meta.in, "_mask_root_res.fasta","") #removes fluff 

#alternatively, pull all the meta files in a subfolder
# meta.in<-list.files(path=meta.in.f, pattern=".csv", full.names=T)
# meta.in<-mixedsort(meta.in)

#CHECK THAT THIS WORKED
if(!all(file.exists(meta.in))){
 print("all meta not found");break
}

#outputs
#pull the folder base (i.e. ft) for output folder
ff.out<-last(unlist(strsplit(trees.in.f,"/")))
f.out<-paste(unlist(strsplit(ff.out,"_"))[[1]], "_ace_out",sep="")
if(!dir.exists(f.out)){dir.create(f.out)}

#output files
#pull the output names from the tree names
base<-str_replace_all(trees.in, trees.in.f, f.out)
base<-str_replace_all(base, "_root_res.fasta","")

all.anc.out<-str_replace_all(base,"timetree.nex","all.anc.csv")
can.anc.out<-str_replace_all(base,"timetree.nex","can.anc.csv")
sublin.out<-str_replace_all(base,"timetree.nex","sublin.csv")
meta.out<-str_replace_all(base,"timetree.nex","meta.csv")
can.states.out<-str_replace_all(base,"timetree.nex","tip.anc.csv")

#### Loop workflow through each tree ####
i=1
for (i in 1:b){
 print(paste("Reading in tree ",i, sep=""))
 #read output from iqtree
 tree<-read.nexus(trees.in[[i]])
 
 #force binary (necessary for recontstr.)
 #note that this is a random resolution and will result in different trees run to run 
 set.seed(111)
 tree<-multi2di(tree)
 
 #need to eliminate branch length==0 with v v small values
 tree$edge.length[tree$edge.length==0]<-max(nodeHeights(tree))*1e-8
 tree$edge.length[tree$edge.length<0]<-max(nodeHeights(tree))*1e-8
 
 #write for future use back into input folder
 # if(!file.exists(paste(trees.in[[i]],".binary.tre",sep=""))){
 write.tree(tree,paste(trees.in[[i]],".binary.tre",sep="")) 
 # }
 
 # if(!file.exists(paste(trees.in[[i]],".binary.nex",sep=""))){
 writeNexus(tree,paste(trees.in[[i]],".binary.nex",sep="")) 
 #for some reason, this file doesn't appear to be binary although the newick above is
 # }
 
 ## Read in subsampled meta 
 print("Reading meta")
 meta<-read.csv(meta.in[[i]])
 
 #change some colnames and reduce size
 #NEED TO CHECK THE FORMAT OF THIS META
 meta<-meta[,c("new.names", "Accession.ID","Variant","date","region","country","division", "Lineage", "Clade")]
 colnames(meta)<-c("tip.label", "GISAID_ID","Variant","date","region","country","division", "Lineage", "GISAID_clade")
 
 #Replace any apostrophes in the newnames (issue fixed in trees, aligns)
 meta$tip.label<-str_replace_all(meta$tip.label, "\\'","")
 #Replace any apostrophes in the country
 meta$country<-str_replace_all(meta$country, "\\'","")
 
 #Replace any ? with _ (that's how fasttree changed tipnames, carried forward)
 meta$tip.label<-str_replace_all(meta$tip.label, "\\?","_")
 
 ## Exclude meta not represented in bootstrap tree (temp outliers, etc)
 # tree$tip.label[which(!tree$tip.label %in% meta.b$tip.label)]
 these<-which(!meta$tip.label %in% tree$tip.label)
 meta.b<-meta
 if (length(these)>0){
   meta.b<-meta.b[-these,]
 }

 print(all(meta.b$tip.label %in% tree$tip.label))
 if(!all(tree$tip.label %in% meta.b$tip.label)){print("tips and meta names no matchy")}
 
 ## Add a state column, which will be country for most but for canadian rows add prov also
 for (j in 1:nrow(meta.b)){
   meta.b$state[j]<-meta.b$country[j]
   if (meta.b$state[j]=="Canada") {meta.b$state[j]<-paste(meta.b$state[j], meta.b$division[j], sep="_")}
 }
 
 #make sure meta.b in same order as the tree
 meta.b<-meta.b[order(match(meta.b[,"tip.label"], tree$tip.label)),]
 
 #### Use trees to reconstruct ancestral state ####
 ## run ace, ML method, equilib
 print("Running ASR (slowww)")
 anc.b<-ace(x=meta.b$state,tree,type="discrete",method="ML",model="ER")
 
 #pull the likelihood of each state at each node into a df
 anc.df.b<-as.data.frame(anc.b$lik.anc)
 
 #### Identify Canadian sublineages ####
 #ie internal nodes with non-Can parental node
 #use the anc object to identify which nodes have majority Can support
 #add node number in order, see http://www.phytools.org/Cordoba2017/ex/8/Anc-states-discrete.html
 print("Finding sublineages")
 n.tip.b<-Ntip(tree)
 anc.df.b$node <- n.tip.b + 1:(n.tip.b-1) 
 
 #pull the majority (Max Likeli state) for each node 
 anc.df.b$node.state<-NA
 anc.df.b$node.lik<-NA
 
 for (j in 1:nrow(anc.df.b)){
   max<-max(anc.df.b[j,1:(ncol(anc.df.b)-3)])
   maxcol<-which(anc.df.b[j,]==max)[1]
   anc.df.b$node.state[j]<-colnames(anc.df.b)[maxcol]
   anc.df.b$node.lik[j]<-anc.df.b[j,maxcol]
 }
 
 #ID the Canadian internal nodes and quantify
 can<-str_which(anc.df.b$node.state,"Canada")
 # length(can); length(tree$tip.label) 
 
 
 ##Query states for canadian nodes parents and desc
 #prime col number for locations
 colz<-ncol(anc.df.b)-3 #exclude node info
 
 #for each of these, look at parent internal node and descendant tips
 anc.df.b$node.par<-NA
 anc.df.b$node.desc<-NA
 
 anc.df.b$par.state<-NA
 anc.df.b$desc.state<-NA
 
 anc.df.b$desc.n<-NA
 anc.df.b$desc.accession<-NA
 
 anc.df.b$par.lik<-NA
 anc.df.b$desc.lik<-NA
 
 #loop through canadian nodes and query parent, descendants, likelihoods
 ####Note could change this to loop through all, but will result in a bigger object
 for (j in 1:length(can)){
   par<-Ancestors(tree, node=anc.df.b$node[can[j]], type="parent") #parent node ID
   desc<-as.numeric(c(Descendants(tree, node=anc.df.b$node[can[j]]))[[1]]) #desc nodes
   desc.c<-paste0(as.character(desc),collapse=", ")
   
   #nodes of desc and par
   anc.df.b$node.desc[can[j]]<-desc.c #descendents as a character string
   anc.df.b$node.par[can[j]]<-par #parent
   
   #parent's geography state
   anc.df.b$par.state[can[j]]<- anc.df.b$node.state [which(anc.df.b$node==par)]
   #descendant geography state
   anc.df.b$desc.state[can[j]]<-paste0(meta.b$state[desc],collapse=", ") #meta sorted the same as tree
   
   #number of descendants
   anc.df.b$desc.n[can[j]]<-length(desc)
   #Accession ID of descendants
   anc.df.b$desc.accession[can[j]]<-paste0(meta.b$GISAID_ID[desc],collapse=", ")
   
   #likelihood of desc and par
   anc.df.b$par.lik[can[j]]<-anc.df.b$node.lik [which(anc.df.b$node==par)]
   anc.df.b$desc.lik[can[j]]<-paste0(anc.df.b$node.lik[desc],collapse=", ")
 }
 
 ## Subset state object to Can states only
 can.anc.b<-anc.df.b[can,]
 
 #check
 # table(can.anc.b$node.state)
 
 #exclude rows if parent node was also Canadian (not an intro to Canada)
 #quantify interprovincial transm overall in subseq. script
 # nrow(can.anc.b)
 if (any(str_detect(can.anc.b$par.state, "Canada"))) {
   can.anc.b<-can.anc.b[-str_which(can.anc.b$par.state, "Canada"),]
 }
 #note that anc.df.b still has all the parental node states involving canada (domestic transmissions)
 # nrow(can.anc.b)
 
 # ADD a new column for whether nested or not, then can exclude on this basis later if desired
 can.anc.b$nested<-NA
 for (j in 1:nrow(can.anc.b)){
   dd<-can.anc.b$node.desc[j]
   if (j>1){
     #if all node descendants are found in any other rows, then nested
     if (mean(str_detect(can.anc.b$node.desc[-j], dd))>0) can.anc.b$nested[j]<-"nested"
   }
 }
 
 #make this into a df with a row for each introduction 
 sum.b<-can.anc.b[,c("node","node.state","node.lik","node.par",
                     "par.state","par.lik",
                     "desc.n","nested",
                     "node.desc","desc.state", "desc.accession","desc.lik")]
 
 #summarize desc state into one character
 for (j in 1:nrow(sum.b)){
   t<-table(unlist(strsplit(sum.b$desc.state[j],split=", ")))
   des<-c()
   for (k in 1:length(t)){
     n<-paste(names(t)[k], ": ",t[k],sep="")
     des<-c(des,n)
   }
   sum.b$desc.state[j]<-paste0(des,collapse="; ")
 }
 
 colnames(sum.b)<-c("Node","Node.Location","Node.Likelihood",
                    "Parent.Node","Parent.Location","Parent.Likelihood",
                    "Number.Descendants","Nested",
                    "Descendant.Node","Descendant.Location","Descendant.AccessionID","Descendant.Likelihood")
 
 #order rows in descending order of n.desc and name them
 sum.b<-sum.b[with(sum.b, order(Number.Descendants, decreasing = TRUE)),]
 
 ## Assign lineages to sum.b sublineages using pango updated Lineage
 # Note: moved sublineage naming into analyze states because we want
 # to know the first sample date (which requires pulling the LSD-inferred dates)
 
 #Give each one a majority lineage designation
 sum.b$Lineage<-NA
 sum.b$Lin.Support<-NA #percent of descendants with same lineage 
 #note could have derivatives inside
 
 #Pull the majority lineage and the support for each type within
 for (j in 1:nrow(sum.b)){
   epis<-unlist(strsplit(sum.b$Descendant.AccessionID[j],split=", "))
   lins<-meta.b$Lineage[meta.b$GISAID_ID %in% epis]
   t<-table(lins)
   all<-sum(t)
   t<-rev(sort(t)) #sort by most frequent
   des<-c()
   for (k in 1:length(t)){
     n<-paste(names(t)[k], ": ",round(t[k]/all*100,digits=2),"%",sep="")
     des<-c(des,n)
   }
   sum.b$Lin.Support[j]<-paste0(des,collapse="; ")
   sum.b$Lineage[j]<-names(rev(sort(t)))[1] 
 }
 
 #### Identify Canadian tips' direct ancestral (parental) state ####
 #note, not the same thing as singletons
 print("Pulling Canadian tips' ancestry")
 #identify Can tips, labs, and states
 can.tips<-match(meta.b$tip.label[meta.b$country=="Canada"], tree$tip.label)
 can.tiplabs<-tree$tip.label[can.tips]
 can.tip.state<-meta.b$state[match(can.tiplabs,meta.b$tip.label)]
 
 #identify parental nodes
 can.par.nodes<-Ancestors(tree,node=can.tips,type="parent")
 
 #start a df for all this
 can.states.df<-data.frame("tip.id"=can.tips, "tip.label"=can.tiplabs, "tip.state"=can.tip.state,"par.id"=can.par.nodes,"par.state"=NA, "par.lik"=NA)
 
 #go through each row, pull the parent state and its likelihood
 for (j in 1:nrow(can.states.df)){
   par<-can.states.df$par.id[j]
   #lookup the node state in anc.df
   can.states.df$par.state[j]<-anc.df.b$node.state[anc.df.b$node==par]
   can.states.df$par.lik[j]<-anc.df.b$node.lik[anc.df.b$node==par]
 }
 #### write these objects #### 
 
 #NO OVERWRITES
 print("Writing objects")
 if (!file.exists(can.anc.out[i])){
   write.csv(anc.df.b,all.anc.out[i],row.names = F) #full object of internal nodes after ASR
   write.csv(can.anc.b,can.anc.out[i],row.names = F) #only Canadian internal nodes
   write.csv(sum.b, sublin.out[i],row.names = F) #each row describes a Canadian sublineage 
   write.csv(meta.b,meta.out[i],row.names = F) #meta corresponding exactly to ASR
   write.csv(can.states.df, can.states.out[i],row.names = F) #Canadian tip states' ancestry
 }
}# final loop closure over boot i