#!/usr/bin/env Rscript
#' ---
#' title: "Downsample by global prev"
#' author: "Angela"
#' date: "21/07/2020"
#' output: github_document
#' ---
#' 
#' OMICRON SPECIFIC: separate individual lineages BA.1 BA.1.1 BA.2
#' 
#' ## Objectives 
#' * Subsample alignment using the global country-specific estimated VOC case counts (from the day of most recent sample) by month
#' * Subsample Canadian sequences proportionally to provinces' burden
#' * Export alignments and corresponding metadata 
#' 
#' ## big changes
#' * use globl regions instead of individual countries
#' * use estimated relative contribution of VOC-specific case counts
#' * mimic temporal distribution of VOC-specific case counts
#' * read in objects exported from 01_summarizeCaseSeqsByGeoLin instead of generating them here
#' 
## ----setup, include=FALSE----------------------------------------------------------------------------------
library(knitr)
library(tidyverse)
library(stringr)
library(ape)
library(coronavirus)
library(ggplot2)
library(RColorBrewer)
library(Biostrings)
library(lubridate)
library(gtools)
library(dplyr)

BiocManager::install("Biostrings")
#' 
#' ## Change inputs 
## ----------------------------------------------------------------------------------------------------------
data.date<-"2022-03-22" # the day of most recent sample
b<-10 #number of bootstraps to subsample
n.samp<-50000 #total number of sequences to subsample
#ideally want to take 50,000/2 = 25000 n can seqs, with remainder global
n.can<-n.samp/2
n.glob<-n.samp-n.can

### DATA INPUT ###
#data is on external harddrive, partitioned into gisaid clades
align.in.folder<-"/Volumes/PhD_2/2022_SC2_Can/20220322_gisaid/00_cleanData/GISAID/sequences/03_cleaner"
align.in.list<-list.files(align.in.folder,full.names = T) %>% mixedsort()

#merged and separated (into global v can) meta
#versions with dates too early removed
case.fold.in<-"Results_PreSubsamp_omi_2023"
canada.meta.in<-paste0(case.fold.in,"/meta.can.csv")
global.meta.in<-paste0(case.fold.in,"/meta.glob.csv")
# global.meta.in<-"CleanMeta_PreSubsamp/global_meta.csv"
# canada.meta.in<-"CleanMeta_PreSubsamp/canada_meta.csv"

#number of partitions for alignments
n.P<-length(align.in.list)

#also wuhan-hu-1 separately
wuh.fas<-"/Volumes/PhD_2/2022_SC2_Can/20220322_gisaid/00_cleanData/GISAID/wuhan-hu-1/EPI_ISL_402125.fasta"


#inputs for monthly VOC-specific case counts and proportions for global and canada
can.monthly.var.in<-paste0(case.fold.in,"/AverageMonthlyVOCcasesByProv.csv")
glob.monthly.var.in<-paste0(case.fold.in,"/AverageMonthlyVOCcasesByCountry.csv")
  
#variant df
variant.df.in<-paste0(case.fold.in,"/variant.df.can.csv")

### OUTPUT ###
#keep subsampled alignments on the harddrive, but meta in the local folder
align.out.folder<-"/Volumes/PhD_2/2022_SC2_Can/20220322_gisaid/01_subsampled/bootsamples_fasta_omi/" 
meta.out.folder<-"bootsamples_meta_omi/"
fig.out.folder<-"Results_Subsamp_omi/"
if(!dir.exists(align.out.folder)){dir.create(align.out.folder)}
if(!dir.exists(meta.out.folder)){dir.create(meta.out.folder)}
if(!dir.exists(fig.out.folder)){dir.create(fig.out.folder)}

#the variants of the analysis have >100 sequences from canada
analysis.vars<-c("BA.1","BA.2","BA.1.1")

#key outputs
meta.out.BA.1<-paste(meta.out.folder,"BA.1_",1:b,".csv",sep="")
fasta.out.BA.1<-paste(align.out.folder,"BA.1_",1:b,".fasta",sep="")

meta.out.BA.2<-paste(meta.out.folder,"BA.2_",1:b,".csv",sep="")
fasta.out.BA.2<-paste(align.out.folder,"BA.2_",1:b,".fasta",sep="")

meta.out.BA.1.1<-paste(meta.out.folder,"BA.1.1_",1:b,".csv",sep="")
fasta.out.BA.1.1<-paste(align.out.folder,"BA.1.1_",1:b,".fasta",sep="")

### set up months for dataframes below ### 
last.mon<-format(as.Date(data.date), "%m")
#For 2021 & 2022 
last.mon.l<-24+as.numeric(last.mon) #starts at 0, so actually leads to 13+..

# mons<-c(0:last.mon) #changed Dec 2019 to 0
mons <- format(ymd(as.Date("2019-12-01")) %m+% months(0:last.mon.l),"%Y-%m")
l.mons<-length(mons)

#2020 months (and 2019)
months.2020<-mons[str_which(mons,"2020|2019")]
l.mons.2020<-length(months.2020)
#2021 months (and 2022)
months.2021<-mons[str_which(mons,"2021|2022")]
l.mons.2021<-length(months.2021)

#' 
#' ## Read in the meta
## ----------------------------------------------------------------------------------------------------------
#fully merged and cleaned meta + metacan
meta.glob<-read.csv(global.meta.in)
meta.can<-read.csv(canada.meta.in)

#mod for meta.glob

omi.ba1<-which(meta.glob$var.WHO=="Omicron" & meta.glob$Lineage=="BA.1")
omi.ba11<-which(meta.glob$var.WHO=="Omicron" & meta.glob$Lineage=="BA.1.1")
omi.ba2<-which(meta.glob$var.WHO=="Omicron" & meta.glob$Lineage=="BA.2")
meta.glob$var.WHO[omi.ba1]<-"BA.1"
meta.glob$var.WHO[omi.ba11]<-"BA.1.1"
meta.glob$var.WHO[omi.ba2]<-"BA.2"
# table(meta.glob$Lineage[meta.glob$var.WHO=="Omicron"]) #BA.1 BA.1.1   BA.2 and 2 None
#remove the lineage==none within omicron, n=2
no.omi<-which(meta.glob$var.WHO=="Omicron" & meta.glob$Lineage=="None")
if(length(no.omi)>0){
  meta.glob<-meta.glob[-no.omi,]
}
print(nrow(meta.glob))


#check the incomplete dates, NA dates
na.dates<-which(is.na(meta.can$date))
for (i in 1:length(na.dates)){
  dt<-unlist(strsplit(unlist(strsplit(meta.can$new.names[na.dates[i]],split="/"))[[2]],"-"))
  if(length(dt)==1){print(paste0("only year ",i))}
  if(length(dt)==2){
  #find days in that month
  dayz<-lubridate::days_in_month(paste(dt[1],dt[2],"01",sep="-")) %>% as.integer
  #Then, sample a day and plop it in the date
  day.samp<-formatC(sample(dayz,size=1), width = 2, format = "d", flag = "0") #make sure day has trailing zero
  meta.can$date[na.dates[i]]<-paste(dt[1],dt[2],day.samp,sep="-") 
  }
}
#check dealt with
length(which(is.na(meta.can$date)))==0

#remove any NA months because only year
# na.month.glob<-which(is.na(meta.glob$month))
na.month.can<-which(is.na(meta.can$month))
if(length(na.month.can)>0){
  meta.can<-meta.can[-na.month.can,]
}


#' 
#' ## Read in the VOC-specific case counts and proportions
## ----------------------------------------------------------------------------------------------------------
#read it
can.monthly.var<-read.csv(can.monthly.var.in)

#only keep analysis vars
can.monthly.var<-can.monthly.var[which(can.monthly.var$var.WHO %in% analysis.vars),]

#remove NA months
no.month<-which(is.na(can.monthly.var$month))
if(length(no.month)>0){
  can.monthly.var<-can.monthly.var[-no.month,]
}

#make a reduced version of the df
can.monthly.var.red<-can.monthly.var[,c("var.WHO","division","month","prop")]
colnames(can.monthly.var.red)[4]<-"probability"

#####GLOBAL####
#repeat for global
glob.monthly.var<-read.csv(glob.monthly.var.in)
# head(glob.monthly.var)

#remove NA months
no.month<-which(is.na(glob.monthly.var$month))
if(length(no.month)>0){
  glob.monthly.var<-glob.monthly.var[-no.month,]
}

#make a reduced version of the df
glob.monthly.var.red<-glob.monthly.var[,c("var.WHO","country","month","prop")]
colnames(glob.monthly.var.red)[4]<-"probability"

### TEMPORAL ###

#total estimated VOC cases per month in all of canada 
can.total.monthly.var<-can.monthly.var %>% filter(division=="Quebec") #arbitrary , contains total
can.total.monthly.var<-can.total.monthly.var[,c("var.WHO","month","total_monthly_VOC_cases")]
#calculate proportion of VOC cases in that month
can.total.var<-can.total.monthly.var %>% group_by(var.WHO) %>% summarize(total_VOC_cases= sum(total_monthly_VOC_cases,na.rm=T)) 
can.total.monthly.var<-can.total.monthly.var%>%left_join(can.total.var,by="var.WHO")
can.total.monthly.var<-can.total.monthly.var%>% mutate(proportion_VOC_cases=total_monthly_VOC_cases/total_VOC_cases)

#actual available VOC seqs availabe per month
meta.can.avail<-meta.can %>% group_by(var.WHO, month) %>% count() %>% as.data.frame()
colnames(meta.can.avail)<-c("var.WHO","month","total_VOC_seqs")
#change NAs to zeroes

#join on the df
can.total.monthly.var<-can.total.monthly.var%>% left_join(meta.can.avail, by=c("var.WHO","month"))
can.total.monthly.var$total_VOC_seqs<-na.replace(can.total.monthly.var$total_VOC_seqs,replace = 0)
# head(can.total.monthly.var)

#extract total estimated VOC cases per month in all of global
glob.total.monthly.var<-glob.monthly.var %>% filter(country=="USA")#arbitrary b/c all have totals
glob.total.monthly.var<-glob.total.monthly.var[,c("var.WHO","month","total_monthly_VOC_cases")]

#calculate proportion of VOC cases in that month
glob.total.var<-glob.total.monthly.var %>% group_by(var.WHO) %>% summarize(total_VOC_cases= sum(total_monthly_VOC_cases,na.rm=T)) 
glob.total.monthly.var<-glob.total.monthly.var%>%left_join(glob.total.var,by="var.WHO")
glob.total.monthly.var<-glob.total.monthly.var%>% mutate(proportion_VOC_cases=total_monthly_VOC_cases/total_VOC_cases)

#actual available VOC seqs availabe per month
meta.glob.avail<-meta.glob %>% group_by(var.WHO, month) %>% count() %>% as.data.frame()
colnames(meta.glob.avail)<-c("var.WHO","month","total_VOC_seqs")
#change NAs to zeroes

#join on the df
glob.total.monthly.var<-glob.total.monthly.var%>% left_join(meta.glob.avail, by=c("var.WHO","month"))
glob.total.monthly.var$total_VOC_seqs<-na.replace(glob.total.monthly.var$total_VOC_seqs,replace = 0)
# head(glob.total.monthly.var)


#' 
#' ## join sampling probabilities on the meta
## ----------------------------------------------------------------------------------------------------------
#Join to meta
meta.can<-meta.can %>% left_join(can.monthly.var.red, by=c("var.WHO","division", "month"))
meta.glob<-meta.glob %>% left_join(glob.monthly.var.red, by=c("var.WHO", "country","month"))

#Change any NA probabilities
na.prob<-which(is.na(meta.can$probability))
if(length(na.prob)>0){
 meta.can$probability[na.prob]<-0.00000000000001
}
#make sure no probabilities are zero
zero.prob<-which(meta.can$probability==0)
if(length(zero.prob)>0){
  meta.can$probability[zero.prob]<-0.00000000000001
}
#and none are negative
neg.prob<-which(meta.can$probability<0)
if(length(neg.prob)>0){
  meta.can$probability[neg.prob]<-0.00000000000001
}


#Change any NA probabilities
na.prob<-which(is.na(meta.glob$probability))
if(length(na.prob)>0){
 meta.glob$probability[na.prob]<-0.00000000000001
}
#make sure no probabilities are zero
zero.prob<-which(meta.glob$probability==0)
if(length(zero.prob)>0){
  meta.glob$probability[zero.prob]<-0.00000000000001
}
#and none are negative
neg.prob<-which(meta.glob$probability<0)
if(length(neg.prob)>0){
  meta.glob$probability[neg.prob]<-0.00000000000001
}

#make sure these df have the same structure
#canada meta has the incomp column, so just add this arbitrarily to global
meta.glob$incomp<-1
#match the colname order to canada
meta.glob<-meta.glob[,match(colnames(meta.can),colnames(meta.glob))]
# all(colnames(meta.glob)==colnames(meta.can))


#' 
#' ### VARIANT SPECIFIC SUBSAMPLING ###
#' 
#' ## read in the variant df
## ----------------------------------------------------------------------------------------------------------
voc.df<-read.csv(variant.df.in)
print("reading in variant df")
voc.2<-voc.df %>% filter(!var.WHO %in% c("Theta","Lambda","GH/490R")) 
var.ord<-voc.2$var.WHO
#restrict to analysis vars with more than 500 Canadian sequences
#sort by analysis order, which is alphabetical
voc.df<-voc.df[match(analysis.vars,voc.df$var.WHO),]
print("restricting to analysis variants")

# each var.WHO has output files x b bootstraps
n.vars<-length(analysis.vars)

# add a column for parents
#Parent lineages are only sampled lightly, before the first sample date of the pango or sibling
# voc.df$var.lin
voc.df$Parent<-c( "B.1.1.529|B.1.1|B.1|B","B.1.1.529|B.1.1|B.1|B","B.1.1.529|B.1.1|B.1|B")

#Note an issue with space in date
voc.df$first.glob.date<-trimws(voc.df$first.glob.date)

#save this again with the parents in there
write.csv(voc.df,paste0(case.fold.in,"/voc.df.parents.csv"))

#' 
#' ## Reduce meta to variants or parents
## ----------------------------------------------------------------------------------------------------------
## Reduce the size of meta and meta can to only include the variants in the analysis
nrow(meta.glob) #7357327
nrow(meta.can) # 274420

#Make a parents global df
all.pars<-unique(unlist(strsplit(as.vector(voc.df$Parent),"\\|")))
meta.glob.par<-meta.glob[which(meta.glob$Lineage %in% all.pars),]

#Reduce to variants of analysis
print("reducing dfs to variants")
meta.glob<-meta.glob[which(meta.glob$var.WHO %in% analysis.vars),]
meta.can<-meta.can[which(meta.can$var.WHO %in% analysis.vars),]
print(nrow(meta.glob)) 
print(nrow(meta.can)) 

#' 
#' 
#' # Subsample Canadian sequences + remainder up to 50000 global sequences 
#' ## allocate sequence sampling distribution over time
## ----------------------------------------------------------------------------------------------------------
# Setup up number of seqs to sample per month table, based on proportion keeping
# In this version, rather than equally distribute sequences, mimic the cases over time
# VOC<-"Delta"
make.sample.t.VOCs.can<-function(VOC){
  t.seq<-can.total.monthly.var[can.total.monthly.var$var.WHO==VOC,-1]
  #remove months with value <1
  under.one<-which(t.seq$total_VOC_seqs<1)
  if(length(under.one)>0){
      t.seq<-t.seq[-under.one,]
  }
  
  #if the total number of sequences available is greater than n.can, take all
  if(sum(t.seq$total_VOC_seqs)<=n.can) {
    t.seq$seq.take<-t.seq$total_VOC_seqs
    return(t.seq[,c("month","total_monthly_VOC_cases","total_VOC_seqs","seq.take")])
    break}
  
  #multiply the proportion of cases by month by the total number of sequences desired
  t.seq$seq.take<-floor(t.seq$proportion_VOC_cases*n.can) #round down to whole integer
  
  #should take double n seqs where sparse and an early month
  too.few<-which(t.seq$seq.take<500 & t.seq$seq.take<t.seq$total_VOC_seqs)
  #only early months (first half of entries)
  too.few<-too.few[which(too.few < nrow(t.seq)/2)] 
  #the months without too few, for later distributing of seqs
  not.too.few<-which(!t.seq$month %in% t.seq$month[too.few])
  #if number of seq available is less than 500, take all, else take 500
  if(length(too.few)>0){
    for(j in 1:length(too.few)){
      avlbl<-t.seq$total_VOC_seqs[too.few[j]]
      ifelse(avlbl<500,
             t.seq$seq.take[too.few[j]]<-t.seq$total_VOC_seqs[too.few[j]],
             t.seq$seq.take[too.few[j]]<-500)
    }
  }
  #redistribute difference (will be a smaller sum) to remaining months
  total.remaining<-n.can - sum(t.seq$seq.take [too.few])
  # sum(t.seq$seq.take)
  #calculate a new proportion of cases among months where there were not too few sequences
  total.cases.remaining<-sum(t.seq$total_monthly_VOC_cases[not.too.few])
  t.seq$proportion_VOC_cases[not.too.few]<-t.seq$total_monthly_VOC_cases[not.too.few]/total.cases.remaining
  # sum(t.seq$proportion_VOC_cases[not.too.few])
  #calculate seqs to take based on those proportions
  t.seq$seq.take[not.too.few]<-floor(t.seq$proportion_VOC_cases[not.too.few]*total.remaining)
  # sum(t.seq$seq.take)
    
  #do any of these desired # of seqs sampled per month exceed what is available?
  sh<-0
  while(sh==0){
    too.many<-which(t.seq$seq.take[not.too.few] > t.seq$total_VOC_seqs[not.too.few]) 
    not.too.many<-which(t.seq$seq.take[not.too.few] <= t.seq$total_VOC_seqs[not.too.few]) 
    if(length(too.many)<1){sh<-1}
    #if any exceed what is available
    if (length(too.many)>=1){ 
      t.seq$seq.take[not.too.few[too.many]]<-t.seq$total_VOC_seqs[not.too.few[too.many]]
      #redistribute difference to remaining months,  proportion of cases they represent among themselves
      total.remaining<-n.can - sum(t.seq$seq.take [not.too.few[too.many]], t.seq$seq.take[too.few])
      total.cases.remaining<-sum(t.seq$total_monthly_VOC_cases[not.too.few[not.too.many]])
      t.seq$proportion_VOC_cases[not.too.few[not.too.many]]<-t.seq$total_monthly_VOC_cases[not.too.few[not.too.many]]/total.cases.remaining
      # sum(t.seq$proportion_VOC_cases[not.too.few[not.too.many]])
      t.seq$seq.take[not.too.few[not.too.many]]<-round(t.seq$proportion_VOC_cases[not.too.few[not.too.many]] * total.remaining,digits=0)
    }
  }
  
  # sum(t.seq$seq.take)

  #if still don't have enough seqs (because of rounding above), add the diff to middle months
  if (sum(t.seq$seq.take)<n.can & length(not.too.many)>0) {
    diff<-as.numeric(n.can)-sum(t.seq$seq.take)
    mid<-floor(length(not.too.many)/2)
    selection<-not.too.few[not.too.many[mid:(mid+diff)]]
    if(any(is.na(selection))){selection<-na.omit(selection) %>% as.vector()}
    t.seq$seq.take[selection]<-t.seq$seq.take[selection]+1
  }
  # sum(t.seq$seq.take)
  return(t.seq[,c("month","total_monthly_VOC_cases","total_VOC_seqs","seq.take")])
}

t.BA.1<-make.sample.t.VOCs.can("BA.1")
t.BA.1.1<-make.sample.t.VOCs.can("BA.1.1")
t.BA.2<-make.sample.t.VOCs.can("BA.2")


# Summarize the sequences to take
t<-bind_rows( t.BA.1,t.BA.1.1, t.BA.2)
t$VOC<-c(rep("BA.1",nrow(t.BA.1)),rep("BA.2",nrow(t.BA.2)), rep("BA.1.1",nrow(t.BA.1.1)))
t$seqpercase<-t$seq.take/t$total_monthly_VOC_cases

t.summary<-t %>% group_by(VOC) %>% summarize(sum.take=sum(seq.take), sum.avail=sum(total_VOC_seqs),sum.cases=sum(total_monthly_VOC_cases))
t<-left_join(t,t.summary, by="VOC")
#calculate proportional densities for each
t<-t %>% mutate(prop.cases=total_monthly_VOC_cases/sum.cases,
                          prop.take=seq.take/sum.take,
                          prop.avail=total_VOC_seqs/sum.avail)

write.csv(t,paste0(fig.out.folder,"CanadaSubSample_table_VOCs.csv"))

#' 
## ----------------------------------------------------------------------------------------------------------
#plot these out 
# run this to cheat!:
t<-read.csv(paste0(fig.out.folder,"CanadaSubSample_table_VOCs.csv"))
t<-t[,-1]
tomi<-t

t.long<-pivot_longer(t,cols=c(10,11,12),names_to = "param",values_to = "value")


plot.can.dis<-ggplot(t.long)+
  geom_density(aes(x=month, y=value,color=param,group=param),stat="identity",
               position=position_jitter(width = 0.1))+
  facet_wrap(~VOC,scales="free_y")+
  theme(axis.text.x=element_text(angle=90,size=rel(0.5)))+
  labs(x=NULL, y="Density")
plot.can.dis
ggsave(paste0(fig.out.folder,"CanadaSubSample_distribution.png"),height=2,width=5)

ggplot(t)+
  geom_line(aes(x=month,y=seqpercase,group=VOC),color="black")+
  facet_wrap(~VOC,scales="free_y")+
  theme(axis.text.x=element_text(angle=90,size=rel(0.5)))+
  labs(x=NULL, y="Sequence per case")

ggsave(paste0(fig.out.folder,"CanadaSubSample_seqpercase.png"),height=6,width=8.5)

#' 
#' ## subsample according to provinces' monthly VOC case proportions to the Canadian total
## ----------------------------------------------------------------------------------------------------------
#input number of bootstraps b, the table of samps per month, the variant
sample.metacan.withT<-function(t, VOC) { 
  #make a reduced df for the variant, ensure dates are after the first global date
  df.seqs<-meta.can[meta.can$var.WHO==VOC, ]
  #build a list, with df for each index across b subsamples
  samp.meta<-replicate(n=b, vector())
  for (i in 1:b){
    can.samp<-vector()
    set.seed(i) #makes reproducible run to run
    #for each month in the table of seqs for a lineage, extract # seqs identified above
    for (j in 1:nrow(t)){
      row.mm<-which(df.seqs$month==t$month[j])
      #take sample of size specified in table t, but if short, will just take that many
      can.samp<-c(can.samp, sample(df.seqs$new.names[row.mm],
                        size=t$seq.take[j], replace=FALSE,
                        prob=df.seqs$probability[row.mm]) )
    }
    #add these to the list item
    samp.meta[[i]]<-df.seqs[which(df.seqs$new.names %in% can.samp),]
    samp.meta[[i]]$replicate<-as.character(i)
    }
  return(samp.meta) # a list of sampled dfs, length=replciates
} #end of function

### Apply function 
samp.can.BA.1<-sample.metacan.withT(t=t.BA.1, "BA.1")
samp.can.BA.2<-sample.metacan.withT(t=t.BA.2, "BA.2")
samp.can.BA.1.1<-sample.metacan.withT(t=t.BA.1.1, "BA.1.1")

#number of seqs taken for each
# sum(t.BA.1$seq.take) #25k
# sum(t.BA.1.1$seq.take) #25k
# sum(t.BA.2$seq.take) #


#' 
#' # Sample remainder up to 50K global
#' ## allocate sequence distribution over time
## ----------------------------------------------------------------------------------------------------------
#Take up to 25000 for nglob or maximum available
#Take either equal number of sequences from global and Canada 
calculate.global.n<-function(VOC){
  n.glob<-n.can #25000
  glob.avail<-length(which(meta.glob$var.WHO==VOC))
  if(glob.avail<n.glob){nglob<-glob.avail}
  return(n.glob)
}

n.glob.BA.1<-calculate.global.n("BA.1")
n.glob.BA.2<-calculate.global.n("BA.2")
n.glob.BA.1.1<-calculate.global.n("BA.1.1")

print(n.glob.BA.1)
print(n.glob.BA.1.1)
print(n.glob.BA.2)

#test for below
# n.glob.voc<-n.glob.BA.1
# VOC<-"BA.1"

# Make a table of sequences to pull for each month 
print("Making sequence tables for global")
make.sample.t.VOCs.glob<-function(n.glob.voc,VOC){
  t.seq<-glob.total.monthly.var[glob.total.monthly.var$var.WHO==VOC,-1]
  #remove months with value <1 or nA
  under.one<-which(t.seq$total_VOC_seqs<1)
  if(length(under.one)>0){
      t.seq<-t.seq[-under.one,]
  }
  na.months<-which(is.na(t.seq$total_VOC_seqs))
  if(length(na.months)>0){
    t.seq<-t.seq[-na.months,]
  }
  #if the total number of sequences available is greater than n.glob.voc, take all
  print(sum(t.seq$total_VOC_seqs))
  print(n.glob.voc)
  if(sum(t.seq$total_VOC_seqs)<=n.glob.voc) {
    t.seq$seq.take<-t.seq$total_VOC_seqs
    return(t.seq[,c("month","total_monthly_VOC_cases","total_VOC_seqs","seq.take")])
    break}
  
  #multiply the proportion of cases by month by the total number of sequences desired
  t.seq$seq.take<-floor(t.seq$proportion_VOC_cases*n.glob.voc) #round down to whole integer
  
  #should take double n seqs where sparse and an early month
  threshold<-500
  if(n.glob.voc<nrow(t.seq)*500){threshold<-floor(n.glob.voc/nrow(t.seq))}
  too.few<-which(t.seq$seq.take<threshold & t.seq$seq.take<=t.seq$total_VOC_seqs)
  #only early months (first half of entries)
  too.few<-too.few[which(too.few < nrow(t.seq)/2)] 
  #the months without too few, for later distributing of seqs
  not.too.few<-which(!t.seq$month %in% t.seq$month[too.few])
  #if number of seq available is less than threshold, take all, else take threshold
  if(length(too.few)>0){
    for(j in 1:length(too.few)){
      avlbl<-t.seq$total_VOC_seqs[too.few[j]]
      ifelse(avlbl<threshold,
             t.seq$seq.take[too.few[j]]<-t.seq$total_VOC_seqs[too.few[j]],
             t.seq$seq.take[too.few[j]]<-threshold)
    }
  }
  #redistribute difference (will be a smaller sum) to remaining months
  total.remaining<-n.glob.voc - sum(t.seq$seq.take [too.few])
  # sum(t.seq$seq.take)
  #calculate a new proportion of cases among months where there were not too few sequences
  total.cases.remaining<-sum(t.seq$total_monthly_VOC_cases[not.too.few])
  t.seq$proportion_VOC_cases[not.too.few]<-t.seq$total_monthly_VOC_cases[not.too.few]/total.cases.remaining
  # sum(t.seq$proportion_VOC_cases[not.too.few])
  #calculate seqs to take based on those proportions
  t.seq$seq.take[not.too.few]<-floor(t.seq$proportion_VOC_cases[not.too.few]*total.remaining)
  # sum(t.seq$seq.take)
    
  #do any of these desired # of seqs sampled per month exceed what is available?
  # sh<-0
  # while(sh==0){
    too.many<-which(t.seq$seq.take[not.too.few] > t.seq$total_VOC_seqs[not.too.few]) 
    not.too.many<-which(t.seq$seq.take[not.too.few] < t.seq$total_VOC_seqs[not.too.few]) 
    # if(length(too.many)<1){sh<-1}
    #if any exceed what is available
    if (length(too.many)>=1){ 
      #take what is available if too many 
      t.seq$seq.take[not.too.few[too.many]]<-t.seq$total_VOC_seqs[not.too.few[too.many]]
      #redistribute difference to remaining months,  proportion of cases they represent among themselves
      total.remaining<-n.glob.voc - sum(t.seq$seq.take [not.too.few[too.many]], t.seq$seq.take[too.few])
      total.cases.remaining<-sum(t.seq$total_monthly_VOC_cases[not.too.few[not.too.many]])
      t.seq$proportion_VOC_cases[not.too.few[not.too.many]]<-t.seq$total_monthly_VOC_cases[not.too.few[not.too.many]]/total.cases.remaining
      # sum(t.seq$proportion_VOC_cases[not.too.few[not.too.many]])
      new.takers<-floor(t.seq$proportion_VOC_cases[not.too.few[not.too.many]] * total.remaining)
      if(total.remaining > sum(new.takers)){total.remaining<-sum(new.takers)}
      new.takers.2<-floor(t.seq$proportion_VOC_cases[not.too.few[not.too.many]] * total.remaining)
      new.takers.2[which(new.takers.2>t.seq$total_VOC_seqs[not.too.few[not.too.many]])]<-t.seq$total_VOC_seqs[which(new.takers.2>t.seq$total_VOC_seqs[not.too.few[not.too.many]])]
      t.seq$seq.take[not.too.few[not.too.many]]<-new.takers.2
    }
  # }
  
  # sum(t.seq$seq.take)

  #if still don't have enough seqs (because of rounding above), add the diff to months with extra seqs
  if (sum(t.seq$total_VOC_seqs)>n.glob.voc & sum(t.seq$seq.take)<n.glob.voc & length(not.too.many)>0) {
    diff<-as.numeric(n.glob.voc)-sum(t.seq$seq.take)
    extra<-which(t.seq$total_VOC_seqs>t.seq$seq.take)
    if(length(extra)<diff){diff<-length(extra)} #if there are more diffs than extra months, just add one to each extra
    t.seq$seq.take[extra[1:diff]]<-t.seq$seq.take[extra[1:diff]]+1
  }
  # sum(t.seq$seq.take)
    
  #any over? just take what is available and call it a day.
  still.too.many<-which(t.seq$seq.take>t.seq$total_VOC_seqs)
  if(length(still.too.many)>0){
    t.seq$seq.take[still.too.many]<-t.seq$total_VOC_seqs[still.too.many]
  }
  return(t.seq[,c("month","total_monthly_VOC_cases","total_VOC_seqs","seq.take")])
}

#apply the function to allocate sequences for each VOC
t.glob.BA.1<-make.sample.t.VOCs.glob(n.glob.BA.1, "BA.1")
t.glob.BA.2<-make.sample.t.VOCs.glob(n.glob.BA.2, "BA.2")
t.glob.BA.1.1<-make.sample.t.VOCs.glob(n.glob.BA.1.1, "BA.1.1")

#export a summary of this distrib (avail vs taken)
t.glob<-bind_rows(t.glob.BA.1,t.glob.BA.2, t.glob.BA.1.1)
t.glob$VOC<-c(rep("BA.1",nrow(t.glob.BA.1)),rep("BA.2",nrow(t.glob.BA.2)),rep("BA.1.1",nrow(t.glob.BA.1.1)))
t.glob$seqpercase<-t.glob$seq.take/t.glob$total_monthly_VOC_cases

t.glob.summary<-t.glob %>% group_by(VOC) %>% summarize(sum.take=sum(seq.take), sum.avail=sum(total_VOC_seqs),sum.cases=sum(total_monthly_VOC_cases))
t.glob<-left_join(t.glob,t.glob.summary, by="VOC")
#calculate proportional densities for each
t.glob<-t.glob %>% mutate(prop.cases=total_monthly_VOC_cases/sum.cases,
                          prop.take=seq.take/sum.take,
                          prop.avail=total_VOC_seqs/sum.avail)
write.csv(t.glob,paste0(fig.out.folder,"GLOBALSubSample_table_VOCs.csv"))
print("t.glob written")
#' 
## ----------------------------------------------------------------------------------------------------------

# run this to cheat!:
t.glob<-read.csv(paste0(fig.out.folder,"GLOBALSubSample_table_VOCs.csv"))
t.glob<-t.glob[,-1]

#plot these out 
t.glob$VOC<-factor(t.glob$VOC,levels=var.ord )
t.glob.long<-pivot_longer(t.glob,cols=c(10,11,12),names_to = "param",values_to = "value")

p.glob.dis<-ggplot(t.glob.long)+
  geom_density(aes(x=month, y=value,color=param,group=param),stat="identity",
               position=position_jitter(width = 0.1))+
  facet_wrap(~VOC,scales="free_y")+
  theme(axis.text.x=element_text(angle=90,size=rel(0.5)))+
  labs(x=NULL, y="Density")
p.glob.dis
ggsave(paste0(fig.out.folder,"GlobalSubSample_distribution.png"),height=2,width=5)

ggplot(t.glob)+
  geom_line(aes(x=month,y=seqpercase,group=VOC),color="black")+
  facet_wrap(~VOC,scales="free_y")+
  theme(axis.text.x=element_text(angle=90,size=rel(0.5)))+
  labs(x=NULL, y="Sequence per case")
ggsave(paste0(fig.out.folder,"CanadaSubSample_seqpercase.png"),height=6,width=8.5)

t.globomi<-t.glob


## Make a grob
library(cowplot)
plot.can.dis.noleg<-plot.can.dis+theme(legend.position = "none")
plot_grid(plot.can.dis.noleg, p.glob.dis, labels=c("A", "B"),
          rel_widths=c(0.8,1),nrow=1)
ggsave(paste0(fig.out.folder,"CanadaGlobal_distribution.png"),width=7,height=2,units="in")
#' 
#' 
#' ## subsample global sequences
## ----------------------------------------------------------------------------------------------------------
print("subsampling global seqs")
sample.metaglob.withT<-function(t, VOC){ #input number of bootstraps replicates, the table of samps per month
  
  #make a reduced df for the variant
  df.seqs<-meta.glob[meta.glob$var.WHO==VOC, ]
  #if any pango lineage ==none (omicron), remove
  no.lin<-which(df.seqs$Lineage=="None")
  if(length(no.lin)>0){df.seqs<-df.seqs[-no.lin,]}
  
  #empty list of dfs, one per replicate
  samp.meta<-replicate(n=b, vector())
  for (i in 1:b){
    glob.samp<-vector()
    set.seed(i) 
    #for each month in the table of seqs for a lineage, extract # seqs identified above
    for (j in 1:nrow(t)){
      row.mm<-which(df.seqs$month==t$month[j])
      #take sample of size specified in table t, but if short, will just take that many
      glob.samp<-c(glob.samp, sample(df.seqs$new.names[row.mm],
                        size=t$seq.take[j], replace=FALSE,
                        prob=df.seqs$probability[row.mm]) )
    }
    #add these to the list item
    samp.meta[[i]]<-df.seqs[which(df.seqs$new.names %in% glob.samp),]
    samp.meta[[i]]$replicate<-as.character(i)
    }
  return(samp.meta) # a list of sampled dfs, length=replciates
} #end of function

### Apply function 
samp.glob.BA.1<-sample.metaglob.withT(t=t.glob.BA.1, "BA.1")
samp.glob.BA.2<-sample.metaglob.withT(t=t.glob.BA.2, "BA.2")
samp.glob.BA.1.1<-sample.metaglob.withT(t=t.glob.BA.1.1, "BA.1.1")

#check it out
# nrow(samp.glob.alpha[[1]])
# nrow(samp.glob.beta[[1]])
# nrow(samp.glob.gamma[[1]])
# nrow(samp.glob.delta[[1]])
# nrow(samp.glob.omicron[[1]])

#Check lineages
# table(samp.glob.alpha[[1]]$Lineage)
# table(samp.glob.beta[[1]]$Lineage)
# table(samp.glob.eta[[1]]$Lineage)
# table(samp.glob.gamma[[1]]$Lineage)
# table(samp.glob.delta[[1]]$Lineage)
# table(samp.glob.omicron[[1]]$Lineage)

#' 
#' #### SAMPLE the parental lineages from global
## ----------------------------------------------------------------------------------------------------------
#only includes parental lineages sampled before first VOC detection
#could also increase probability of selecting sequences from same geo as VOC early detection?
#also, for B context, always add wuhan-hu-1 (can still root on it that way)
# VOC<-"Other"
print("sample parental lineages")
sample.metaglob.parent<-function(VOC){ 
  #a reduced df for the variant parents          
  parents<-unlist(strsplit(voc.df$Parent[voc.df$var.WHO==VOC],"\\|"))
  par.rows<-c()
  for (i in 1:length(parents)){
    par.glob<-which(meta.glob.par$Lineage==parents[i] & 
                    meta.glob.par$date < voc.df$first.glob.date[voc.df$var.WHO==VOC]) #BEFORE earliest date of variant
    if(VOC=="Other"){par.glob<-which(meta.glob.par$Lineage==parents[i])} #a workaround for other, which has an early global
    par.rows<-c(par.rows,par.glob)
  }
  df.seqs.par<-meta.glob.par[par.rows,]
  
  # table(df.seqs.par$Lineage)
  samp.meta<-replicate(n=b, vector())
  for (i in 1:b){
    glob.samp<-vector()
    set.seed(i) 
    #for each unique parental lineage, take ten sequences
    parz<-unique(df.seqs.par$Lineage)
    for (j in 1:length(parz)){
      #Only include 1 B sequence, wuhan-hu-1 
      if (parz[j]=="B"){
        B.samp<-df.seqs.par$new.names[which(df.seqs.par$Accession.ID=="EPI_ISL_402125")]
        glob.samp<-c(glob.samp,B.samp); next
      }
      #check there are at least sz=ten parent seqs available, otherwise take what's available (i.e. b.1.617)
      sz<-10
      avl<-length(df.seqs.par$new.names[df.seqs.par$Lineage==parz[j]])
      if (avl<sz){sz<-avl}
      #Otherwise, take ten sequences fully at random
      par.samp<-sample(df.seqs.par$new.names[df.seqs.par$Lineage==parz[j]],
                       size=sz, replace=F)
      glob.samp<-c(glob.samp,par.samp)
    }
    
    #add these to the list item
    samp.meta[[i]]<-df.seqs.par[which(df.seqs.par$new.names %in% glob.samp),]
    samp.meta[[i]]$replicate<-as.character(i)
    }
  return(samp.meta) # a list of sampled dfs, length=replciates
} #end of function

### Apply function 
samp.glob.par.BA.1<-sample.metaglob.parent("BA.1")
samp.glob.par.BA.2<-sample.metaglob.parent("BA.2")
samp.glob.par.BA.1.1<-sample.metaglob.parent("BA.1.1")

# nrow(samp.glob.par.alpha[[1]])
# nrow(samp.glob.par.gamma[[1]])
# nrow(samp.glob.par.delta[[1]])
# nrow(samp.glob.par.omicron[[1]])
# table(samp.glob.par.alpha[[1]]$Lineage)
# table(samp.glob.par.gamma[[1]]$Lineage)
# table(samp.glob.par.delta[[1]]$Lineage)
# table(samp.glob.par.omicron[[1]]$Lineage)

#' 
#' ## Merge the Canada and global subsamples (and parent) together 
## ----------------------------------------------------------------------------------------------------------
#function that operates over three lists to bind dataframes 
bind.dfs.in.lists<-function(list1,list2,list3,b){
  newlist<-replicate(b, vector())
  for (i in 1:b){
    #make sure all the same class of date, character
    list1[[i]]$date<-as.character(list1[[i]]$date)
    list2[[i]]$date<-as.character(list2[[i]]$date)
    list3[[i]]$date<-as.character(list3[[i]]$date)

    #bind rows with identical columns
    newlist[[i]]<-dplyr::bind_rows(list1[[i]],list2[[i]])
    newlist[[i]]<-dplyr::bind_rows(newlist[[i]], list3[[i]])
    print(nrow(newlist[[i]])) #should==n.samp=50000
  }
  return(newlist) #return a single list of length b 
}

samp.full.BA.1<-bind.dfs.in.lists(samp.glob.par.BA.1, samp.glob.BA.1, samp.can.BA.1, b)
samp.full.BA.2<-bind.dfs.in.lists(samp.glob.par.BA.2, samp.glob.BA.2, samp.can.BA.2, b)
samp.full.BA.1.1<-bind.dfs.in.lists(samp.glob.par.BA.1.1, samp.glob.BA.1.1, samp.can.BA.1.1, b)

# #check lineages
# nrow(samp.full.alpha[[1]])
# table(samp.full.alpha[[1]]$Lineage)
# table(samp.full.beta[[1]]$Lineage)
# table(samp.full.gamma[[1]]$Lineage)
# table(samp.full.delta[[1]]$Lineage)
# table(samp.full.omicron[[1]]$Lineage)

#' 
#' ## WRITE the metadata
## ----------------------------------------------------------------------------------------------------------
for (i in 1:b){
  write.csv(samp.full.BA.1[[i]],meta.out.BA.1[i])
  write.csv(samp.full.BA.2[[i]],meta.out.BA.2[i])
  write.csv(samp.full.BA.1.1[[i]],meta.out.BA.1.1[i])
}

nrow(samp.full.BA.1[[i]])
head(samp.full.BA.1[[i]]$new.names)

#' 
## ----------------------------------------------------------------------------------------------------------
#fully compiled and cleaned meta + metacan
# write.csv(meta.glob,"~/../../Volumes/PhD_2/2022_SC2_Can/20220322_gisaid/00_cleanData/GISAID/global_meta_presamp_omi.csv")
# write.csv(meta.can,"~/../../Volumes/PhD_2/2022_SC2_Can/20220322_gisaid/00_cleanData/GISAID/canada_meta_presamp_omi.csv")

#' 
#' ## Clear the memory
## ----------------------------------------------------------------------------------------------------------
rm(meta.glob)
rm(meta.can)

rm(meta.can.avail)
rm(meta.glob.par)
rm(meta.glob.avail)

rm(samp.glob.BA.1)
rm(samp.glob.BA.1.1)
rm(samp.glob.BA.2)

rm(samp.can.BA.1)
rm(samp.can.BA.1.1)
rm(samp.can.BA.2)

rm(samp.glob.par.BA.1)
rm(samp.glob.par.BA.1.1)
rm(samp.glob.par.BA.2)

#read these in for each loop to reduce mem
rm(samp.full.BA.1)
rm(samp.full.BA.2)
rm(samp.full.BA.1.1)

#' 
#' ## Generate bootstrap alignments from bootstrap subsampled metadata
## ----------------------------------------------------------------------------------------------------------
#make alignments using sequence names in the bootstrap dataframes from both the global and canada samples
#DONT overwrite things all willy nilly
# if(file.exists(fasta.out.alpha[1])){print("CAUTION - overwriting existing files")}

#each one starts with wuhan-hu-1 (part of parent subsampling)
print("make subsamp'd alignments")
wu.fas<-readDNAStringSet(wuh.fas)
for (i in 1:b){
  writeXStringSet(wu.fas, fasta.out.BA.1[i],format="fasta",append=F)
  writeXStringSet(wu.fas, fasta.out.BA.2[i],format="fasta",append=F)
  writeXStringSet(wu.fas, fasta.out.BA.1.1[i],format="fasta",append=F)
}

#loop through alignment partitions, pulling only what's needed, concatenating/appending onto files
for (k in 1:n.P){ #n.P = number of partitions
  
  #Read in the partition
  print(paste("Reading in ",align.in.list[k],sep=""))
  align<-readDNAStringSet(align.in.list[k])

  print("Sampling and writing for BA.1")
  for (i in 1:b){
    samp.full.BA.1<-read.csv(meta.out.BA.1[i])
    align.BA.1<-align [which(names(align) %in% samp.full.BA.1$new.names)]
    writeXStringSet(align.BA.1, fasta.out.BA.1[i],format="fasta",append=T)
    rm(align.BA.1)
    rm(samp.full.BA.1)
  }

  print("Sampling and writing for BA.2")
  for (i in 1:b){
    samp.full.BA.2<-read.csv(meta.out.BA.2[i])
    align.BA.2<-align [which(names(align) %in% samp.full.BA.2$new.names)]
    writeXStringSet(align.BA.2, fasta.out.BA.2[i],format="fasta",append=T)
    rm(align.BA.2)
    rm(samp.full.BA.2)
  }


  print("Sampling and writing for BA.1.1")
  for (i in 1:b){
    samp.full.BA.1.1<-read.csv(meta.out.BA.1.1[i])
    align.BA.1.1<-align [which(names(align) %in% samp.full.BA.1.1$new.names)]
    writeXStringSet(align.BA.1.1, fasta.out.BA.1.1[i],format="fasta",append=T)
    rm(align.BA.1.1)
    rm(samp.full.BA.1.1)
  }

  print("Done sampling for this partition")
  #clear the mem
  rm(align)
} # END Of n.P partition loop
print("Sampling complete")
#clear the memory


##Re-import and export alignments for bootstraps in proper format
#export  fasta, SLOWWW
for (i in 1:b){
  print("Re-writing fasta")
  ## Workaround to fix the gaps added by Biostrings::writeXstringset
  #read in and overwrite the previous files using the ape commands
  write.FASTA(read.FASTA(fasta.out.BA.1[i],type = "DNA"), fasta.out.BA.1[i],append = F)
  write.FASTA(read.FASTA(fasta.out.BA.2[i],type = "DNA"), fasta.out.BA.2[i],append = F)
  write.FASTA(read.FASTA(fasta.out.BA.1.1[i],type = "DNA"), fasta.out.BA.1.1[i],append = F)
}
 
