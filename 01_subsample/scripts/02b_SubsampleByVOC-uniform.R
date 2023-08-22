#!/usr/bin/env Rscript
#' ---
#' title: "Downsample by global prev"
#' author: "Angela"
#' date: "21/07/2020"
#' output: github_document
#' ---
#' 
#' UNIFORM TEMPORAL SAMPLING INSTEAD OF PROPORTIONAL TO CASES BY MONTH - for comparison of improvment of correlations
#' use the scripts from the old builds (2020 analysis)
#' change output folder 
#' 
#' ## Objectives 
#' * Subsample Canadian sequences proportionally to provinces' burden
#' * Export alignments and corresponding metadata 
#' 
#' ## big changes
#' * use globl regions instead of individual countries
#' * use estimated relative contribution of VOC-specific case counts for probabilities
#' * DO NOT mimic temporal distribution of VOC-specific case counts
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
canada.meta.in<-"Results_PreSubsamp/meta.can.csv"
global.meta.in<-"Results_PreSubsamp/meta.glob.csv"
# global.meta.in<-"CleanMeta_PreSubsamp/global_meta.csv"
# canada.meta.in<-"CleanMeta_PreSubsamp/canada_meta.csv"

#number of partitions for alignments
n.P<-length(align.in.list)

#also wuhan-hu-1 separately
wuh.fas<-"/Volumes/PhD_2/2022_SC2_Can/20220322_gisaid/00_cleanData/GISAID/wuhan-hu-1/EPI_ISL_402125.fasta"

#inputs for monthly VOC-specific case counts and proportions for global and canada
can.monthly.var.in<-"Results_PreSubsamp/AverageMonthlyVOCcasesByProv.csv"
glob.monthly.var.in<-"Results_PreSubsamp/AverageMonthlyVOCcasesByCountry.csv"
  
#variant df
variant.df.in<-"Results_PreSubsamp/variant.df.can.csv"

### OUTPUT ###
#keep subsampled alignments on the harddrive, but meta in the local folder
align.out.folder<-"/Volumes/PhD_2/2022_SC2_Can/20220322_gisaid/01_subsampled/bootsamples_fasta_unif/" 
meta.out.folder<-"bootsamples_meta_unif/"
fig.out.folder<-"Results_Subsamp_unif/"
if(!dir.exists(align.out.folder)){dir.create(align.out.folder)}
if(!dir.exists(meta.out.folder)){dir.create(meta.out.folder)}
if(!dir.exists(fig.out.folder)){dir.create(fig.out.folder)}

#the variants of the analysis have >100 sequences from canada
analysis.vars<-c("Other","Alpha","Beta","Delta","Epsilon","Eta","Gamma","Iota","Kappa","Mu","Omicron","Zeta")

#key outputs
meta.out.other<-paste(meta.out.folder,"other_",1:b,".csv",sep="")
fasta.out.other<-paste(align.out.folder,"other_",1:b,".fasta",sep="")

meta.out.alpha<-paste(meta.out.folder,"alpha_",1:b,".csv",sep="")
fasta.out.alpha<-paste(align.out.folder,"alpha_",1:b,".fasta",sep="")

meta.out.beta<-paste(meta.out.folder,"beta_",1:b,".csv",sep="")
fasta.out.beta<-paste(align.out.folder,"beta_",1:b,".fasta",sep="")

meta.out.delta<-paste(meta.out.folder,"delta_",1:b,".csv",sep="")
fasta.out.delta<-paste(align.out.folder,"delta_",1:b,".fasta",sep="")

meta.out.epsilon<-paste(meta.out.folder,"epsilon_",1:b,".csv",sep="")
fasta.out.epsilon<-paste(align.out.folder,"epsilon_",1:b,".fasta",sep="")

meta.out.eta<-paste(meta.out.folder,"eta_",1:b,".csv",sep="")
fasta.out.eta<-paste(align.out.folder,"eta_",1:b,".fasta",sep="")

meta.out.gamma<-paste(meta.out.folder,"gamma_",1:b,".csv",sep="")
fasta.out.gamma<-paste(align.out.folder,"gamma_",1:b,".fasta",sep="")

meta.out.iota<-paste(meta.out.folder,"iota_",1:b,".csv",sep="")
fasta.out.iota<-paste(align.out.folder,"iota_",1:b,".fasta",sep="")

meta.out.kappa<-paste(meta.out.folder,"kappa_",1:b,".csv",sep="")
fasta.out.kappa<-paste(align.out.folder,"kappa_",1:b,".fasta",sep="")

meta.out.mu<-paste(meta.out.folder,"mu_",1:b,".csv",sep="")
fasta.out.mu<-paste(align.out.folder,"mu_",1:b,".fasta",sep="")

meta.out.omicron<-paste(meta.out.folder,"omicron_",1:b,".csv",sep="")
fasta.out.omicron<-paste(align.out.folder,"omicron_",1:b,".fasta",sep="")

meta.out.zeta<-paste(meta.out.folder,"zeta_",1:b,".csv",sep="")
fasta.out.zeta<-paste(align.out.folder,"zeta_",1:b,".fasta",sep="")

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
head(glob.total.monthly.var)


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
meta.glob$incomp<-0
#match the colname order to canada
meta.glob<-meta.glob[,match(colnames(meta.can),colnames(meta.glob))]
# all(colnames(meta.glob)==colnames(meta.can))


#' 
#' ### VARIANT SPECIFIC SUBSAMPLING ###
#' 
#' ## read in the variant df
## ----------------------------------------------------------------------------------------------------------
voc.df<-read.csv(variant.df.in)

voc.2<-voc.df %>% filter(!var.WHO %in% c("Theta","Lambda","GH/490R")) 
var.ord<-voc.2$var.WHO
#restrict to analysis vars with more than 500 Canadian sequences
#sort by analysis order, which is alphabetical
voc.df<-voc.df[match(analysis.vars,voc.df$var.WHO),]

# each var.WHO has output files x b bootstraps
n.vars<-length(analysis.vars)

# add a column for parents
#Parent lineages are only sampled lightly, before the first sample date of the pango or sibling
# voc.df$var.lin
voc.df$Parent<-c("B","B.1.1|B.1|B","B.1|B", "B.1.617|B.1|B", 
                 "B.1|B","B.1|B","B.1.1.28|B.1.1|B.1|B",
                  "B.1|B", "B.1.617|B.1|B","B.1|B",
                 "B.1.1.529|B.1.1|B.1|B","B.1.1.28|B.1.1|B.1|B")

#Note an issue with space in date
voc.df$first.glob.date<-trimws(voc.df$first.glob.date)

#save this again with the parents in there
write.csv(voc.df,"Results_PreSubsamp/voc.df.parents.csv")

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
meta.glob<-meta.glob[which(meta.glob$var.WHO %in% analysis.vars),]
meta.can<-meta.can[which(meta.can$var.WHO %in% analysis.vars),]
nrow(meta.glob) 
nrow(meta.can) 

#' 
#' 
#' # Subsample Canadian sequences + remainder up to 50000 global sequences 
#' ## allocate sequence sampling distribution over time
## ----------------------------------------------------------------------------------------------------------
# In this version, UNIFORMLY/equally distribute sequences as much as possible
# VOC<-"Alpha"

make.sample.t.VOCs.can.unif<-function(VOC){
  t.seq<-can.total.monthly.var[can.total.monthly.var$var.WHO==VOC,-1]
  #remove months with value <1
  under.one<-which(t.seq$total_VOC_seqs<1)
  if(length(under.one)>0){
      t.seq<-t.seq[-under.one,]
  }
  
  #if the total number of sequences available is greater than n.can, take all
  if(sum(t.seq$total_VOC_seqs)<=n.can) {
    t.seq$seq.take<-t.seq$total_VOC_seqs
    print(sum(t.seq$seq.take,na.rm = T))
    return(t.seq[,c("month","total_monthly_VOC_cases","total_VOC_seqs","seq.take")])
    break}
  
  #distribute the desired sequences between the number of months with sequences
  mm<-trunc(n.can/nrow(t.seq))
  t.seq$seq.take<-mm #default take
  
  #check which are "short months" where all seqs will be taken
  too.few<-which(t.seq$total_VOC_seqs<mm)
  short<-t.seq$month[too.few]
  t.seq$seq.take[too.few]<-t.seq$total_VOC_seqs[too.few] #take all where too few
  sh<-0
  while(sh==0){ #until sh=1, optimize number of seqs for months to be evenly distributed where possible
    t.sh<-t.seq[which(!t.seq$month %in% short), ]
    #calculate month-specific mm.2, the number of seqs to take each month for non-short months
    mm.2<-trunc((n.can- sum(t.seq$seq.take[t.seq$month %in% short],na.rm=T)) / nrow(t.sh) )
    #which are still short?
    st.short<-t.sh$month[which(t.sh$total_VOC_seqs < mm.2)]
    st.okay<-t.sh$month[which(t.sh$total_VOC_seqs >= mm.2)]
    t.seq$seq.take[t.seq$month%in%st.short]<-t.seq$total_VOC_seqs[t.seq$month%in%st.short]
    t.seq$seq.take[t.seq$month%in%st.okay]<-mm.2
    #if all remaining months have that many seqs, great, sample that many
    if(length(st.short)==0){sh<-1;break} #if none are short, break the loop
    if(length(st.short)>0){
      short<-c(short,st.short)} #if at least one still short, add it to short list and repeat
  }
  
  #if still don't have enough seqs (because of truncating above), add the diff
  if (sum(t.seq$seq.take,na.rm = T)<n.can) {
    diff<-as.numeric(n.can)-sum(t.seq$seq.take)
    t.seq$seq.take[which(!t.seq$month%in%short)[1:diff]]<-t.seq$seq.take[which(!t.seq$month%in%short)[1:diff]] + 1
  }

  print(sum(t.seq$seq.take))
  
  return(t.seq[,c("month","total_monthly_VOC_cases","total_VOC_seqs","seq.take")])
}

t.other<-make.sample.t.VOCs.can.unif("Other")
t.alpha<-make.sample.t.VOCs.can.unif("Alpha")
t.beta<-make.sample.t.VOCs.can.unif("Beta")
t.delta<-make.sample.t.VOCs.can.unif("Delta")
t.epsilon<-make.sample.t.VOCs.can.unif("Epsilon")
t.eta<-make.sample.t.VOCs.can.unif("Eta")
t.gamma<-make.sample.t.VOCs.can.unif("Gamma")
t.kappa<-make.sample.t.VOCs.can.unif("Kappa")
t.iota<-make.sample.t.VOCs.can.unif("Iota")
t.mu<-make.sample.t.VOCs.can.unif("Mu")
t.omicron<-make.sample.t.VOCs.can.unif("Omicron")
t.zeta<-make.sample.t.VOCs.can.unif("Zeta")

# Summarize the sequences to take
t<-bind_rows(t.other, t.alpha, t.beta, t.delta, t.epsilon, t.eta, t.gamma, t.kappa, t.iota, t.mu, t.omicron, t.zeta)
t$VOC<-c(rep("Other",nrow(t.other)), rep("Alpha",nrow(t.alpha)), rep("Beta",nrow(t.beta)),  rep("Delta",nrow(t.delta)),
         rep("Epsilon",nrow(t.epsilon)), rep("Eta",nrow(t.eta)),  rep("Gamma",nrow(t.gamma)),
         rep("Iota",nrow(t.iota)), rep("Kappa",nrow(t.kappa)),  rep("Mu",nrow(t.mu)),
         rep("Omicron",nrow(t.omicron)),rep("Zeta",nrow(t.zeta)))
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
t$VOC<-factor(t$VOC,levels=var.ord )
t.long<-pivot_longer(t,cols=c(10,11,12),names_to = "param",values_to = "value")

ggplot(t.long)+
  geom_density(aes(x=month, y=value,color=param,group=param),stat="identity")+
  facet_wrap(~VOC,scales="free_y")+
  theme(axis.text.x=element_text(angle=90,size=rel(0.5)))+
  labs(x=NULL, y="Density")
ggsave(paste0(fig.out.folder,"CanadaSubSample_distribution.png"),height=6,width=8.5)

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
samp.can.other<-sample.metacan.withT(t=t.other, "Other")
samp.can.alpha<-sample.metacan.withT(t=t.alpha, "Alpha")
samp.can.beta<-sample.metacan.withT(t=t.beta, "Beta")
samp.can.delta<-sample.metacan.withT(t=t.delta, "Delta")
samp.can.epsilon<-sample.metacan.withT(t=t.epsilon, "Epsilon")
samp.can.eta<-sample.metacan.withT(t=t.eta, "Eta")
samp.can.gamma<-sample.metacan.withT(t=t.gamma, "Gamma")
samp.can.iota<-sample.metacan.withT(t=t.iota, "Iota")
samp.can.kappa<-sample.metacan.withT(t=t.kappa, "Kappa")
samp.can.mu<-sample.metacan.withT(t=t.mu, "Mu")
samp.can.omicron<-sample.metacan.withT(t=t.omicron, "Omicron")
samp.can.zeta<-sample.metacan.withT(t=t.zeta, "Zeta")

#number of seqs taken for each
# nrow(samp.can.alpha[[1]])==sum(t.alpha$seq.take)
# sum(t.alpha$seq.take) #25k
# sum(t.beta$seq.take) #1414 (all)
# sum(t.gamma$seq.take) #13734 (all)
# sum(t.delta$seq.take) #25k
# sum(t.omicron$seq.take) #25k

#Check sampled lineages
# table(samp.can.alpha[[1]]$Lineage)
# table(samp.can.beta[[1]]$Lineage)
# table(samp.can.gamma[[1]]$Lineage)
# table(samp.can.delta[[1]]$Lineage)
# table(samp.can.omicron[[1]]$Lineage)

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
n.glob.other<-calculate.global.n("Other")
n.glob.alpha<-calculate.global.n("Alpha")
n.glob.beta<-calculate.global.n("Beta")
n.glob.delta<-calculate.global.n("Delta")
n.glob.epsilon<-calculate.global.n("Epsilon")
n.glob.eta<-calculate.global.n("Eta")
n.glob.gamma<-calculate.global.n("Gamma")
n.glob.iota<-calculate.global.n("Iota")
n.glob.kappa<-calculate.global.n("Kappa")
n.glob.mu<-calculate.global.n("Mu")
n.glob.omicron<-calculate.global.n("Omicron")
n.glob.zeta<-calculate.global.n("Zeta")

# Make a table of sequences to pull for each month 
n.glob.voc<-n.glob.alpha
VOC<-"Alpha"

make.sample.t.VOCs.glob.unif<-function(n.glob.voc,VOC){
    t.seq<-glob.total.monthly.var[glob.total.monthly.var$var.WHO==VOC,-1]
    #remove months with value <1
    under.one<-which(t.seq$total_VOC_seqs<1)
    if(length(under.one)>0){
      t.seq<-t.seq[-under.one,]
    }
    
    #if the total number of sequences available is greater than n.glob.voc, take all
    if(sum(t.seq$total_VOC_seqs)<=n.glob.voc) {
      t.seq$seq.take<-t.seq$total_VOC_seqs
      print(sum(t.seq$seq.take,na.rm = T))
      return(t.seq[,c("month","total_monthly_VOC_cases","total_VOC_seqs","seq.take")])
      break}
    
    #distribute the desired sequences between the number of months with sequences
    mm<-trunc(n.glob.voc/nrow(t.seq))
    t.seq$seq.take<-mm #default take
    
    #check which are "short months" where all seqs will be taken
    too.few<-which(t.seq$total_VOC_seqs<mm)
    short<-t.seq$month[too.few]
    t.seq$seq.take[too.few]<-t.seq$total_VOC_seqs[too.few] #take all where too few
    sh<-0
    while(sh==0){ #until sh=1, optimize number of seqs for months to be evenly distributed where possible
      t.sh<-t.seq[which(!t.seq$month %in% short), ]
      #calculate month-specific mm.2, the number of seqs to take each month for non-short months
      mm.2<-trunc((n.glob.voc- sum(t.seq$seq.take[t.seq$month %in% short],na.rm=T)) / nrow(t.sh) )
      #which are still short?
      st.short<-t.sh$month[which(t.sh$total_VOC_seqs < mm.2)]
      st.okay<-t.sh$month[which(t.sh$total_VOC_seqs >= mm.2)]
      t.seq$seq.take[t.seq$month%in%st.short]<-t.seq$total_VOC_seqs[t.seq$month%in%st.short]
      t.seq$seq.take[t.seq$month%in%st.okay]<-mm.2
      #if all remaining months have that many seqs, great, sample that many
      if(length(st.short)==0){sh<-1;break} #if none are short, break the loop
      if(length(st.short)>0){
        short<-c(short,st.short)} #if at least one still short, add it to short list and repeat
    }    

    #if still don't have enough seqs (because of truncating above), add the diff
    if (sum(t.seq$seq.take)<n.glob.voc) {
      diff<-as.numeric(n.glob.voc)-sum(t.seq$seq.take)
      t.seq$seq.take[which(!t.seq$month%in%short)[1:diff]]<-t.seq$seq.take[which(!t.seq$month%in%short)[1:diff]] + 1
    }
    
    print(sum(t.seq$seq.take))
    
    return(t.seq[,c("month","total_monthly_VOC_cases","total_VOC_seqs","seq.take")])
}

#apply the function to allocate sequences for each VOC
t.glob.other<-make.sample.t.VOCs.glob.unif(n.glob.other, "Other")
t.glob.alpha<-make.sample.t.VOCs.glob.unif(n.glob.alpha, "Alpha")
t.glob.beta<-make.sample.t.VOCs.glob.unif(n.glob.beta, "Beta")
t.glob.delta<-make.sample.t.VOCs.glob.unif(n.glob.delta, "Delta")
t.glob.epsilon<-make.sample.t.VOCs.glob.unif(n.glob.epsilon, "Epsilon")
t.glob.eta<-make.sample.t.VOCs.glob.unif(n.glob.eta, "Eta")
t.glob.gamma<-make.sample.t.VOCs.glob.unif(n.glob.gamma, "Gamma")
t.glob.iota<-make.sample.t.VOCs.glob.unif(n.glob.iota, "Iota")
t.glob.kappa<-make.sample.t.VOCs.glob.unif(n.glob.kappa, "Kappa")
t.glob.mu<-make.sample.t.VOCs.glob.unif(n.glob.mu, "Mu")
t.glob.omicron<-make.sample.t.VOCs.glob.unif(n.glob.omicron, "Omicron")
t.glob.zeta<-make.sample.t.VOCs.glob.unif(n.glob.zeta, "Zeta")

#export a summary of this distrib (avail vs taken)
t.glob<-bind_rows(t.glob.other, t.glob.alpha, t.glob.beta, t.glob.delta, t.glob.epsilon, t.glob.eta,
                  t.glob.gamma, t.glob.iota, t.glob.kappa, t.glob.mu, t.glob.omicron, t.glob.zeta)
t.glob$VOC<-c(rep("Other",nrow(t.glob.other)), rep("Alpha",nrow(t.glob.alpha)), rep("Beta",nrow(t.glob.beta)),
              rep("Delta",nrow(t.glob.delta)), rep("Epsilon",nrow(t.glob.epsilon)),rep("Eta",nrow(t.glob.eta)),
              rep("Gamma",nrow(t.glob.gamma)), rep("Iota",nrow(t.glob.iota)),rep("Kappa",nrow(t.glob.kappa)),
              rep("Mu",nrow(t.glob.mu)),rep("Omicron",nrow(t.glob.omicron)), rep("Zeta",nrow(t.glob.zeta)))
t.glob$seqpercase<-t.glob$seq.take/t.glob$total_monthly_VOC_cases

t.glob.summary<-t.glob %>% group_by(VOC) %>% summarize(sum.take=sum(seq.take), sum.avail=sum(total_VOC_seqs),sum.cases=sum(total_monthly_VOC_cases))
t.glob<-left_join(t.glob,t.glob.summary, by="VOC")
#calculate proportional densities for each
t.glob<-t.glob %>% mutate(prop.cases=total_monthly_VOC_cases/sum.cases,
                          prop.take=seq.take/sum.take,
                          prop.avail=total_VOC_seqs/sum.avail)
write.csv(t.glob,paste0(fig.out.folder,"GLOBALSubSample_table_VOCs.csv"))

#' 
## ----------------------------------------------------------------------------------------------------------
#plot these out 
t.glob$VOC<-factor(t.glob$VOC,levels=var.ord )
t.glob.long<-pivot_longer(t.glob,cols=c(10,11,12),names_to = "param",values_to = "value")

ggplot(t.glob.long)+
  geom_density(aes(x=month, y=value,color=param,group=param),stat="identity")+
  facet_wrap(~VOC,scales="free_y")+
  theme(axis.text.x=element_text(angle=90,size=rel(0.5)))+
  labs(x=NULL, y="Density")
ggsave(paste0(fig.out.folder,"GlobalSubSample_distribution.png"),height=6,width=8.5)

ggplot(t.glob)+
  geom_line(aes(x=month,y=seqpercase,group=VOC),color="black")+
  facet_wrap(~VOC,scales="free_y")+
  theme(axis.text.x=element_text(angle=90,size=rel(0.5)))+
  labs(x=NULL, y="Sequence per case")
ggsave(paste0(fig.out.folder,"CanadaSubSample_seqpercase.png"),height=6,width=8.5)

#' 
#' 
#' ## subsample global sequences
## ----------------------------------------------------------------------------------------------------------
sample.metaglob.withT<-function(t, VOC){ #input number of bootstraps replicates, the table of samps per month
  
  #make a reduced df for the variant
  df.seqs<-meta.glob[meta.glob$var.WHO==VOC, ]
  
  #if any pango lineage ==none (omicron), remove
  # no.lin<-which(df.seqs$Lineage=="None")
  # if(length(no.lin)>0){df.seqs<-df.seqs[-no.lin,]}
  
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
samp.glob.other<-sample.metaglob.withT(t=t.glob.other, "Other")
samp.glob.alpha<-sample.metaglob.withT(t=t.glob.alpha, "Alpha")
samp.glob.beta<-sample.metaglob.withT(t=t.glob.beta, "Beta")
samp.glob.delta<-sample.metaglob.withT(t=t.glob.delta, "Delta")
samp.glob.epsilon<-sample.metaglob.withT(t=t.glob.epsilon, "Epsilon")
samp.glob.eta<-sample.metaglob.withT(t=t.glob.eta, "Eta")
samp.glob.gamma<-sample.metaglob.withT(t=t.glob.gamma, "Gamma")
samp.glob.iota<-sample.metaglob.withT(t=t.glob.iota, "Iota")
samp.glob.kappa<-sample.metaglob.withT(t=t.glob.kappa, "Kappa")
samp.glob.mu<-sample.metaglob.withT(t=t.glob.mu, "Mu")
samp.glob.omicron<-sample.metaglob.withT( t=t.glob.omicron, "Omicron")
samp.glob.zeta<-sample.metaglob.withT( t=t.glob.zeta, "Zeta")

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
samp.glob.par.other<-sample.metaglob.parent("Other")
samp.glob.par.alpha<-sample.metaglob.parent("Alpha")
samp.glob.par.beta<-sample.metaglob.parent("Beta")
samp.glob.par.delta<-sample.metaglob.parent("Delta")
samp.glob.par.epsilon<-sample.metaglob.parent("Epsilon")
samp.glob.par.eta<-sample.metaglob.parent("Eta")
samp.glob.par.gamma<-sample.metaglob.parent("Gamma")
samp.glob.par.iota<-sample.metaglob.parent("Iota")
samp.glob.par.kappa<-sample.metaglob.parent("Kappa")
samp.glob.par.mu<-sample.metaglob.parent("Mu")
samp.glob.par.omicron<-sample.metaglob.parent( "Omicron")
samp.glob.par.zeta<-sample.metaglob.parent( "Zeta")

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

samp.full.other<-bind.dfs.in.lists(samp.glob.par.other, samp.glob.other, samp.can.other, b)
samp.full.alpha<-bind.dfs.in.lists(samp.glob.par.alpha, samp.glob.alpha, samp.can.alpha, b)
samp.full.beta<-bind.dfs.in.lists(samp.glob.par.beta, samp.glob.beta, samp.can.beta, b)
samp.full.delta<-bind.dfs.in.lists(samp.glob.par.delta, samp.glob.delta, samp.can.delta, b)
samp.full.epsilon<-bind.dfs.in.lists(samp.glob.par.epsilon, samp.glob.epsilon, samp.can.epsilon, b)
samp.full.eta<-bind.dfs.in.lists(samp.glob.par.eta, samp.glob.eta, samp.can.eta, b)
samp.full.gamma<-bind.dfs.in.lists(samp.glob.par.gamma, samp.glob.gamma, samp.can.gamma, b)
samp.full.iota<-bind.dfs.in.lists(samp.glob.par.iota, samp.glob.iota, samp.can.iota, b)
samp.full.kappa<-bind.dfs.in.lists(samp.glob.par.kappa, samp.glob.kappa, samp.can.kappa, b)
samp.full.mu<-bind.dfs.in.lists(samp.glob.par.mu, samp.glob.mu, samp.can.mu, b)
samp.full.omicron<-bind.dfs.in.lists(samp.glob.par.omicron, samp.glob.omicron, samp.can.omicron, b)
samp.full.zeta<-bind.dfs.in.lists(samp.glob.par.zeta, samp.glob.zeta, samp.can.zeta, b)


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
  write.csv(samp.full.other[[i]],meta.out.other[i])
  write.csv(samp.full.alpha[[i]],meta.out.alpha[i])
  write.csv(samp.full.beta[[i]],meta.out.beta[i])
  write.csv(samp.full.delta[[i]],meta.out.delta[i])
  write.csv(samp.full.epsilon[[i]],meta.out.epsilon[i])
  write.csv(samp.full.eta[[i]],meta.out.eta[i])
  write.csv(samp.full.gamma[[i]],meta.out.gamma[i])
  write.csv(samp.full.iota[[i]],meta.out.iota[i])
  write.csv(samp.full.kappa[[i]],meta.out.kappa[i])
  write.csv(samp.full.mu[[i]],meta.out.mu[i])
  write.csv(samp.full.omicron[[i]],meta.out.omicron[i])
  write.csv(samp.full.zeta[[i]],meta.out.zeta[i])
}


#' 
## ----------------------------------------------------------------------------------------------------------
#fully compiled and cleaned meta + metacan
write.csv(meta.glob,"~/../../Volumes/PhD_2/2022_SC2_Can/20220322_gisaid/00_cleanData/GISAID/global_meta_presamp.csv")
write.csv(meta.can,"~/../../Volumes/PhD_2/2022_SC2_Can/20220322_gisaid/00_cleanData/GISAID/canada_meta_presamp.csv")

#' 
#' ## Clear the memory
## ----------------------------------------------------------------------------------------------------------
rm(meta.glob)
rm(meta.can)

#' 
#' ## Generate bootstrap alignments from bootstrap subsampled metadata
## ----------------------------------------------------------------------------------------------------------
#make alignments using sequence names in the bootstrap dataframes from both the global and canada samples
#DONT overwrite things all willy nilly
if(file.exists(fasta.out.alpha[1])){print("CAUTION - overwriting existing files")}

#each one starts with wuhan-hu-1 (part of parent subsampling)
wu.fas<-readDNAStringSet(wuh.fas)
for (i in 1:b){
  writeXStringSet(wu.fas, fasta.out.other[i],format="fasta",append=F)
  writeXStringSet(wu.fas, fasta.out.alpha[i],format="fasta",append=F)
  writeXStringSet(wu.fas, fasta.out.beta[i],format="fasta",append=F)
  writeXStringSet(wu.fas, fasta.out.delta[i],format="fasta",append=F)
  writeXStringSet(wu.fas, fasta.out.epsilon[i],format="fasta",append=F)
  writeXStringSet(wu.fas, fasta.out.eta[i],format="fasta",append=F)
  writeXStringSet(wu.fas, fasta.out.gamma[i],format="fasta",append=F)
  writeXStringSet(wu.fas, fasta.out.iota[i],format="fasta",append=F)
  writeXStringSet(wu.fas, fasta.out.kappa[i],format="fasta",append=F)
  writeXStringSet(wu.fas, fasta.out.mu[i],format="fasta",append=F)
  writeXStringSet(wu.fas, fasta.out.omicron[i],format="fasta",append=F)
  writeXStringSet(wu.fas, fasta.out.zeta[i],format="fasta",append=F)
}

#loop through alignment partitions, pulling only what's needed, concatenating/appending onto files
for (k in 1:n.P){ #n.P = number of partitions
  
  #Read in the partition
  print(paste("Reading in ",align.in.list[k],sep=""))
  align<-readDNAStringSet(align.in.list[k])

  #sample VOC and boots one by one to not exceed memory
  print("Sampling and writing for other")
  for (i in 1:b){
    samp.full.other<-read.csv(meta.out.other[i])
    align.other<-align [which(names(align) %in% samp.full.other$new.names)]
    writeXStringSet(align.other, fasta.out.other[i],format="fasta",append=T)
  }
  rm(align.other)
  rm(samp.full.other)
  
  #sample VOC and boots one by one to not exceed memory
  print("Sampling and writing for Alpha")
  for (i in 1:b){
    samp.full.alpha<-read.csv(meta.out.alpha[i])
    align.alpha<-align [which(names(align) %in% samp.full.alpha$new.names)]
    writeXStringSet(align.alpha, fasta.out.alpha[i],format="fasta",append=T)
  }
  rm(align.alpha)
  rm(samp.full.alpha)

  print("Sampling and writing for beta")
  for (i in 1:b){
    samp.full.beta<-read.csv(meta.out.beta[i])
    align.beta<-align [which(names(align) %in% samp.full.beta$new.names)]
    writeXStringSet(align.beta, fasta.out.beta[i],format="fasta",append=T)
  }
  rm(align.beta)
  rm(samp.full.beta)


  print("Sampling and writing for Delta")
  for (i in 1:b){
    samp.full.delta<-read.csv(meta.out.delta[i])
    align.delta<-align [which(names(align) %in% samp.full.delta$new.names)]
    writeXStringSet(align.delta, fasta.out.delta[i],format="fasta",append=T)
  }
  rm(align.delta)
  rm(samp.full.delta)

  print("Sampling and writing for epsilon")
  for (i in 1:b){
    samp.full.epsilon<-read.csv(meta.out.epsilon[i])
    align.epsilon<-align [which(names(align) %in% samp.full.epsilon$new.names)]
    writeXStringSet(align.epsilon, fasta.out.epsilon[i],format="fasta",append=T)
  }
  rm(align.epsilon)
  rm(samp.full.epsilon)

  print("Sampling and writing for eta")
  for (i in 1:b){
    samp.full.eta<-read.csv(meta.out.eta[i])
    align.eta<-align [which(names(align) %in% samp.full.eta$new.names)]
    writeXStringSet(align.eta, fasta.out.eta[i],format="fasta",append=T)
  }
  rm(align.eta)
  rm(samp.full.eta)

  print("Sampling and writing for Gamma")
  for (i in 1:b){
    samp.full.gamma<-read.csv(meta.out.gamma[i])
    align.gamma<-align [which(names(align) %in% samp.full.gamma$new.names)]
    writeXStringSet(align.gamma, fasta.out.gamma[i],format="fasta",append=T)
  }
  rm(align.gamma)
  rm(samp.full.gamma)

  print("Sampling and writing for iota")
  for (i in 1:b){
    samp.full.iota<-read.csv(meta.out.iota[i])
    align.iota<-align [which(names(align) %in% samp.full.iota$new.names)]
    writeXStringSet(align.iota, fasta.out.iota[i],format="fasta",append=T)
  }
  rm(align.iota)
  rm(samp.full.iota)

  print("Sampling and writing for kappa")
  for (i in 1:b){
    samp.full.kappa<-read.csv(meta.out.kappa[i])
    align.kappa<-align [which(names(align) %in% samp.full.kappa$new.names)]
    writeXStringSet(align.kappa, fasta.out.kappa[i],format="fasta",append=T)
  }
  rm(align.kappa)
  rm(samp.full.kappa)

  print("Sampling and writing for mu")
  for (i in 1:b){
    samp.full.mu<-read.csv(meta.out.mu[i])
    align.mu<-align [which(names(align) %in% samp.full.mu$new.names)]
    writeXStringSet(align.mu, fasta.out.mu[i],format="fasta",append=T)
  }
  rm(align.mu)
  rm(samp.full.mu)

  print("Sampling and writing for Omicron")
  for (i in 1:b){
    samp.full.omicron<-read.csv(meta.out.omicron[i])
    align.omicron<-align [which(names(align) %in% samp.full.omicron$new.names)]
    writeXStringSet(align.omicron, fasta.out.omicron[i],format="fasta",append=T)
  }
  rm(align.omicron)
  rm(samp.full.omicron)

  print("Sampling and writing for zeta")
  for (i in 1:b){
    samp.full.zeta<-read.csv(meta.out.zeta[i])
    align.zeta<-align [which(names(align) %in% samp.full.zeta$new.names)]
    writeXStringSet(align.zeta, fasta.out.zeta[i],format="fasta",append=T)
  }
  rm(align.zeta)
  rm(samp.full.zeta)
  
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
  write.FASTA(read.FASTA(fasta.out.other[i],type = "DNA"), fasta.out.other[i],append = F)
  write.FASTA(read.FASTA(fasta.out.alpha[i],type = "DNA"), fasta.out.alpha[i],append = F)
  write.FASTA(read.FASTA(fasta.out.beta[i],type = "DNA"), fasta.out.beta[i],append = F)
  write.FASTA(read.FASTA(fasta.out.delta[i],type = "DNA"), fasta.out.delta[i],append = F)
  write.FASTA(read.FASTA(fasta.out.epsilon[i],type = "DNA"), fasta.out.epsilon[i],append = F)
  write.FASTA(read.FASTA(fasta.out.eta[i],type = "DNA"), fasta.out.eta[i],append = F)
  write.FASTA(read.FASTA(fasta.out.gamma[i],type = "DNA"), fasta.out.gamma[i],append = F)
  write.FASTA(read.FASTA(fasta.out.iota[i],type = "DNA"), fasta.out.iota[i],append = F)
  write.FASTA(read.FASTA(fasta.out.kappa[i],type = "DNA"), fasta.out.kappa[i],append = F)
  write.FASTA(read.FASTA(fasta.out.mu[i],type = "DNA"), fasta.out.mu[i],append = F)
  write.FASTA(read.FASTA(fasta.out.omicron[i],type = "DNA"), fasta.out.omicron[i],append = F)
  write.FASTA(read.FASTA(fasta.out.zeta[i],type = "DNA"), fasta.out.zeta[i],append = F)
}

#' 
#' 
#' ## check meta look like approp lineages
## ----------------------------------------------------------------------------------------------------------
  # #sample VOC and boots one by one to not exceed memory
  # print("Sampling and writing for Alpha")
  # for (i in 1:b){
  #   samp.full.alpha<-read.csv(meta.out.alpha[i])
  #   print(table(samp.full.alpha$Lineage))
  # }
  # rm(samp.full.alpha)
  # 
  # #sample VOC and boots one by one to not exceed memory
  # print("Sampling and writing for beta")
  # for (i in 1:b){
  #   samp.full.beta<-read.csv(meta.out.beta[i])
  #   print(table(samp.full.beta$Lineage))
  # }
  # rm(samp.full.beta)
  # 
  # 
  # print("Sampling and writing for Gamma")
  # for (i in 1:b){
  #   samp.full.gamma<-read.csv(meta.out.gamma[i])
  #   print(table(samp.full.gamma$Lineage))
  # }
  # rm(samp.full.gamma)
  # 
  # print("Sampling and writing for Delta")
  # for (i in 1:b){
  #   samp.full.delta<-read.csv(meta.out.delta[i])    
  #   print(table(samp.full.delta$Lineage))
  # 
  # }
  # rm(samp.full.delta)
  # 
  # print("Sampling and writing for Omicron")
  # for (i in 1:b){
  #   samp.full.omicron<-read.csv(meta.out.omicron[i])
  #     print(table(samp.full.omicron$Lineage))
  # 
  # }
  # rm(samp.full.omicron)
  # 

 # print("Sampling and writing for other")
 # for (i in 1:b){
 #   samp.full.other<-read.csv(meta.out.other[i])
 #     print(table(samp.full.other$Lineage))
 # 
 # }
 # rm(samp.full.other)


#' 
#' 
