library(rlang)
library(ggplot2)
library(dplyr)
library(tidyr)

# Add MMETSP names
mmetspid <- read.csv("mmetsp_all_lineages.csv",header=F,stringsAsFactors = F)

## DSYB
# Read in table keep columns of interest
DSYB_reblast<-read.delim(file="DSYB_reblast.txt",header=F,stringsAsFactors = F)[,c(1:3,11)]
colnames(DSYB_reblast)<-c("qseqid","sseqid","pident","evalue")
DSYB_reblast_taxa<-read.csv("DSYB_reblasttaxa.csv",stringsAsFactors = F)
if(all(DSYB_reblast_taxa$accessid == DSYB_reblast$sseqid)){
  DSYB_reblast$kingdom <- DSYB_reblast_taxa$superkingdom
  #Some of the protein ID's are not currently available in the taxonomzir package, I manually checked that these are all Bacteria
  DSYB_reblast$kingdom[which(is.na(DSYB_reblast$kingdom)==T)] <- "Bacteria" }
for (i in 1: dim(DSYB_reblast)[1]){
  DSYB_reblast$isolate[i] <- mmetspid[grep(DSYB_reblast$qseqid[i], mmetspid[,1]),2]}
# Simplify to just euk or prok
DSYB_reblast$kingdom[which(DSYB_reblast$kingdom == "Bacteria" | DSYB_reblast$kingdom == "Archaea")] <- "Prokaryota"
# Filter by evalue
DSYB_reblast<-DSYB_reblast[which(DSYB_reblast$evalue<1e-30),]
# Create df that sums the euks, proks for each unique search id
DSYB_summary <- DSYB_reblast %>%
  group_by(qseqid, isolate, kingdom) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = kingdom, values_from = count) %>%
  mutate(Eukaryota = replace_na(Eukaryota, 0)) %>%
  mutate(Prokaryota = replace_na(Prokaryota, 0)) %>%
  mutate(perPro = Prokaryota/(Eukaryota + Prokaryota))
# Add median per id
tmp <- DSYB_reblast %>%
  group_by(qseqid) %>%
  summarise(m = median(pident))
DSYB_summary$med <- tmp$m
# Remove hits that with perpro < 97%
DSYB_contam<-DSYB_summary[which(DSYB_summary$perPro >=0.97),]
ggplot(DSYB_reblast,aes(x=pident,fill=kingdom))+geom_density()+ggtitle("DSYB")+xlim(c(20,100))
range(DSYB_reblast$pident[which(DSYB_reblast$kingdom=="Eukaryota")])
range(DSYB_reblast$pident[which(DSYB_reblast$kingdom=="Prokaryota")])


## TpMT2
# Read in table keep columns of interest
TpMT2_reblast<-read.delim(file="TpMT2_reblast.txt",header=F,stringsAsFactors = F)
colnames(TpMT2_reblast)<-c("qseqid","sseqid","pident","evalue")
TpMT2_reblast_taxa<-read.csv("TpMT2_reblasttaxa.csv",stringsAsFactors = F)
if(all(TpMT2_reblast_taxa$accessid == TpMT2_reblast$sseqid)){
  TpMT2_reblast$kingdom <- TpMT2_reblast_taxa$superkingdom
  #Some of the protein ID's are not currently available in the taxonomzir package, I manually checked that these are all Bacteria
  TpMT2_reblast$kingdom[which(is.na(TpMT2_reblast$kingdom)==T)] <- "Bacteria"}
for (i in 1: dim(TpMT2_reblast)[1]){
  TpMT2_reblast$isolate[i] <- mmetspid[grep(TpMT2_reblast$qseqid[i], mmetspid[,1]),2]}
# Simplify to just euk or prok
TpMT2_reblast$kingdom[which(TpMT2_reblast$kingdom == "Bacteria" | TpMT2_reblast$kingdom == "Archaea")] <- "Prokaryota"
# Filter by evalue
TpMT2_reblast<-TpMT2_reblast[which(TpMT2_reblast$evalue<1e-30),]
# Create df that sums the euks, proks for each unique search id
TpMT2_summary <- TpMT2_reblast %>%
  group_by(qseqid, isolate, kingdom) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = kingdom, values_from = count) %>%
  mutate(Eukaryota = replace_na(Eukaryota, 0)) %>%
  mutate(Prokaryota = replace_na(Prokaryota, 0)) %>%
  mutate(perPro = Prokaryota/(Eukaryota + Prokaryota))
# Add median per id
tmp <- TpMT2_reblast %>%
  group_by(qseqid) %>%
  summarise(m = median(pident))
TpMT2_summary$med <- tmp$m

# Remove hits that with perpro < 97%
TpMT2_contam<-TpMT2_summary[which(TpMT2_summary$perPro >=0.97),]
ggplot(TpMT2_reblast,aes(x=pident,fill=kingdom))+geom_density()+ggtitle("TpMT2")+xlim(c(20,100))
range(TpMT2_reblast$pident[which(TpMT2_reblast$kingdom=="Eukaryota")])
range(TpMT2_reblast$pident[which(TpMT2_reblast$kingdom=="Prokaryota")])
