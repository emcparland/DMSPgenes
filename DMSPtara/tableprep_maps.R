

setwd("~/DMSPgenes/DMSPtara/")
make_maps<-0
#Load required libraries
library(reshape2)
library(ggplot2)
#library(pheatmap)
library(RColorBrewer)
library(vsn)
library(hexbin)
library(varhandle)
library(maps)
library(plyr)
library(gridExtra)
library(rlang)

## load tables
gene <- "DSYB"
if (gene == "DSYB"){
  env_csv <- "environmental_parameters_DSYB.csv"
  aln_tsv <- "OceanGeneAtlas_20201123-22_51_DSYB_aln_result.tsv"
  abund_csv <- "abundance_matrix_DSYB.csv"}
if (gene =="TpMT2"){
  env_csv <- "environmental_parameters_TpMT2.csv"
  aln_tsv <- "OceanGeneAtlas_20201123-22_55_TpMT2_aln_result.tsv"
  abund_csv <- "abundance_matrix_TpMT2.csv"}

# create env data 
env<-read.delim(file = env_csv, sep="\t", header=TRUE,skip=0)
# create aln table
aln<-read.delim(file = aln_tsv,sep="\t",header = TRUE,skip=8)
# create abund table
abund <- read.delim(file = abund_csv, sep="\t",header=TRUE,skip=1)
## Create taxa tab
taxa <- abund[,1:2]
abund[,1:2] <- NULL

## filter results
# remove the extra number to make subject id searchable
tmp <- unlist(strsplit(as.character(aln$subject_id), split = "_"))
tmp <- tmp[which(tmp == "MATOU-v1")+1]
aln$subject_id <- paste0("MATOU-v1_",tmp)

# keep only hits with evalue < 1e-20
idx <- which(aln$e_value > 1e-30)

for (i in 1: length(idx)){
  j <- which(as.character(taxa$Gene_ID) %in% aln$subject_id[idx[i]])
  if (!is_empty(j)){abund <- abund[-j,]; taxa <- taxa[-j,]}
}

aln <- aln[-(idx),]

## Filter by depth
env <- env[which(env$depth == "SRF"),]
abund <- abund[,which(colnames(abund) %in% env$sample_ID)]

## Filter out fraction >180
env <- env[which(!env$size.fraction == "[180-2000µm]"),]
abund <- abund[,which(colnames(abund) %in% env$sample_ID)]

#env <- env[which(!env$size.fraction == "[20-180µm]"),]
#abund <- abund[,which(colnames(abund) %in% env$sample_ID)]

## Remove metazoa
idx <- grep("Metazoa", taxa$taxonomy)
if(is_empty(idx)==F){taxa <- taxa[-c(idx),]
abund <- abund[-c(idx),]}

## Remove bacteria
idx <- grep("Bacteria", taxa$taxonomy)
if(is_empty(idx)==F){taxa <- taxa[-c(idx),]
abund <- abund[-c(idx),]}

## Drop levels
env$depth <- droplevels(env$depth)
env$size.fraction <- droplevels(env$size.fraction)
all(env$sample_ID == colnames(abund))

## Sum the size fractions together at each station (=0.8-180um)
locals<-unique(sort(env$station))

a <- which(env$station == locals[1])
b <- which(colnames(abund) == env$sample_ID[a])

new_env <- as.data.frame(env[a[1],])
rownames(new_env) <- env$station[a[1]]

new_abund <- as.data.frame(rowSums(abund[,b]))
colnames(new_abund) <- env$station[a[1]]


for (i in 2:length(locals)){
  a <- which(env$station == locals[i])
  b <- which(colnames(abund) %in% env$sample_ID[a])
  
  if (length(b)==1) {
    new_env[i,] <- env[a,]
    rownames(new_env)[i] <- as.character(env$station[a])
    
    new_abund[,i] <- abund[,b]
    colnames(new_abund)[i] <- as.character(env$station[a])}
  
  if (length(b) > 1){
    new_env[i,] <- env[a[1],]
    rownames(new_env)[i] <- as.character(env$station[a[1]])
    
    new_abund[,i] <- rowSums(abund[,b])
    colnames(new_abund)[i] <- as.character(env$station[a[1]])}
}
rm(env,abund); abund <- new_abund
env <- new_env; rm(new_abund,new_env)
all(colnames(abund)==rownames(env))

if (gene == "DSYB"){
  DSYB_abund <- abund
  DSYB_env <- env
  DSYB_taxa <- taxa
}
if (gene == "TpMT2"){
  TpMT2_abund <- abund
  TpMT2_env <- env
  TpMT2_taxa <- taxa
}

rm(abund,aln,env,taxa)

################################################################################
### combine metaT DSYB and TpMT2 results for plotting
# Create a df of stations common (or unique to) each gene
stn_tot_abund <- as.data.frame(unique(sort(c(colnames(DSYB_abund),colnames(TpMT2_abund)))))
colnames(stn_tot_abund) <- "station"
stn_tot_abund$DSYB <- 0
stn_tot_abund$DSYB_hapto <- 0
stn_tot_abund$DSYB_chloro <- 0
stn_tot_abund$DSYB_pelago <- 0
stn_tot_abund$DSYB_bac <- 0
stn_tot_abund$DSYB_dino <- 0
stn_tot_abund$DSYB_other <- 0

stn_tot_abund$TpMT2 <- 0
stn_tot_abund$TpMT2_hapto <- 0
stn_tot_abund$TpMT2_chloro <- 0
stn_tot_abund$TpMT2_pelago <- 0
stn_tot_abund$TpMT2_bac <- 0
stn_tot_abund$TpMT2_dino <- 0
stn_tot_abund$TpMT2_other <-0

# index for DSYBgroups
hidx <- grep("Haptophyceae",DSYB_taxa$taxonomy)
cidx<-grep("Viridiplantae",DSYB_taxa$taxonomy)
pidx<-grep("Pelagophyceae",DSYB_taxa$taxonomy)
bidx<-grep("Bacillariophyta",DSYB_taxa$taxonomy)
didx<-grep("Dinophyceae",DSYB_taxa$taxonomy)
idx_used<-c(hidx,cidx,pidx,bidx,didx)
idx<-(!c(1:dim(DSYB_taxa)[1]) %in% idx_used)
idx_Other<-which(idx==TRUE)

#index for TpMT2 groups
hidx_tpmt <- grep("Haptophyceae",TpMT2_taxa$taxonomy)
cidx_tpmt<-grep("Viridiplantae",TpMT2_taxa$taxonomy)
pidx_tpmt<-grep("Pelagophyceae",TpMT2_taxa$taxonomy)
bidx_tpmt<-grep("Bacillariophyta",TpMT2_taxa$taxonomy)
didx_tpmt<-grep("Dinophyceae",TpMT2_taxa$taxonomy)
idx_used<-c(hidx_tpmt,cidx_tpmt,pidx_tpmt,bidx_tpmt,didx_tpmt)
idx<-(!c(1:dim(TpMT2_taxa)[1]) %in% idx_used)
idx_Other_tpmt<-which(idx==TRUE)

for (i in 1:dim(stn_tot_abund)[1]){
  # 
  idx <- which(colnames(DSYB_abund) == stn_tot_abund$station[i])
  stn_tot_abund$DSYB[i] <- sum(DSYB_abund[,idx])
  stn_tot_abund$DSYB_hapto[i] <- sum(DSYB_abund[hidx,idx])
  stn_tot_abund$DSYB_chloro[i] <- sum(DSYB_abund[cidx,idx])
  stn_tot_abund$DSYB_pelago[i] <- sum(DSYB_abund[pidx,idx])
  stn_tot_abund$DSYB_bac[i] <- sum(DSYB_abund[bidx,idx])
  stn_tot_abund$DSYB_dino[i] <- sum(DSYB_abund[didx,idx])
  stn_tot_abund$DSYB_other[i] <- sum(DSYB_abund[idx_Other,idx])
  
  idx <- which(colnames(TpMT2_abund) == stn_tot_abund$station[i])
  stn_tot_abund$TpMT2[i] <- sum(TpMT2_abund[,idx])
  stn_tot_abund$TpMT2_hapto[i] <- sum(TpMT2_abund[hidx_tpmt,idx])
  stn_tot_abund$TpMT2_chloro[i] <- sum(TpMT2_abund[cidx_tpmt,idx])
  stn_tot_abund$TpMT2_pelago[i] <- sum(TpMT2_abund[pidx_tpmt,idx])
  stn_tot_abund$TpMT2_bac[i] <- sum(TpMT2_abund[bidx_tpmt,idx])
  stn_tot_abund$TpMT2_dino[i] <- sum(TpMT2_abund[didx_tpmt,idx])
  stn_tot_abund$TpMT2_other[i] <- sum(TpMT2_abund[idx_Other_tpmt,idx])
}
# sum of tot DSYB + tot TpMT2
idx <- which(colnames(stn_tot_abund) == "DSYB" | colnames(stn_tot_abund) == "TpMT2")
stn_tot_abund$Total <- rowSums(stn_tot_abund[,idx])

# Make an env df that will match the stn tot abund df
tot_env <- as.data.frame(stn_tot_abund$station)
colnames(tot_env) <- "station"

# combine the two env df and then break down to one that matches tot abund df
tmp <- as.data.frame(rbind(DSYB_env,TpMT2_env))
idx <- which(stn_tot_abund$station[1] == tmp$station)
tot_env <- tmp[idx[1],]
for (i in 2: dim(stn_tot_abund)[1]){
  idx <- which(stn_tot_abund$station[i] == tmp$station)
  tot_env[i,] <- tmp[idx[1],]
}

################################################################################
### Map making ###
zoom=FALSE
LABELS <- 0 # set to 0 for no labels

# Map: DYSB+TpMT2 @ all stns, color = DYSB/total
df<-data.frame(stn_tot_abund$Total)
df$Color_all<-stn_tot_abund$DSYB/stn_tot_abund$Total
df$Lat<-tot_env$latitude; df$Lon<-tot_env$longitude;df$samples<-tot_env$station
colnames(df)[1:2]<-c("Class_all","Color_all")

ASSIGN_MAP_TITLE="DYSB+TpMT2 All"
CIRCLE_LOWER_LIM=0 #min(df$Class_all)
CIRCLE_UPPER_LIM=max(df$Class_all)
source("map_plotting.R")

# Map: Dino Totals
df<-data.frame(stn_tot_abund$DSYB_dino+stn_tot_abund$TpMT2_dino)
df$Color_all<-stn_tot_abund$DSYB_dino/(stn_tot_abund$DSYB_dino+stn_tot_abund$TpMT2_dino)
df$Lat<-tot_env$latitude; df$Lon<-tot_env$longitude;df$samples<-tot_env$station
colnames(df)[1:2]<-c("Class_all","Color_all")

ASSIGN_MAP_TITLE="Dino"
source("map_plotting.R")

# Map: Hapto Totals
df<-data.frame(stn_tot_abund$DSYB_hapto+stn_tot_abund$TpMT2_hapto)
df$Color_all<-stn_tot_abund$DSYB_hapto/(stn_tot_abund$DSYB_hapto+stn_tot_abund$TpMT2_hapto)
df$Lat<-tot_env$latitude; df$Lon<-tot_env$longitude;df$samples<-tot_env$station
colnames(df)[1:2]<-c("Class_all","Color_all")
ASSIGN_MAP_TITLE="Hapto"
source("map_plotting.R")

# Map: Bacill Totals
df<-data.frame(stn_tot_abund$DSYB_bac+stn_tot_abund$TpMT2_bac)
df$Color_all<-stn_tot_abund$DSYB_bac/(stn_tot_abund$DSYB_bac+stn_tot_abund$TpMT2_bac)
df$Lat<-tot_env$latitude; df$Lon<-tot_env$longitude;df$samples<-tot_env$station
colnames(df)[1:2]<-c("Class_all","Color_all")
ASSIGN_MAP_TITLE="Bacillariophyta"
source("map_plotting.R")

# Map: Pelago/Chloroph Totals 
df<-data.frame(stn_tot_abund$DSYB_pelago+stn_tot_abund$TpMT2_pelago+
                 stn_tot_abund$DSYB_chloro+stn_tot_abund$TpMT2_chloro)
df$Color_all<-(stn_tot_abund$DSYB_pelago+stn_tot_abund$DSYB_chloro)/(stn_tot_abund$DSYB_pelago+stn_tot_abund$TpMT2_pelago+
                                           stn_tot_abund$DSYB_chloro+stn_tot_abund$TpMT2_chloro)
df$Lat<-tot_env$latitude; df$Lon<-tot_env$longitude;df$samples<-tot_env$station
colnames(df)[1:2]<-c("Class_all","Color_all")
ASSIGN_MAP_TITLE="Chlorophyta/Pelagophyta"
source("map_plotting.R")

######
CIRCLE_LOWER_LIM=0 
CIRCLE_UPPER_LIM=1

# Map: Dino Percent
df<-data.frame((stn_tot_abund$DSYB_dino+stn_tot_abund$TpMT2_dino)/stn_tot_abund$Total)
df$Color_all<-stn_tot_abund$DSYB_dino/(stn_tot_abund$DSYB_dino+stn_tot_abund$TpMT2_dino)
df$Lat<-tot_env$latitude; df$Lon<-tot_env$longitude;df$samples<-tot_env$station
colnames(df)[1:2]<-c("Class_all","Color_all")

ASSIGN_MAP_TITLE="Dino"
source("map_plotting.R")

# Map: Hapto Percent
df<-data.frame((stn_tot_abund$DSYB_hapto+stn_tot_abund$TpMT2_hapto)/stn_tot_abund$Total)
df$Color_all<-stn_tot_abund$DSYB_hapto/(stn_tot_abund$DSYB_hapto+stn_tot_abund$TpMT2_hapto)
df$Lat<-tot_env$latitude; df$Lon<-tot_env$longitude;df$samples<-tot_env$station
colnames(df)[1:2]<-c("Class_all","Color_all")
ASSIGN_MAP_TITLE="Hapto"
source("map_plotting.R")

# Map: Bacill Percent
df<-data.frame((stn_tot_abund$DSYB_bac+stn_tot_abund$TpMT2_bac)/stn_tot_abund$Total)
df$Color_all<-stn_tot_abund$DSYB_bac/(stn_tot_abund$DSYB_bac+stn_tot_abund$TpMT2_bac)
df$Lat<-tot_env$latitude; df$Lon<-tot_env$longitude;df$samples<-tot_env$station
colnames(df)[1:2]<-c("Class_all","Color_all")
ASSIGN_MAP_TITLE="Bacillariophyta"
source("map_plotting.R")

# Map: Pelago/Chloroph Percent
df<-data.frame((stn_tot_abund$DSYB_pelago+stn_tot_abund$TpMT2_pelago+
                 stn_tot_abund$DSYB_chloro+stn_tot_abund$TpMT2_chloro)/stn_tot_abund$Total)
df$Color_all<-(stn_tot_abund$DSYB_pelago+stn_tot_abund$DSYB_chloro)/(stn_tot_abund$DSYB_pelago+stn_tot_abund$TpMT2_pelago+
                                                                       stn_tot_abund$DSYB_chloro+stn_tot_abund$TpMT2_chloro)
df$Lat<-tot_env$latitude; df$Lon<-tot_env$longitude;df$samples<-tot_env$station
colnames(df)[1:2]<-c("Class_all","Color_all")
ASSIGN_MAP_TITLE="Chlorophyta/Pelagophyta"
source("map_plotting.R")

# Map: Tot chla
df<-data.frame(tot_env$Chlorophyll_A..mg.m..3.)
df$Color_all<-0.5
df$Lat<-tot_env$latitude; df$Lon<-tot_env$longitude;df$samples<-tot_env$station
colnames(df)[1:2]<-c("Class_all","Color_all")

CIRCLE_LOWER_LIM=0.001 #min(df$Class_all)
CIRCLE_UPPER_LIM=5 #max(df$Class_all)
ASSIGN_MAP_TITLE="Chla (mg/m3)"
source("map_plotting.R")


# Map: Tot 
df<-data.frame(stn_tot_abund$Total)
df$Color_all<-1
df$Lat<-tot_env$latitude; df$Lon<-tot_env$longitude;df$samples<-tot_env$station
colnames(df)[1:2]<-c("Class_all","Color_all")

CIRCLE_LOWER_LIM=0 #min(df$Class_all)
CIRCLE_UPPER_LIM=max(df$Class_all)
ASSIGN_MAP_TITLE="DSYB"
source("map_plotting.R")





##########################
# Looking at how Chloro/Pelago and Diatoms switch
df <- data.frame((stn_tot_abund$TpMT2_chloro+stn_tot_abund$TpMT2_pelago)/stn_tot_abund$TpMT2)
colnames(df) <- "perCP"
df$perbac <- stn_tot_abund$TpMT2_bac/stn_tot_abund$TpMT2
df$chla <- tot_env$NO3..µmol.l.

ggplot(df,aes(x=perCP,y=perbac,color=chla))+geom_point()+
  scale_color_gradientn(colors=rainbow(5),limits=c(0,0.5))

