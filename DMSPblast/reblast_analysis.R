library(rlang)
library(ggplot2)
accession_search<-read.delim(file="search2.txt",header=F)
colnames(accession_search)<-c("accession","speciesname")

kingdom_search<-read.delim(file="search1.txt",header=F)
colnames(kingdom_search)<-c("speciesname","kingdom")
kingdom_search$kingdom<-as.character(kingdom_search$kingdom)
idx<-grep("bacteria|archae|firmicutes|high GC Gram+|planctomycetes",kingdom_search$kingdom)
kingdom_search$kingdom[idx]<-"Prokaryote"
idx<-grep("ants|aphids|ascomycetes|bats|basidiomycetes|beetles|bony fishes|choanoflagellates|crustaceans|cryptomonads|diatoms|dinoflagellates|eudicots|eukaryotes|flies|frogs & toads|gastropods|green algae|haptophytes|hemichordates|hydrozoans|lizards|moths|pelagophytes|priapulids|red algae|rodents|sea anemones|segmented worms|snakes|soft corals|species|springtails|stony corals|termites|wasps",kingdom_search$kingdom)
kingdom_search$kingdom[idx]<-"Eukaryote"

DSYB_reblast_all<-read.delim(file="DSYB_reblast_all.txt",header=F,sep = " ")
idx<-c(1,3,5,13)
DSYB_reblast_all<-DSYB_reblast_all[,idx]
colnames(DSYB_reblast_all)<-c("qseqid","sseqid","pident","evalue")
DSYB_reblast_all$kingdom<-NA
DSYB_reblast_all<-DSYB_reblast_all[which(DSYB_reblast_all$evalue<1e-30),]

for (i in 1: dim(DSYB_reblast_all)[1]){
  accid<-unlist(strsplit(as.character(DSYB_reblast_all$sseqid[i]),"\\|"))[4]
  species<-accession_search$speciesname[grep(accid,accession_search$accession)]
  
  kingdom<-kingdom_search$kingdom[grep(species,kingdom_search$speciesname)]
  
  if(!is_empty(kingdom)){DSYB_reblast_all$kingdom[i]<-kingdom}
  if(length(kingdom)>1){print(c(i,kingdom))}
}

DSYB_uniq<-as.data.frame(unique(sort(DSYB_reblast_all$qseqid)))
colnames(DSYB_uniq)<-"MMETSPname"
DSYB_uniq$eukno<-NA
DSYB_uniq$prokno<-NA
DSYB_uniq$perpro<-NA
for (i in 1:dim(DSYB_uniq)[1]){
  tmp<-DSYB_reblast_all$kingdom[grep(DSYB_uniq$MMETSPname[i],DSYB_reblast_all$qseqid)]
  DSYB_uniq$eukno[i]<-as.numeric(length(which(tmp=="Eukaryote")))
  DSYB_uniq$prokno[i]<-as.numeric(length(which(tmp=="Prokaryote")))
  DSYB_uniq$perpro[i]<-(DSYB_uniq$prokno[i]/(DSYB_uniq$eukno[i]+DSYB_uniq$prokno[i]))*100
  
}  

# Remove hits that with perpro > 95%
contam<-as.character(DSYB_uniq$MMETSPname[which(DSYB_uniq$perpro>97)])
idx<-which(DSYB_reblast_all$qseqid==contam[1]|DSYB_reblast_all$qseqid==contam[2]|
             DSYB_reblast_all$qseqid==contam[3]|DSYB_reblast_all$qseqid==contam[4])
DSYB_reblast_all<-DSYB_reblast_all[-idx,]
ggplot(DSYB_reblast_all,aes(x=pident,fill=kingdom))+geom_density()+ggtitle("DSYB")




TpMT2_reblast_all<-read.delim(file="TpMT2_reblast_all.txt",header=F)
idx<-c(1,3,5,13)
TpMT2_reblast_all<-TpMT2_reblast_all[,idx]
colnames(TpMT2_reblast_all)<-c("qseqid","sseqid","pident","evalue")
TpMT2_reblast_all$kingdom<-NA
TpMT2_reblast_all<-TpMT2_reblast_all[which(TpMT2_reblast_all$evalue<1e-30),]

for (i in 1: dim(TpMT2_reblast_all)[1]){
  accid<-unlist(strsplit(as.character(TpMT2_reblast_all$sseqid[i]),"\\|"))[4]
  species<-accession_search$speciesname[grep(accid,accession_search$accession)]
  
  kingdom<-kingdom_search$kingdom[grep(species,kingdom_search$speciesname)]
  
  if(!is_empty(kingdom)){TpMT2_reblast_all$kingdom[i]<-kingdom}
  if(length(kingdom)>1){print(c(i,kingdom))}
}

TpMT2_uniq<-as.data.frame(unique(sort(TpMT2_reblast_all$qseqid)))
colnames(TpMT2_uniq)<-"MMETSPname"
TpMT2_uniq$eukno<-NA
TpMT2_uniq$prokno<-NA
TpMT2_uniq$perpro<-NA
for (i in 1:dim(TpMT2_uniq)[1]){
  tmp<-TpMT2_reblast_all$kingdom[grep(TpMT2_uniq$MMETSPname[i],TpMT2_reblast_all$qseqid)]
  TpMT2_uniq$eukno[i]<-as.numeric(length(which(tmp=="Eukaryote")))
  TpMT2_uniq$prokno[i]<-as.numeric(length(which(tmp=="Prokaryote")))
  TpMT2_uniq$perpro[i]<-(TpMT2_uniq$prokno[i]/(TpMT2_uniq$eukno[i]+TpMT2_uniq$prokno[i]))*100
  
}  

# Remove hits that with perpro > 97%
#contam<-as.character(TpMT2_uniq$MMETSPname[which(TpMT2_uniq$perpro>95)])
#idx<-which(TpMT2_reblast_all$qseqid==contam[1]|TpMT2_reblast_all$qseqid==contam[2]|
#             TpMT2_reblast_all$qseqid==contam[3]|TpMT2_reblast_all$qseqid==contam[4])
#TpMT2_reblast_all<-TpMT2_reblast_all[-idx,]
ggplot(TpMT2_reblast_all,aes(x=pident,fill=kingdom))+geom_density()+ggtitle("TpMT2")
