
library(ggplot2)
library(cowplot)
library(dplyr)
library(reshape2)
library(rlang)
library(Hmisc)

setwd("~/DMSPgenes/DMSPblast/")

pheno_full<-read.csv("PhenotypeList.csv")
colnames(pheno_full)<-c("DMSPType","Kingdom","Group","Phylum","Class","Species",
                        "Strain","IntraDMSP","Genus","Generic")

# Create presence or absence table
pheno_full$HiDP=0; pheno_full$LoDP=0; pheno_full$det=0; pheno_full$bd=0;
for (i in 1:length(pheno_full$DMSPType)){
  if (pheno_full$DMSPType[i] == "High") {pheno_full$HiDP[i] = 1}
  if (pheno_full$DMSPType[i] == "Low") {pheno_full$LoDP[i] = 1}
  if (pheno_full$DMSPType[i] == "Det") {pheno_full$det[i] = 1}
  if (pheno_full$DMSPType[i] == "BD") {pheno_full$bd[i] = 1}
}
#pheno_full$Generic[which(pheno_full$Generic=="Prokaryote")]<-"Prokaryotes"

pheno_table<- pheno_full %>%
  select(Genus,HiDP,LoDP,det,bd) %>%
  group_by(Genus) %>%
  summarise_each(funs(sum)) %>%
  mutate(tot=HiDP+LoDP+det+bd)

pheno_table_sum <- pheno_full  %>%
  select(Generic,HiDP,LoDP,det,bd) %>%
  group_by(Generic) %>%
  summarise_each(funs(sum))

pheno_forplot<-pheno_table %>%
  select(HiDP,LoDP,det,bd,Genus) %>%
  #mutate(det_bd=det+bd) %>%
  select(HiDP,LoDP,det,bd,Genus) %>%
  mutate_at(vars(-Genus),funs(./pheno_table$tot)) %>%
  melt(id="Genus")
pheno_forplot$Genus<-factor(pheno_forplot$Genus,levels=c("Alphaproteobacteria","Gammaproteobacteria",
                                     "Prochlorococcus","Synechococcus","Trichodesmium","Other Cyanobacteria","Other Prokaryotes",
                                     "Alexandrium","Heterocapsa",
                                     "Karenia","Prorocentrum","Scrippsiella","Symbiodinium","Other Dinophyta",
                                     "Calcidiscus","Chrysochromulina","Coccolithus","Emiliania","Gephyrocapsa","Imantonia",
                                     "Isochrysis","Pavlova","Phaeocystis","Prymnesium","Other Haptophyta",
                                     "Bathycoccus","Mantoniella","Micromonas","Pycnococcus","Pyramimonas",
                                     "Tetraselmis","Other Chlorophyta",
                                     "Chaetoceros","Coscinodiscus","Cylindrotheca","Fragilariopsis",
                                     "Pseudo-nitzschia","Skeletonema","Thalassiosira","Other Bacillariophyta",
                                     "Ochromonas","Pelagococcus","Pelagomonas","Other Pelagophyta","Other Stramenopiles",
                                     "Cryptomonas","Cryptophyta",
                                     "Bigelowiella","Rhizaria","Rhodophyta","Excavata",
                                     "Glaucophyta","Euglenozoa","Ciliophora","Other"))

bub1<-
  ggplot(pheno_forplot,aes(x=variable,y=Genus,size=ifelse(value==0, NA, value),color=variable))+
  geom_point(alpha=1)+
  guides(color=FALSE)+
  scale_size(range=c(0.0001, 4),breaks=c(0.25,0.5,1.0),name="% of Total")+
  scale_x_discrete(position="top",labels=c("HiDP","LoDP","Detect","Det+BD"))+
  scale_y_discrete(limits=rev(levels(pheno_forplot$Genus)))+
  #xlab("Phenotype")+
  scale_color_manual(values=c("red2","royalblue","black","grey"))+
  theme_minimal()+
  theme(axis.text.y = element_blank(),
        text = element_text(size=14),
        axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=12),
        legend.position = "none",
        panel.background = element_blank(),
        plot.margin = margin(0,0,0,0,"cm"))

pheno_table$Genus<-factor(pheno_table$Genus,levels=c("Alphaproteobacteria","Gammaproteobacteria",
                                                     "Prochlorococcus","Synechococcus","Trichodesmium","Other Cyanobacteria","Other Prokaryotes",
                                                     "Alexandrium","Heterocapsa",
                                                     "Karenia","Prorocentrum","Scrippsiella","Symbiodinium","Other Dinophyta",
                                                     "Calcidiscus","Chrysochromulina","Coccolithus","Emiliania","Gephyrocapsa","Imantonia",
                                                     "Isochrysis","Pavlova","Phaeocystis","Prymnesium","Other Haptophyta",
                                                     "Bathycoccus","Mantoniella","Micromonas","Pycnococcus","Pyramimonas",
                                                     "Tetraselmis","Other Chlorophyta",
                                                     "Chaetoceros","Coscinodiscus","Cylindrotheca","Fragilariopsis",
                                                     "Pseudo-nitzschia","Skeletonema","Thalassiosira","Other Bacillariophyta",
                                                     "Ochromonas","Pelagococcus","Pelagomonas","Other Pelagophyta","Other Stramenopiles",
                                                     "Cryptomonas","Cryptophyta",
                                                     "Bigelowiella","Rhizaria","Rhodophyta","Excavata",
                                                     "Glaucophyta","Euglenozoa","Ciliophora","Other"))


vbar1<-
  ggplot(pheno_table,aes(x=factor(Genus,levels=rev(levels(factor(Genus)))),y=tot))+
  geom_col(fill="slategrey",position="dodge")+
  scale_y_reverse(limits=c(25,0))+
  xlab("")+
  coord_flip()+
  theme_minimal()+
  theme(#axis.title = element_blank(),
        text = element_text(size=14),
        axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        panel.grid = element_blank(),
        plot.margin = margin(0,0,0,0,"cm"))

pheno_tot<-as.data.frame(pheno_table %>%
                          select("HiDP","LoDP","det","bd") %>%
                          mutate(det_bd=det+bd) %>%
                          summarise_each(funs(sum)) %>%
                          select("HiDP","LoDP","det","det_bd") %>%
                          t())
                          
colnames(pheno_tot)<-"variable"
pheno_tot$Names<-rownames(pheno_tot)
pheno_tot$Names<-factor(pheno_tot$Names,levels=c("HiDP","LoDP","det","det_bd"))
hbar1<-
  ggplot(pheno_tot,aes(x=Names,y=variable,fill=Names))+
  geom_col(position="dodge",width=0.25)+
  labs(title="Phenotype")+
  ylim(c(0,160))+
  theme_minimal()+
  scale_fill_manual(values=c("red2","royalblue","black","grey"))+
  theme(axis.title=element_blank(),
        text = element_text(size=14),
        axis.text.x=element_blank(),
        axis.line.y = element_line(),
        axis.ticks.y = element_line(),
        legend.position = "none",
        panel.grid = element_blank(),
        plot.margin = margin(0,0,0,0,"cm"))



geno_full<-read.delim("genekey_final.txt",header = TRUE)
mmetsp_line<-read.csv("mmetsp_all_lineages.csv",header=FALSE)
colnames(mmetsp_line)<-c("mmetspid","strain","","genus","family","order","class",
                         "rank","rank","rank","rank","kingdom","Generic1","Specific","Generic")
prok_line<-read.delim("prok_line.txt",header=TRUE)
geno_full$order<-"NA"
geno_full$generic<-"NA"
for(i in 1:dim(geno_full)[1]){
  idx<-grep(geno_full$MMETSPstrain[i],mmetsp_line$strain)[1]
  
  if (!is.na(idx)){
  geno_full$order[i]<-as.character(mmetsp_line$Specific[idx])
  geno_full$generic[i]<-as.character(mmetsp_line$Generic[idx])}
  
  if(is.na(idx)){
    tmp<-as.character(lapply(strsplit(as.character(geno_full$MMETSPstrain[i]), split="-"), "[",1))
    idx<-grep(tmp,mmetsp_line$order)[1]
    if(!is.na(idx)){
    geno_full$order[i]<-as.character(mmetsp_line$Specific[idx])
    geno_full$generic[i]<-as.character(mmetsp_line$Generic[idx])
    }}
  
  if(is.na(idx)){
    idx<-grep(geno_full$MMETSPstrain[i],prok_line$Species)[1]
    geno_full$order[i]<-as.character(prok_line$Genus[idx])
    geno_full$generic[i]<-as.character(prok_line$Generic[idx])
    }
}
idx<- which(geno_full$order=="Actinobacteria"); geno_full$order[idx]<-"Other Prokaryotes"
idx<- which(geno_full$order=="Firmicutes"); geno_full$order[idx]<-"Other Prokaryotes"
idx<- which(geno_full$order=="Acidobacteria"); geno_full$order[idx]<-"Other Prokaryotes"
idx<- which(geno_full$order=="Deltaproteobacteria"); geno_full$order[idx]<-"Other Prokaryotes"
idx <- which(geno_full$MMETSPstrain=="Nostoc sp."); geno_full$order[idx]<-"Other Cyanobacteria"

# Add columns for no genes and DSYB+TpMT2 genotype
geno_full$TpMT2_DSYB<-0
geno_full$None<-0
for(i in 1:dim(geno_full)[1]){
  # Add column for No hit
  if(geno_full$DSYB[i]==0 & geno_full$TpMT2[i]==0){geno_full$None[i]=1}
  
  # Add column for TpMT2 + DSYB
  if(geno_full$DSYB[i] ==1 & geno_full$TpMT2[i]==1){geno_full$TpMT2_DSYB[i] = 1}
  if(geno_full$TpMT2_DSYB[i] ==1){geno_full[i,1:2]=0}
}

# Here I'm adding in numbers based on refseq prok totals
a<-c(0,0,"","Alphaproteobacteria","Prokaryote",0,78)
c<-c(0,0,"","Other Cyanobacteria","Prokaryote",0,11)
p<-c(0,0,"","Prochlorococcus","Prokaryote",0,7)
s<-c(0,0,"","Synechococcus","Prokaryote",0,0)
t<-c(0,0,"","Trichodesmium","Prokaryote",0,0)
g<-c(0,0,"","Gammaproteobacteria","Prokaryote",0,183)
o<-c(0,0,"","Other Prokaryotes","Prokaryote",0,873)

geno_full$MMETSPstrain<-as.character(geno_full$MMETSPstrain)
geno_full<- rbind(geno_full,a,c,p,s,t,g,o)
geno_full$DSYB<-as.numeric(geno_full$DSYB);geno_full$TpMT2<-as.numeric(geno_full$TpMT2)
geno_full$TpMT2_DSYB<-as.numeric(geno_full$TpMT2_DSYB);geno_full$None<-as.numeric(geno_full$None)

# Remove contaminants based on reblast-analysis
for (i in 1: dim(DSYB_contam)[1]){
  idx <- grep(DSYB_contam$isolate[i], geno_full$MMETSPstrain)
  print(idx)
  geno_full$DSYB[idx] <- 0
  geno_full$None[idx] <- 1
}

for (i in 1: dim(TpMT2_contam)[1]){
  idx <- grep(TpMT2_contam$isolate[i], geno_full$MMETSPstrain)
  print(idx)
  geno_full$TpMT2[idx] <- 0
  geno_full$None[idx] <- 1
}

geno_table<- geno_full %>%
  select(order,DSYB,TpMT2,TpMT2_DSYB,None) %>%
  group_by(order) %>%
  summarise_each(funs(sum)) %>%
  mutate(tot=DSYB+TpMT2+TpMT2_DSYB+None)

geno_table_sum <- geno_full %>%
  select(generic,DSYB,TpMT2,TpMT2_DSYB,None) %>%
  group_by(generic) %>%
  summarise_each(funs(sum)) %>%
  mutate(tot=DSYB+TpMT2+TpMT2_DSYB+None)

geno_forplot<-geno_table %>%
  select(DSYB,TpMT2,TpMT2_DSYB,None,order) %>%
  #mutate(TpMT2_tot=TpMT2+TpMT2_DSYB) %>%
  select(DSYB,TpMT2,TpMT2_DSYB,None,order) %>%
  mutate_at(vars(-order),funs(./geno_table$tot)) %>%
  melt(id="order")
geno_forplot$order<-factor(geno_forplot$order,levels=c("Alphaproteobacteria","Gammaproteobacteria",
                                                       "Prochlorococcus","Synechococcus","Trichodesmium","Other Cyanobacteria","Other Prokaryotes",
                                                       "Alexandrium","Heterocapsa",
                                                       "Karenia","Prorocentrum","Scrippsiella","Symbiodinium","Other Dinophyta",
                                                       "Calcidiscus","Chrysochromulina","Coccolithus","Emiliania","Gephyrocapsa","Imantonia",
                                                       "Isochrysis","Pavlova","Phaeocystis","Prymnesium","Other Haptophyta",
                                                       "Bathycoccus","Mantoniella","Micromonas","Pycnococcus","Pyramimonas",
                                                       "Tetraselmis","Other Chlorophyta",
                                                       "Chaetoceros","Coscinodiscus","Cylindrotheca","Fragilariopsis",
                                                       "Pseudo-nitzschia","Skeletonema","Thalassiosira","Other Bacillariophyta",
                                                       "Ochromonas","Pelagococcus","Pelagomonas","Other Pelagophyta","Other Stramenopiles",
                                                       "Cryptomonas","Cryptophyta",
                                                       "Bigelowiella","Rhizaria","Rhodophyta","Excavata",
                                                       "Glaucophyta","Euglenozoa","Ciliophora","Other"))

bub2<-
  ggplot(geno_forplot,aes(x=variable,y=order,size=ifelse(value==0, NA, value),color=variable))+
  geom_point(alpha=1)+
  guides(color=FALSE)+
  scale_size(range=c(0.0001, 4),breaks=c(0.25,0.5,1.0),name="% of Total")+
  scale_x_discrete(position="top")+
  scale_y_discrete(limits=rev(levels(geno_forplot$order)))+
  #xlab("Genotype")+
  ylab("order")+
  scale_color_manual(values=c("red2","royalblue","slateblue1","grey"))+
  theme_minimal()+
  theme(axis.text.y = element_blank(),
        text = element_text(size=14),
        axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=12),
        legend.position = "none",
        panel.background = element_blank(),
        plot.margin = margin(0,0,0,0,"cm"))


geno_table$order<-factor(geno_table$order,levels=c("Alphaproteobacteria","Gammaproteobacteria",
                                                   "Prochlorococcus","Synechococcus","Trichodesmium","Other Cyanobacteria","Other Prokaryotes",
                                                   "Alexandrium","Heterocapsa",
                                                   "Karenia","Prorocentrum","Scrippsiella","Symbiodinium","Other Dinophyta",
                                                   "Calcidiscus","Chrysochromulina","Coccolithus","Emiliania","Gephyrocapsa","Imantonia",
                                                   "Isochrysis","Pavlova","Phaeocystis","Prymnesium","Other Haptophyta",
                                                   "Bathycoccus","Mantoniella","Micromonas","Pycnococcus","Pyramimonas",
                                                   "Tetraselmis","Other Chlorophyta",
                                                   "Chaetoceros","Coscinodiscus","Cylindrotheca","Fragilariopsis",
                                                   "Pseudo-nitzschia","Skeletonema","Thalassiosira","Other Bacillariophyta",
                                                   "Ochromonas","Pelagococcus","Pelagomonas","Other Pelagophyta","Other Stramenopiles",
                                                   "Cryptomonas","Cryptophyta",
                                                   "Bigelowiella","Rhizaria","Rhodophyta","Excavata",
                                                   "Glaucophyta","Euglenozoa","Ciliophora","Other"))
  
  

vbar2<-
    ggplot(geno_table,aes(x=factor(order,levels=rev(levels(factor(order)))),y=tot))+
    geom_col(fill="slategrey",position="dodge")+
    #scale_y_discrete(limits=rev(levels(pheno_forplot$Genus)))+
    #scale_y_reverse()+
    coord_flip()+ylim(c(0,160))+
    xlab(" ")+
    theme_minimal()+
    theme(axis.title.y = element_blank(),
          axis.text.y=element_blank(),
          text = element_text(size=14),
          axis.line.x = element_line(),
          axis.ticks.x = element_line(),
          panel.grid = element_blank(),
          plot.margin = margin(0,0,0,0,"cm"))


geno_tot<-as.data.frame(geno_table %>%
  select("DSYB","TpMT2","TpMT2_DSYB","None") %>%
  mutate(TpMT2_tot=TpMT2+TpMT2_DSYB) %>%
  summarise_each(funs(sum)) %>%
  select("DSYB","TpMT2","TpMT2_tot","None") %>%
  t())
colnames(geno_tot)<-"variable"
geno_tot$Names<-rownames(geno_tot)
geno_tot$Names<-factor(geno_tot$Names,levels=c("DSYB","TpMT2","TpMT2_tot","None"))

hbar2<-
  ggplot(geno_tot,aes(x=Names,y=variable,fill=Names))+
  geom_col(position="dodge",width=0.25)+
  labs(title="Genotype")+
  ylim(c(0,160))+
  theme_minimal()+
  scale_fill_manual(values=c("red2","royalblue","slateblue1","grey"))+
  theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        legend.position = "none",
        axis.line.y = element_line(),
        text = element_text(size=14),
        axis.ticks.y = element_line(),
        panel.grid = element_blank(),
        plot.margin = margin(0,0,0,0,"cm"))

## Table 1
rows <- c("Prokaryote","Dinophyta","Haptophyta","Chlorophyta","Bacillariophyta",
          "Pelagophyta","Stramenopiles","Cryptophyta","Rhizaria","Other")

table1 <- bind_cols(pheno_table_sum,geno_table_sum[,2:6]) %>%
  # phenotype
  mutate(perHi=(HiDP/(HiDP+LoDP))*100) %>%
  mutate(perLo=(LoDP/(HiDP+LoDP))*100) %>%
  mutate(nPheno=HiDP+LoDP) %>%
  mutate(perdet=(det/(det+bd))*100) %>%
  mutate(perbd=(bd/(det+bd))*100) %>%
  mutate(ndet=det+bd) %>%
  mutate(totPheno=HiDP+LoDP+det+bd) %>%
  # genotype
  mutate(perDSYB=(DSYB/tot)*100) %>%
  mutate(perTpMT2=(TpMT2/tot)*100) %>%
  mutate(perTpMT2_DSYB=(TpMT2_DSYB/tot)*100) %>%
  mutate(perNone=(None/tot)*100) %>%
  mutate(totGeno=tot) %>%
  select(Generic,totPheno,ndet,perdet,perbd,nPheno,perHi,perLo,totGeno,perDSYB,perTpMT2,perTpMT2_DSYB,perNone) %>%
  slice(match(rows, Generic))
  
write.csv(x=table1,"Table1.csv")




## Correlation of pheno and geno
pgcorr<-as.data.frame(geno_full$MMETSPstrain)
colnames(pgcorr)<-"strain"
pgcorr$HiDP<-NA
pgcorr$LoDP<-NA
pgcorr$DSYB<-NA
pgcorr$TpMT2orBoth<-NA

for (i in 1:dim(pgcorr)[1]){
  tmp<-as.character(pgcorr$strain[i])
  a<-unlist(strsplit(tmp,split="-"))[1]
  b<-unlist(strsplit(tmp,split="-"))[2]
  idx<-grep(paste(a,b),pheno_full$Species)
  if(is_empty(idx)){
    idx<-grep(paste(a),pheno_full$Species)
  }
  if(!is_empty(idx)){
    
    if(pheno_full$DMSPType[idx]=="High"){pgcorr$HiDP[i]<-1}
    if(pheno_full$DMSPType[idx]=="High"){pgcorr$LoDP[i]<-0}
    
    if(pheno_full$DMSPType[idx]=="Low"){pgcorr$HiDP[i]<-0}
    if(pheno_full$DMSPType[idx]=="Low"){pgcorr$LoDP[i]<-1}
    
    if(geno_full$DSYB[i]==1 ){pgcorr$DSYB[i]<-1}
    if(geno_full$DSYB[i]==1 ){pgcorr$TpMT2orBoth[i]<-0}
    
    if(geno_full$TpMT2[i]==1 | geno_full$TpMT2_DSYB[i]==1){pgcorr$TpMT2orBoth[i]<-1}
    if(geno_full$TpMT2[i]==1 | geno_full$TpMT2_DSYB[i]==1){pgcorr$DSYB[i]<-0}
  }
}
# Remove species with no phenotype measured
pgcorr<-pgcorr[which(!is.na(pgcorr$HiDP) & !is.na(pgcorr$LoDP)),]
# Remove genotype for None
pgcorr<-pgcorr[which(!is.na(pgcorr$DSYB) & !is.na(pgcorr$TpMT2orBoth)),]
rcorr(as.matrix(pgcorr[,2:5]))

