setwd("C:/Users/jduan")
library(seqinr)
library(tidyverse)
library(ggplot2)
library(dplyr)
library("readxl")
library("writexl")


SEQ <- read.alignment("5PSD_master_040422_tilda_uniq.fas",format="fasta") #import the file

########################GET DELETION RULER POSITIONS


seqname<-c()
Dstart<-c()
Dend<-c()
Dlen<-c()
hxb2_Dstart<-c()
hxb2_Dend<-c()

hxb2<-c()
for(n in 1:2157){
  hxb2 <-c(hxb2,substr(SEQ[[3]][1],n,n))}#split hxb2 into individual nucleotides

for(r in 2:length(SEQ[[2]])){
  nseq<-c()
  delstart<-c()
  delend<-c()
  dlen<-c()
  hxb2_delstart<-c()
  hxb2_delend<-c()
  for(n in 1:2157){
    nseq<-c(nseq,substr(SEQ[[3]][r],n,n))}
  for(i in 2:2157){
    if(identical(nseq[i],"~")==TRUE){
      x=0
      while(identical(nseq[i-1],"~")==FALSE){
        if(identical(nseq[i+x],"~")==TRUE&&(identical(nseq[i+x+1],"~")==FALSE)){
          delstart<-c(delstart,i) #ruler of del start
          delend<-c(delend,i+x) #ruler of del end
          
          c=0
          for(n in 1:i){ #go back to ruler position of del start
            if (hxb2[n]=="~"){
              c=c+1 #count number of spacers
            }
          }
          if(identical(hxb2[i],nseq[i])==TRUE){ #if it's a deletion start directly after a spacer
            hxb2_delstart<-c(hxb2_delstart,n-c+1)
          }
          else{
            hxb2_delstart<-c(hxb2_delstart,n-c) #delete the number of spacers to get hxb2 position of del start
          }
          
          p=0
          for(k in 1:(i+x)){
            if (hxb2[k]=="~"){
              p=p+1
            }
          }
          if(identical(hxb2[i],nseq[i])==TRUE){
            hxb2_delend<-c(hxb2_delend,k-p+1)
          }
          else{
            hxb2_delend<-c(hxb2_delend,k-p) #delete the number of spacers to get hxb2 position of del end
          }
          
          
          if(identical(hxb2[i],nseq[i])==TRUE){ #if it's a deletion start directly after a spacer
            dlen<-c(dlen,(k-p)-(n-c+1)+1) 
          }
          else{
            dlen<-c(dlen,(k-p)-(n-c)+1) #get the length of deletion using hxb2 coordinates
          }
          break
        }
          
        else if(identical(nseq[i+x],"~")==TRUE&&(identical(hxb2[i+x],nseq[i+x])==TRUE)&&(identical(nseq[i+x+1],"~")==FALSE)){
          break
        }
        else{
          x=x+1
        }
      }
    }
  }
  Dstart<-c(Dstart,delstart)
  Dend<-c(Dend,delend)
  Dlen<-c(Dlen,dlen)
  hxb2_Dstart<-c(hxb2_Dstart,hxb2_delstart)
  hxb2_Dend<-c(hxb2_Dend,hxb2_delend)
  seqname<-c(seqname,rep(SEQ[[2]][r],length(delstart)))
}


pre_table<-tibble(seqname,Dstart,Dend,Dlen,hxb2_Dstart,hxb2_Dend)

Del_Start<-c()
Del_End<-c()
Del_len<-c()
HXB2_Delstart<-c()
HXB2_Delend<-c()
SeqName<-c()
for(i in 1:length(Dend)){
  if(identical(pre_table$Dend[i],2157)==FALSE){
    Del_Start<-c(Del_Start,pre_table$Dstart[i])
    Del_End<-c(Del_End,pre_table$Dend[i])
    Del_len<-c(Del_len,pre_table$Dlen[i])
    HXB2_Delstart<-c(HXB2_Delstart,pre_table$hxb2_Dstart[i])
    HXB2_Delend<-c(HXB2_Delend,pre_table$hxb2_Dend[i])
    SeqName<-c(SeqName,pre_table$seqname[i])
  }
}


ID<-c()
for(s in 1:length(SeqName)){
  for(t in 2:length(SEQ[[2]]))
    if(SeqName[s]==SEQ[[2]][t]){
      ID<-c(ID,t)
    }
}

table<-tibble(ID,SeqName,Del_Start,Del_End,Del_len,HXB2_Delstart,HXB2_Delend)

remove<-c()
for(w in 1:length(table$ID)){
  if(table$Del_Start[w]>=1155&&table$Del_End[w]<=1237){#remove complicated region
    remove<-c(remove,w)
  }
  if(table$Del_Start[w]>=1955&&table$Del_End[w]<=1986){#remove complicated region
    remove<-c(remove,w)
  }
  if(table$Del_len[w]==0){
    remove<-c(remove,w)
  }
}
for(r in 1:length(table$ID)){ #for the number of rows in the table (aka number of microhomology seqs)
  input=table$Del_len[r]
  if(table$Del_len[r]<10){
    fsh=(input%%3)==0
    if(fsh==TRUE){
      remove<-c(remove,r)
    }
  }
}

for(r in 1:length(table$ID)){
  if(table$Del_len[r]<3){ #1-nt and 2-nt dels are not counted (a hotspot of 2-nt del is at 759nt, a non-coding position)
    remove<-c(remove,r)
  }
}

table<-table[-c(remove),]


write_xlsx(table,"C:\\Users\\jduan\\5PSD_master_041922_tilda_uniq_allDel_cleared.xlsx")
