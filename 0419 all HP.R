setwd("C:/Users/jduan")
library(seqinr)
library(tidyverse)
library(ggplot2)
library(dplyr)
library("readxl")
library("writexl")


SEQ <- read.alignment("5PSD_master_040422_tilda_uniq.fas",format="fasta") #import the file

table<-read_xlsx("C:\\Users\\jduan\\5PSD_master_041922_tilda_uniq_allDel_cleared.xlsx")

poly_5end<-c()
poly_5end_nt<-c()
poly_3end<-c()
poly_3end_nt<-c()

for(r in 1:length(table$ID)){#for each deletion
  nseq<-c()
  for(n in 1:2157){
    nseq<-c(nseq,substr(SEQ[[3]][table$ID[r]],n,n))} #get the sequence of the deletion
  
  Cstart<-table$Del_Start[r] #deletion of this round
  for (i in 1:15){ #sliding window of 1nts to 10nts starting from each end of the del junction towards 5'
    if (nseq[Cstart-i]==nseq[Cstart-i-1]){
      if(nseq[Cstart-i]!="~"){
        next
      }
      else{
        break
      }
    }
    else{
      if(Cstart>=760&&Cstart<=768){#if it's in that AAAAATTTT region
        if(nseq[Cstart-i-2]==nseq[Cstart-i-3]){
          if(nseq[Cstart-i-2]!="~"){
            next
          }
          else{
            break
          }
        }
        else{
          break
        }
      }
      else{
        break
      }
    }
  }
  if(i>1){
    while(nseq[Cstart-i]!=nseq[Cstart-i+1]){
      i=i+1
    }
  }
  poly_5end<-c(poly_5end,i)
  poly_5end_nt<-c(poly_5end_nt,toString(nseq[(Cstart-1):(Cstart-i)]))
  
  Cend<-table$Del_End[r] #deletion of this round
  for (j in 1:15){ #sliding window of 1nts to 10nts starting from each end of the del junction towards 3'
    if (identical(nseq[Cend+j],nseq[Cend+j+1])==TRUE){
      if(nseq[Cend+j]!="~"){
        next
      }
      else{
        break
      }
    }
    else{
      if(Cend>=760&&Cend<=768){#if it's in that AAAAATTTT region
        if(nseq[Cend+j+2]==nseq[Cend+j+3]){
          if(nseq[Cend+j+2]!="~"){
            next
          }
          else{
            break
          }
        }
        else{
          break
        }
      }
      else{
        break
      }
    }
  }
  if(j>1){
    while(nseq[Cend+j]!=nseq[Cend+j-1]){
      j=j-1
    }
  }
  poly_3end<-c(poly_3end,j)
  poly_3end_nt<-c(poly_3end_nt,toString(nseq[(Cend+1):(Cend+j)]))
}

HP_table<-tibble(table$ID,table$SeqName,table$Del_len,table$Del_Start,table$HXB2_Delstart,poly_5end,poly_5end_nt,table$Del_End,table$HXB2_Delend,poly_3end,poly_3end_nt)  

write_xlsx(HP_table,"C:\\Users\\jduan\\5PSD_master_041922_tilda_uniq_allHP.xlsx")
