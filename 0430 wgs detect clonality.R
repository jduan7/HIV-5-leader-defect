setwd("C:/Users/jduan")
library(seqinr)
library(tidyverse)
library(ggplot2)
library(dplyr)
library("readxl")
library("writexl")

table<-read_xlsx("C:\\Users\\jduan\\wgs_041922_allDel_cleared.xlsx")
Seq<-c(table$SeqName[1])
S_ID<-c(table$ID[1])
del_junction<-c(table$Del_End[1],table$Del_len[1])
Del_junc<-c()

for(v in 2:length(table$SeqName)){
  if(table$SeqName[v]==table$SeqName[v-1]){
    del_junction<-c(del_junction,table$Del_End[v],table$Del_len[v])
    if(v==length(table$SeqName)){
      Del_junc<-c(Del_junc,toString(del_junction))
    }
  }
  else{
    Del_junc<-c(Del_junc,toString(del_junction))
    del_junction<-c(table$Del_End[v],table$Del_len[v])
    Seq<-c(Seq,table$SeqName[v])
    S_ID<-c(S_ID,table$ID[v])
    if(v==length(table$SeqName)){
      Del_junc<-c(Del_junc,toString(del_junction))
    }
  }
}
S_table<-tibble(S_ID,Seq,Del_junc)
write_xlsx(S_table,"C:\\Users\\jduan\\wgs_0430_S_table_uniq.xlsx")
