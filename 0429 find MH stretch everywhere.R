setwd("C:/Users/jduan")
library(seqinr)
library(tidyverse)
library(ggplot2)
library(dplyr)
library("readxl")
library("writexl")


SEQ <- read.alignment("ConsBandHXB2.fas",format="fasta") #import the file

table<-read_xlsx("C:\\Users\\jduan\\5PSD_master_041922_MH stretch counts.xlsx")

consensus<-c()
for(n in 638:2085){ #consensus coordinates 580:2026, ruler 638:2026
  consensus <-c(consensus,substr(SEQ[[3]][2],n,n))}

consensus<-consensus[-grep("-",consensus)]


STR<-c()
NUM_STR<-c()
freq_everywhere<-c()

for(i in 1:length(table$stretch)){
  split<-c(unlist(strsplit(table$stretch[i],"")))
  num_str=0
  for(j in 1:length(consensus)){
    len=length(split)
    if(identical(consensus[j:(j+(len-1))],split)==TRUE){
      num_str=num_str+1
    }
  }
  STR<-c(STR,table$stretch[i])
  NUM_STR<-c(NUM_STR,num_str)
  #freq_everywhere<-c(freq_everywhere,num_str/(length(consensus)))
}

freq_in_allMH<-table$counts/sum(table$counts)

STR_table<-tibble(STR,table$counts,freq_in_allMH,NUM_STR,freq_everywhere)
write_xlsx(STR_table,"C:\\Users\\jduan\\5PSD_master_042922_tilda_uniq_STReverywhere.xlsx")

  

  
  