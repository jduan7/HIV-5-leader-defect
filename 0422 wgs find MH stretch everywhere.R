setwd("C:/Users/jduan")
library(seqinr)
library(tidyverse)
library(ggplot2)
library(dplyr)
library("readxl")
library("writexl")



SEQ <- read.alignment("ConsBandHXB2.fas",format="fasta") #import the file
table<-read_xlsx("C:\\Users\\jduan\\wgs_042222_MH stretch counts.xlsx")

consensus<-c()
for(n in 706:9645){ #consensus coordinates 648:9559, wgs ruler 708:10287
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
  freq_everywhere<-c(freq_everywhere,num_str/length(consensus))
}

freq_in_allMH<-table$counts/sum(table$counts)

STR_table<-tibble(STR,table$counts,freq_in_allMH,NUM_STR,freq_everywhere)

write_xlsx(STR_table,"C:\\Users\\jduan\\wgs_042222_STReverywhere.xlsx")









#SEQ <- read.alignment("wgs_clean_unique_0130.fas",format="fasta") #import the file
#table<-read_xlsx("C:\\Users\\jduan\\wgs_042222_MH stretch counts.xlsx")

#seq_name<-c()
#LEN<-c()
#STR<-c()
#NUM_STR<-c()
#freq_everywhere<-c()

#for(r in 2:length(SEQ[[2]])){
#  nseq<-c()
#  for(n in 620:10355){
#    nseq<-c(nseq,substr(SEQ[[3]][r],n,n))}
#  nseq<-nseq[-grep("~",nseq)]
  
#  for(i in 1:length(table$stretch)){
#    split<-c(unlist(strsplit(table$stretch[i],"")))
#    num_str=0
#   for(j in 1:length(nseq)){
#     len=length(split)
#      if(identical(nseq[j:(j+(len-1))],split)==TRUE){
#        num_str=num_str+1
#      }
#    }
#    seq_name<-c(seq_name,SEQ[[2]][r])
#    LEN<-c(LEN,length(nseq))
#    STR<-c(STR,table$stretch[i])
#    NUM_STR<-c(NUM_STR,num_str)
#    freq_everywhere<-c(freq_everywhere,num_str/length(nseq))
#  }
#}


#STR_table<-tibble(seq_name,STR,NUM_STR,LEN,freq_everywhere)
#write_xlsx(STR_table,"C:\\Users\\jduan\\wgs_042222_STReverywhere.xlsx")

#str_avg<-c()
#for(k in 1:length(table$stretch)){
#  str_freq=0
#  for(l in 1:length(STR_table$STR)){
#    if(STR_table$STR[l]==table$stretch[k]){ #1-nt and 2-nt dels are not counted (a hotspot of 2-nt del is at 759nt, a non-coding position)
#      str_freq=str_freq+STR_table$freq_everywhere[l]
#    }
#  }
#  str_avg<-c(str_avg,str_freq/(length(SEQ[[2]])-1))
#}

#STR_avg_table<-tibble(table$stretch,str_avg)


#write_xlsx(STR_avg_table,"C:\\Users\\jduan\\wgs_042222_STR_avg.xlsx")
