setwd("C:/Users/jduan")
library(seqinr)
library(tidyverse)
library(ggplot2)
library(dplyr)
library("readxl")
library("writexl")
library(forcats)


SEQ <- read.alignment("wgs_clean_unique_0130.fas",format="fasta") #import the file

consensus<-c()
for(n in 1:10322){ #hxb2 coordinates 1:9672
  consensus <-c(consensus,substr(SEQ[[3]][1],n,n))}


table<-read_xlsx("C:\\Users\\jduan\\wgs_042622_STReverywhere.xlsx")


#to plot stretch everywhere
seq_name<-c()
id<-c()
STR<-c()
pos<-c()
type<-c()

for(i in 1:length(table$stretch)){
  split<-c(unlist(strsplit(table$stretch[i],"")))
  for(j in 708:10287){
    len=length(split)
    if(consensus[j]=="~"){
      next
    }
    else{
      str<-c(consensus[j])
      x=1
      str_len=1
      while(str_len<len){
        if(consensus[j+x]=="~"){
          x=x+1
        }
        else{
          str<-c(str,consensus[j+x])
          x=x+1
          str_len=str_len+1
        }
      }
      if(identical(str,split)==TRUE){
        spacer=sum(consensus[1:j+x]=="~",na.rm=TRUE)
        pos<-c(pos,j+x-1-spacer) #record the position of the 3'end of the stretch
        STR<-c(STR,table$stretch[i])
        type<-c(type,'everywhere')
        seq_name<-c(seq_name,"consensus")
        id<-c(id,"1")
      }
    }
  }
}


#to plot MH

MH_table<-read_xlsx("C:\\Users\\jduan\\wgs_041922_allMH_updated.xlsx") #this should have the "for genome plot" in the first sheet

for(t in 1:length(table$stretch)){
  for(k in 1:length(MH_table$Mseq)){
    if(table$stretch[t]==MH_table$Mseq[k]){
      STR<-c(STR,table$stretch[t])
      type<-c(type,'5PSD MH')
      seq_name<-c(seq_name,MH_table$Mseqname[k])
      id<-c(id,MH_table$Mid[k])
      if(MH_table$Mdir[k]=="reverse"){#if it's aligned in the reverse direction
        pos<-c(pos,MH_table$M_hxb2delend[k]) #record the position of the 3'end of the stretch
      }
      else{#if it's aligned in the forward direction
        pos<-c(pos,MH_table$M_hxb2delend[k]+MH_table$Mlen[k])
      }
    }
  }
}


selected_STR_table<-tibble(id,seq_name,pos,STR,type)


write_xlsx(selected_STR_table,"C:\\Users\\jduan\\wgs_042922_tilda_uniq_selected STR table.xlsx")

selected_STR_table %>%
  mutate(STR=fct_reorder(STR,pos,.fun='length')) %>%
  ggplot(aes(x=pos, y=STR,col=type))+geom_point()+labs(title="Selected MH Stretch Distribution",y="Selected MH Stretch",x="HXB2 Coordinates")


ggplot(selected_STR_table,aes(x=pos, y=STR,col=type))+geom_point()+labs(title="Selected MH Stretch Distribution",y="Selected MH Stretch",x="Consensus B Coordinates")
