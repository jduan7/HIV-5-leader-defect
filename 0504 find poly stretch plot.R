setwd("C:/Users/jduan")
library(seqinr)
library(tidyverse)
library(ggplot2)
library(dplyr)
library("readxl")
library("writexl")
library(forcats)


SEQ <- read.alignment("5PSD_master_040422_tilda_uniq.fas",format="fasta") #import the file

HXB2<-c()
for(n in 1:2170){ #HXB2 coordinates 580:2026, ruler 638:2026
  HXB2 <-c(HXB2,substr(SEQ[[3]][1],n,n))}


table<-read_xlsx("C:\\Users\\jduan\\5PSD_master_050922_tilda_uniq_STReverywhere_poly.xlsx") #this one now has selected STR on the first sheet


#to plot stretch everywhere
seq_name<-c()
id<-c()
STR<-c()
pos<-c()
type<-c()

for(i in 1:length(table$stretch)){
  split<-c(unlist(strsplit(table$stretch[i],"")))
  for(j in 638:2157){
    len=length(split)
    if(HXB2[j]=="~"){
      next
    }
    else{
      str<-c(HXB2[j])
      x=1
      str_len=1
      while(str_len<len){
        if(HXB2[j+x]=="~"){
          x=x+1
        }
        else{
          str<-c(str,HXB2[j+x])
          x=x+1
          str_len=str_len+1
        }
      }
      if(identical(str,split)==TRUE){
        spacer=sum(HXB2[1:j+x]=="~",na.rm=TRUE)
        pos<-c(pos,j+x-1-spacer) #record the position of the 3'end of the stretch
        #STR<-c(STR,table$stretch[i])
        type<-c(type,'everywhere')
        seq_name<-c(seq_name,"hxb2")
        id<-c(id,"1")
        if(i<8){
          STR<-c(STR,table$stretch[i])
        }
        else{
          STR<-c(STR,'AAAAATTTT region')
        }
      }
    }
  }
}



#to plot HP

poly_table<-read_xlsx("C:\\Users\\jduan\\5PSD_master_042922_tilda_uniq_checkedHP.xlsx") #this should have the "for genome plot" in the first sheet

for(t in 1:7){
  for(k in 1:length(poly_table$stretch)){
    if(table$stretch[t]==poly_table$stretch[k]){
      STR<-c(STR,table$stretch[t])
      type<-c(type,'5PSD poly')
      seq_name<-c(seq_name,poly_table$Mseqname[k])
      id<-c(id,poly_table$Mid[k])
      pos<-c(pos,poly_table$M_hxb2delend[k]+poly_table$corrected_poly_len[k])
    }
  }
}

for(t in 8:length(table$stretch)){
  for(k in 1:length(poly_table$stretch)){
    if(table$stretch[t]==poly_table$stretch[k]){
      STR<-c(STR,'AAAAATTTT region')
      type<-c(type,'5PSD poly')
      seq_name<-c(seq_name,poly_table$Mseqname[k])
      id<-c(id,poly_table$Mid[k])
      pos<-c(pos,poly_table$M_hxb2delend[k]+poly_table$corrected_poly_len[k])
    }
  }
}


selected_STR_table<-tibble(id,seq_name,pos,STR,type)


#largest to smallest order
selected_STR_table %>%
  mutate(STR=fct_reorder(STR,pos,.fun='length')) %>%
  ggplot(aes(x=pos, y=STR,col=type))+geom_point()+labs(title="Selected HP Stretch Distribution",y="Selected HP Stretch",x="HXB2 Coordinates")




write_xlsx(selected_STR_table,"C:\\Users\\jduan\\5PSD_master_050922_tilda_uniq_selected STR table.xlsx")


#plot everything except AAAAATTTT region
wo_AT<-selected_STR_table[-grep('AAAAATTTT region',selected_STR_table$STR),]
wo_AT %>%
  mutate(STR=fct_reorder(STR,pos,.fun='length')) %>%
  ggplot(aes(x=pos, y=STR,col=type))+geom_point()+labs(title="Selected HP Stretch Distribution",y="Selected HP Stretch",x="HXB2 Coordinates")

#plot just AAAAATTTT region
w_AT<-selected_STR_table[grep('AAAAATTTT region',selected_STR_table$STR),]
w_AT %>%
  mutate(STR=fct_reorder(STR,pos,.fun='length')) %>%
  ggplot(aes(x=pos, y=STR,col=type))+geom_point()+labs(title="Selected HP Stretch Distribution",y="Selected HP Stretch",x="HXB2 Coordinates")

#zoom in on AAAAATTTT region
w_AT %>%
  mutate(STR=fct_reorder(STR,pos,.fun='length')) %>%
  ggplot(aes(x=pos, y=STR,col=type))+geom_point()+labs(title="Selected HP Stretch Distribution",y="Selected HP Stretch",x="HXB2 Coordinates")+scale_x_continuous(breaks=c(seq(755,770)),limits=c(755,770))
