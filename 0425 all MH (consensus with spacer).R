setwd("C:/Users/jduan")
library(seqinr)
library(tidyverse)
library(ggplot2)
library(dplyr)
library("readxl")
library("writexl")



########################GET DELETION RULER POSITIONS

table<-read_xlsx("C:\\Users\\jduan\\5PSD_master_041922_tilda_uniq_allDel_cleared.xlsx")

##############################Get all MH


SEQ <- read.alignment("5PSD_master_041222_tilda_uniq_consensus.fas",format="fasta") #import the file

Mdelstart<-c() #list of deletion starts that have microhomology
Mdelend<-c() #list of deletion ends that have microhomology
Mlen<-c() #list of microhomology length
Mid<-c() #list of Gid that has microhomology
Mseq<-c() #list of sequences with microhomology
Mseqname<-c() #list of names for sequences with microhomology at del junction
Mdir<-c()#list of directions for sequences with microhomology
M_hxb2delstart<-c()
M_hxb2delend<-c()
M_dlen<-c()
M_ConsDelstart<-c()
M_ConsDelend<-c()

consensus<-c()
for(n in 1:2157){
  consensus <-c(consensus,substr(SEQ[[3]][1],n,n))}#split hxb2 into individual nucleotides

for(y in 1:length(table[[1]])){#for all the entries in the poly-table
  m=table[[1]][y]
  S=table[[3]][y] #ruler del start
  E=table[[4]][y] #ruler del end
  nseq<-c() #for extracting the sequence string into individual nucleotides
  spacer_S=0
  spacer_E=0
  for(n in 1:2157){
    nseq<-c(nseq,substr(SEQ[[3]][m],n,n))} #m is the respective seq in the fas file
  
  
  for (i in 0:30){ #sliding window of 1nt to 15nts starting from each end of the del junction towards 5'
    if (consensus[S+i-spacer_S]==nseq[E+1+i-spacer_E]&&nseq[E+1+i-spacer_E]!="~"){
      next
    }
    else{
      if((nseq[E+1+i-spacer_E]=="~")&&(consensus[E+1+i-spacer_E]==nseq[E+1+i-spacer_E])){
        spacer_S=spacer_S+1
        next
      }
      else if((consensus[S+i-spacer_S]=="~")&&(consensus[S+i-spacer_S]==nseq[S+i-spacer_S])){
        spacer_E=spacer_E+1
        next
      }
      
      else{
        mlen=i-spacer_E-spacer_S
        if(mlen>=2){
          print("There is forward microhomology.")
          Mdelstart<-c(Mdelstart,S) #gives ruler position of the deletion junction that has microhomology
          Mlen<-c(Mlen,mlen)
          Mdelend<-c(Mdelend,E)
          Mid<-c(Mid,m)
          Mseq<-c(Mseq,toString(nseq[(E+1):(E+i-spacer_E)]))
          Mseqname<-c(Mseqname,table$SeqName[y])
          Mdir<-c(Mdir,"forward")
          
          M_ConsDelstart<-c(M_ConsDelstart,S-sum(consensus[1:S]=="~",na.rm=TRUE))
          M_ConsDelend<-c(M_ConsDelend,E-sum(consensus[1:E]=="~",na.rm=TRUE))
          
          M_hxb2delstart<-c(M_hxb2delstart,table$HXB2_Delstart[y])
          M_hxb2delend<-c(M_hxb2delend,table$HXB2_Delend[y])
          M_dlen<-c(M_dlen,table$Del_len[y])
          break
        }
        else{
          print("There is no forward microhomology.")
          break
        }
      }
    }
  }
  
  
  
  
  for (i in 0:30){ #sliding window of 1nt to 15nts starting from each end of the del junction towards 3'
    if (consensus[E-i+spacer_E]==nseq[S-1-i+spacer_S]&&nseq[S-1-i+spacer_S]!="~"){
      next
    }
    else{
      if((nseq[S-i-1+spacer_S]=="~")&&(consensus[S-i-1+spacer_S]==nseq[S-i-1+spacer_S])){
        spacer_E=spacer_E+1
        next
      }
      else if((consensus[E-i+spacer_E]=="~")&&(consensus[E-i+spacer_E]==nseq[E-i+spacer_E])){
        spacer_S=spacer_S+1
        next
      }
      
      else{
        mlen=i-spacer_E-spacer_S
        if(mlen>=2){
          print("There is reverse microhomology.")
          Mdelstart<-c(Mdelstart,S) #gives ruler position of the deletion junction that has microhomology
          Mlen<-c(Mlen,mlen)
          Mdelend<-c(Mdelend,E)
          Mid<-c(Mid,m)
          Mseq<-c(Mseq,toString(nseq[(S-i+spacer_S):(S-1)])) #flip it
          Mseqname<-c(Mseqname,table$SeqName[y])
          Mdir<-c(Mdir,"reverse")
          
          M_ConsDelstart<-c(M_ConsDelstart,S-sum(consensus[1:S]=="~",na.rm=TRUE))
          M_ConsDelend<-c(M_ConsDelend,E-sum(consensus[1:E]=="~",na.rm=TRUE))
          
          M_hxb2delstart<-c(M_hxb2delstart,table$HXB2_Delstart[y])
          M_hxb2delend<-c(M_hxb2delend,table$HXB2_Delend[y])
          M_dlen<-c(M_dlen,table$Del_len[y])
          break
        }
        else{
          print("There is no reverse microhomology.")
          break
        }
      }
    }
  }
}

all_MH<-tibble(Mid,Mseqname,Mdelstart,Mdelend,Mlen,Mdir,Mseq,M_ConsDelstart,M_ConsDelend,M_dlen,M_hxb2delstart,M_hxb2delend)      

write_xlsx(all_MH,"C:\\Users\\jduan\\5PSD_master_042522_tilda_uniq_consensus_WithSpacer_allMH.xlsx")
