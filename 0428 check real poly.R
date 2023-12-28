setwd("C:/Users/jduan")
library(seqinr)
library(tidyverse)
library(ggplot2)
library(dplyr)
library("readxl")
library("writexl")


SEQ <- read.alignment("5PSD_master_041222_tilda_uniq_consensus.fas",format="fasta") #import the file

table<-read_xlsx("C:\\Users\\jduan\\5PSD_master_041922_tilda_uniq_allHP.xlsx") #raw data

###################to get just the ones under 3nts at 3'end
remove<-c()
for (u in 1:length(table$ID)){
  if(table$poly_3end[u]>=3){
    remove<-c(remove,u)
  }
}
table<-table[-c(remove),]

consensus<-c()
for(n in 1:2157){
  consensus <-c(consensus,substr(SEQ[[3]][1],n,n))}#split hxb2 into individual nucleotides




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
poly_MH_5end_nt<-c()
poly_MH_3end_nt<-c()
poly_MH_type<-c()
poly_MH<-c()
poly_len<-c()

for(y in 1:length(table[[1]])){#for all the entries in the poly-table
  m=table[[1]][y]
  S=table[[4]][y] #ruler del start
  E=table[[8]][y] #ruler del end
  nseq<-c() #for extracting the sequence string into individual nucleotides
  spacer_S=0
  spacer_E=0
  for(n in 1:2157){
    nseq<-c(nseq,substr(SEQ[[3]][m],n,n))} #m is the respective seq in the fas file
  
  #############only need to check for reverse direction MH (b/c we only care about those that can be concatenated to the 3'end)
  
  for (i in 0:30){ #sliding window of 1nt to 15nts starting from each end of the del junction, going 3' to 5'
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
        if(mlen>=2){#if MH is present
          print("There is reverse microhomology.")
          Mdelstart<-c(Mdelstart,S) #gives ruler position of the deletion junction that has microhomology
          Mlen<-c(Mlen,mlen)
          Mdelend<-c(Mdelend,E)
          Mid<-c(Mid,m)
          Mseq<-c(Mseq,toString(nseq[(S-i+spacer_S):(S-1)])) #flip it
          Mseqname<-c(Mseqname,table$SeqName[y])
          Mdir<-c(Mdir,"reverse")
          
          M_hxb2delstart<-c(M_hxb2delstart,table$HXB2_Delstart[y])
          M_hxb2delend<-c(M_hxb2delend,table$HXB2_Delend[y])
          M_dlen<-c(M_dlen,table$Del_len[y])
          
          poly_MH_5end_nt<-c(poly_MH_5end_nt,table$poly_5end_nt[y])
          poly_MH_3end_nt<-c(poly_MH_3end_nt,table$poly_3end_nt[y])
          
          
          if(table$poly_5end[y]<=mlen){#if poly is within MH, use poly length
            if((unlist(strsplit(table$poly_5end_nt[y],", "))[1])==(unlist(strsplit(table$poly_3end_nt[y],", "))[1])){#can only concatenate if it's MH and it's same poly-n
              if(table$poly_5end[y]+table$poly_3end[y]>=3){
                poly_len<-c(poly_len,table$poly_5end[y]+table$poly_3end[y])
                poly_MH_type<-c(poly_MH_type,"poly within MH")
                poly_MH<-c(poly_MH,toString(table$poly_5end_nt[y],table$poly_3end_nt[y]))
              }
              else{
                poly_len<-c(poly_len,NA)
                poly_MH_type<-c(poly_MH_type,"no poly")
                poly_MH<-c(poly_MH,NA)
              }
            }
            else{#cannot concatenate
              if(table$poly_5end[y]>=3){
                poly_len<-c(poly_len,table$poly_5end[y])
                poly_MH_type<-c(poly_MH_type,"poly within MH")
                poly_MH<-c(poly_MH,table$poly_5end_nt[y])
              }
              else{
                poly_len<-c(poly_len,NA)
                poly_MH_type<-c(poly_MH_type,"no poly")
                poly_MH<-c(poly_MH,NA)
              }
            }
          }
          else{#if poly extends beyond MH, use mlen
            if((unlist(strsplit(table$poly_5end_nt[y],", "))[1])==(unlist(strsplit(table$poly_3end_nt[y],", "))[1])){#can only concatenate if it's MH and it's same poly-n
              if(mlen+table$poly_3end[y]>=3){
                poly_len<-c(poly_len,mlen+table$poly_3end[y])
                poly_MH_type<-c(poly_MH_type,"poly-MH")
                poly_MH<-c(poly_MH,toString(toString(nseq[(S-i+spacer_S):(S-1)]),table$poly_3end_nt[y]))
              }
              else{
                poly_len<-c(poly_len,NA)
                poly_MH_type<-c(poly_MH_type,"no poly")
                poly_MH<-c(poly_MH,NA)
              }
            }
            else{#cannot concatenate
              if(mlen>=3){
                poly_len<-c(poly_len,mlen)
                poly_MH_type<-c(poly_MH_type,"poly-MH")
                poly_MH<-c(poly_MH,toString(nseq[(S-i+spacer_S):(S-1)]))
              }
              else{
                poly_len<-c(poly_len,NA)
                poly_MH_type<-c(poly_MH_type,"no poly")
                poly_MH<-c(poly_MH,NA)
              }
            }
          }
          
          break
        }
        else{#MH is not present
          print("There is no reverse microhomology. No true poly.")
          break
        }
      }
    }
  }
}


######################################

table<-read_xlsx("C:\\Users\\jduan\\5PSD_master_041922_tilda_uniq_allHP.xlsx") #do it all over again

###################to get just the ones equal to or more than 3nts at 3'end
remove<-c()
for (u in 1:length(table$ID)){
  if(table$poly_3end[u]<3){
    remove<-c(remove,u)
  }
}
table<-table[-c(remove),]


for(y in 1:length(table[[1]])){#for all the entries in the poly-table
  m=table[[1]][y]
  S=table[[4]][y] #ruler del start
  E=table[[8]][y] #ruler del end
  nseq<-c() #for extracting the sequence string into individual nucleotides
  spacer_S=0
  spacer_E=0
  for(n in 1:2157){
    nseq<-c(nseq,substr(SEQ[[3]][m],n,n))} #m is the respective seq in the fas file
  
  #############only need to check for reverse direction MH (b/c we only care about those that can be concatenated to the 3'end)
  
  for (i in 0:30){ #sliding window of 1nt to 15nts starting from each end of the del junction, going 3' to 5'
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
        if(mlen>=2){#if MH is present
          print("There is reverse microhomology.")
          Mdelstart<-c(Mdelstart,S) #gives ruler position of the deletion junction that has microhomology
          Mlen<-c(Mlen,mlen)
          Mdelend<-c(Mdelend,E)
          Mid<-c(Mid,m)
          Mseq<-c(Mseq,toString(nseq[(S-i+spacer_S):(S-1)])) #flip it
          Mseqname<-c(Mseqname,table$SeqName[y])
          Mdir<-c(Mdir,"reverse")
          
          M_hxb2delstart<-c(M_hxb2delstart,table$HXB2_Delstart[y])
          M_hxb2delend<-c(M_hxb2delend,table$HXB2_Delend[y])
          M_dlen<-c(M_dlen,table$Del_len[y])
          
          poly_MH_5end_nt<-c(poly_MH_5end_nt,table$poly_5end_nt[y])
          poly_MH_3end_nt<-c(poly_MH_3end_nt,table$poly_3end_nt[y])
          
          
          if(table$poly_5end[y]<=mlen){#if poly is within MH, use poly length
            if((unlist(strsplit(table$poly_5end_nt[y],", "))[1])==(unlist(strsplit(table$poly_3end_nt[y],", "))[1])){#can only concatenate if it's MH and it's same poly-n
              poly_len<-c(poly_len,table$poly_5end[y]+table$poly_3end[y])
              poly_MH_type<-c(poly_MH_type,"poly-MH")
              poly_MH<-c(poly_MH,toString(table$poly_5end_nt[y],table$poly_3end_nt[y]))
            }
            else{#cannot concatenate
              poly_MH<-c(poly_MH,toString(toString(nseq[(S-i+spacer_S):(S-1)]),table$poly_3end_nt[y]))
              poly_len<-c(poly_len,table$poly_3end[y])
              poly_MH_type<-c(poly_MH_type,"5end MH proximal to 3end poly")
            }
          }
          else{#if poly extends beyond MH, use mlen
            if((unlist(strsplit(table$poly_5end_nt[y],", "))[1])==(unlist(strsplit(table$poly_3end_nt[y],", "))[1])){#can only concatenate if it's MH and it's same poly-n
              poly_len<-c(poly_len,mlen+table$poly_3end[y])
              poly_MH_type<-c(poly_MH_type,"poly-MH")
              poly_MH<-c(poly_MH,toString(toString(nseq[(S-i+spacer_S):(S-1)]),table$poly_3end_nt[y]))
            }
            else{#cannot concatenate
              poly_MH<-c(poly_MH,toString(toString(nseq[(S-i+spacer_S):(S-1)]),table$poly_3end_nt[y]))
              poly_len<-c(poly_len,table$poly_3end[y])
              poly_MH_type<-c(poly_MH_type,"5end MH proximal to 3end poly")
            }
          }
          
          break
        }
        else{#5end MH is not present
          print("There is no reverse microhomology. No true 5end poly. Only true 3end poly")
          
          m=table[[1]][y]
          S=table[[4]][y] #ruler del start
          E=table[[8]][y] #ruler del end
          nseq<-c() #for extracting the sequence string into individual nucleotides
          spacer_S=0
          spacer_E=0
          for(n in 1:2157){
            nseq<-c(nseq,substr(SEQ[[3]][m],n,n))}
          
          for(j in 0:30){
            if (consensus[S+j-spacer_S]==nseq[E+1+j-spacer_E]&&nseq[E+1+j-spacer_E]!="~"){
              next
            }
            else{
              if((nseq[E+1+j-spacer_E]=="~")&&(consensus[E+1+j-spacer_E]==nseq[E+1+j-spacer_E])){
                spacer_S=spacer_S+1
                next
              }
              else if((consensus[S+j-spacer_S]=="~")&&(consensus[S+j-spacer_S]==nseq[S+j-spacer_S])){
                spacer_E=spacer_E+1
                next
              }
              
              else{
                mlen=j-spacer_E-spacer_S
                if(mlen>=2){
                  print("There is forward microhomology.")
                  Mdelstart<-c(Mdelstart,S) #gives ruler position of the deletion junction that has microhomology
                  Mlen<-c(Mlen,mlen)
                  Mdelend<-c(Mdelend,E)
                  Mid<-c(Mid,m)
                  Mseq<-c(Mseq,toString(nseq[(E+1):(E+j-spacer_E)]))
                  Mseqname<-c(Mseqname,table$SeqName[y])
                  Mdir<-c(Mdir,"forward")
                  
                  M_hxb2delstart<-c(M_hxb2delstart,table$HXB2_Delstart[y])
                  M_hxb2delend<-c(M_hxb2delend,table$HXB2_Delend[y])
                  M_dlen<-c(M_dlen,table$Del_len[y])
                  
                  poly_MH_5end_nt<-c(poly_MH_5end_nt,table$poly_5end_nt[y])
                  poly_MH_3end_nt<-c(poly_MH_3end_nt,table$poly_3end_nt[y])
                  
                  if(mlen<=table$poly_3end[y]){
                    poly_MH<-c(poly_MH,toString(nseq[(E+1):(E+j-spacer_E)]))
                    poly_MH_type<-c(poly_MH_type,"MH within poly on 3end")
                    poly_len<-c(poly_len,table$poly_3end[y])
                  }
                  else{
                    poly_MH<-c(poly_MH,table$poly_3end_nt[y])
                    poly_MH_type<-c(poly_MH_type,"poly within MH on 3end")
                    poly_len<-c(poly_len,table$poly_3end[y])
                  }
                  break
                }
                else{
                  print("There is no forward microhomology.")
                  
                  Mdelstart<-c(Mdelstart,S) #gives ruler position of the deletion junction that has microhomology
                  Mlen<-c(Mlen,NA)
                  Mdelend<-c(Mdelend,E)
                  Mid<-c(Mid,m)
                  Mseq<-c(Mseq,NA)
                  Mseqname<-c(Mseqname,table$SeqName[y])
                  Mdir<-c(Mdir,"forward")
                  
                  M_hxb2delstart<-c(M_hxb2delstart,table$HXB2_Delstart[y])
                  M_hxb2delend<-c(M_hxb2delend,table$HXB2_Delend[y])
                  M_dlen<-c(M_dlen,table$Del_len[y])
                  
                  poly_MH_5end_nt<-c(poly_MH_5end_nt,table$poly_5end_nt[y])
                  poly_MH_3end_nt<-c(poly_MH_3end_nt,table$poly_3end_nt[y])
                  
                  poly_MH_type<-c(poly_MH_type,"3end poly")
                  poly_MH<-c(poly_MH,NA)
                  poly_len<-c(poly_len,table$poly_3end[y])
                  break
                }
              }
            }
          }
          break
          
          
          
          
          
          
          #Mdelstart<-c(Mdelstart,S) #gives ruler position of the deletion junction that has microhomology
          #Mlen<-c(Mlen,mlen)
          #Mdelend<-c(Mdelend,E)
          #Mid<-c(Mid,m)
          #Mseq<-c(Mseq,NA) #flip it
          #Mseqname<-c(Mseqname,table$SeqName[y])
          #Mdir<-c(Mdir,"reverse")
          
          #M_hxb2delstart<-c(M_hxb2delstart,table$HXB2_Delstart[y])
          #M_hxb2delend<-c(M_hxb2delend,table$HXB2_Delend[y])
          #M_dlen<-c(M_dlen,table$Del_len[y])
          
          #poly_MH_5end_nt<-c(poly_MH_5end_nt,table$poly_5end_nt[y])
          #poly_MH_3end_nt<-c(poly_MH_3end_nt,table$poly_3end_nt[y])
          
          #poly_len<-c(poly_len,table$poly_3end[y])
          #poly_MH<-c(poly_MH,NA) #waiting to be further analyzed
          #poly_MH_type<-c(poly_MH_type,NA) #waiting to be further analyzed
          #break
        }
      }
    }
  }
}


###########################################








###########################################

real_poly<-tibble(Mid,Mseqname,Mdelstart,Mdelend,Mlen,Mdir,Mseq,poly_MH_5end_nt,poly_MH_3end_nt,poly_len,poly_MH_type,poly_MH,M_hxb2delstart,M_hxb2delend,M_dlen)

write_xlsx(real_poly,"C:\\Users\\jduan\\5PSD_master_042922_tilda_uniq_checkedHP.xlsx")
