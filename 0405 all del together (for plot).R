######all del plotted together
setwd("C:/Users/jduan")
library(seqinr)
library(tidyverse)
library(ggplot2)
library(dplyr)
library("writexl")


#1 get all the sequences
SEQ <- read.alignment("5PSD_master_040522_hxb2.fas",format="fasta") #import the file
names<-c()
AllSeq<-c()
for(i in 1:SEQ[[1]]){ #for all the sequences
  names<-c(names,SEQ[[2]][i]) #put names into a list
  AllSeq<-c(AllSeq,SEQ[[3]][i]) #put each sequence as an element in a list
}

archGstart<-c() #build an archive
archGend<-c()
archNames<-c() 
archLength<-c()
namesA<-c()
Gid<-c()
source<-c()
sourceA<-c()
sourceB<-c()
##################################################################################
all_Gpos<-c() #all gap positions for all sequences
all_Gid<-c() #each sequence gets an id, corresponding to the y-value on graph
all_Gset<-c() #Gpos and Gid matched
all_ThreePrime<-c()
all_listA<-c() #list of 3' and associated sequence name, later sorting it and adding id in descending order
all_listB<-c() #list of positions with sequence names, later aligning to sequence names in sorted listA
all_namesA<-c()
all_namesB<-c()#list of sequence names

#2 start a for loop
for(i in 2:SEQ[[1]]){
  seq<-c(AllSeq[i]) #each round of loop, the individual seq is a new one
  nseq<-c() #extracting the sequence string into individual nucleotides
  for(n in 1:2157){ #here, 2000 can be changed to other numbers according to how much we like
    nseq<-c(nseq,substr(AllSeq[i],n,n))
  }
  
  
  #3 get the start and end of each deletion for the sequence in this round of loop
  
  start<-c() #list of del start positions
  end<-c() #list of del end positions
  p=0
  for(m in 1:length(nseq)){
    p=p+1 #add up 1 on position number for each round
    if(nseq[m]=="-"){
      if (m==1){
        if(nseq[m]=="-"){
          if(nseq[m+1]!="-"){
            start<-c(start,p)
            end<-c(end,p)
          }
          else{
            start<-c(start,p)
          }
        }
      }
      else if (m==length(nseq)){
        if(nseq[m]=="-"){
          end<-c(end,p) 
        }
      }
      else if (m!=1 & m!=length(nseq)){
        if(nseq[m-1]=="t"|nseq[m-1]=="c"|nseq[m-1]=="g"|nseq[m-1]=="a"|nseq[m-1]=="~"){
          if(nseq[m+1]=="t"|nseq[m+1]=="c"|nseq[m+1]=="g"|nseq[m+1]=="a"|nseq[m+1]=="~"){
            start<-c(start,p)
            end<-c(end,p)
          }
          else {
            start<-c(start,p)
          }
        } 
        else if(nseq[m+1]=="t"|nseq[m+1]=="c"|nseq[m+1]=="g"|nseq[m+1]=="a"|nseq[m+1]=="~") {
          end<-c(end,p)
        }
      }
    }
    else if (nseq[m]=="~"){ 
      p=p-1
    }
  }
  
  
  #4 gives the start and end of the DEL regions
  Gstart<-c(start[1])
  Gend<-c()
  for (w in 1:length(start)){
    if(w==length(start)){
      Gend<-c(Gend,end[w])
    }
    else{
      if(start[w+1]-end[w]!=1){
        Gstart<-c(Gstart,start[w+1])
        Gend<-c(Gend,end[w])
      }
    }
  }
  
  
  #6 adding source
  if(i>1 & i<59){
    sourceA<-c(sourceA,'FS')
  }
  else if(i>58 & i<67){
    sourceA<-c(sourceA,'unpFS')
  }
  else if(i==67){
    sourceA<-c(sourceA,'SIM2')
  }
  else if(i>67 & i<75){
    sourceA<-c(sourceA,'HIE')
  }
  else if(i>74 & i<78){
    sourceA<-c(sourceA,'HOR')
  }
  else if(i>77 & i<104){
    sourceA<-c(sourceA,'ANT')
  }
  else if(i>103 & i<113){
    sourceA<-c(sourceA,'BRU')
  }
  else if(i>112 & i<120){
    sourceA<-c(sourceA,'HOX')
  }
  else if(i>119 & i<213){
    sourceA<-c(sourceA,'GAE')
  }
  else if(i>212 & i<216){
    sourceA<-c(sourceA,'MEN')
  }
  else if(i>215 & i<471){
    sourceA<-c(sourceA,'LIC')
  }
  else if(i==471){
    sourceA<-c(sourceA,'PAT')
  }
  else if(i>471 & i<510){
    sourceA<-c(sourceA,'COL')
  }
  else if(i>509 & i<576){
    sourceA<-c(sourceA,'PIN')
  }
  
  all_ThreePrime<-c(all_ThreePrime,Gend[length(Gend)])
  all_namesA<-c(all_namesA,names[i])
  all_listA<-data.frame(sourceA,all_namesA,all_ThreePrime) #have to use data.frame because tibble wouldn't allow direct access
  
  #5
  for(u in 2:length(Gstart)){#for all sections aka pairs of del starts & end
    all_Gpos<-c(all_Gpos,Gstart[u],Gend[u],NA) #pair up del start and end, put NA to be distinguished from the next DEL section of this sequence
    archGstart<-c(archGstart,Gstart[u]) #for archive
    archGend<-c(archGend,Gend[u])
    archLength<-c(archLength,Gend[u]-Gstart[u]+1)
  }
  archNames<-c(archNames,rep(names[i],length(Gstart)-1)) #excluding the first deletion
  all_namesB<-c(all_namesB,rep(names[i],(length(Gstart)-1)*3)) #let every listed deletion position has its appropriate sequence name
  all_listB<-data.frame(all_namesB,all_Gpos) #put the positions and their names together into a list
  
}


all_listA<-all_listA%>%arrange(desc(all_ThreePrime))%>%add_column(id=1:length(all_ThreePrime)) #sort listA based on 3' and add id from largest to smallest, so that the longest deletion (farthest 3' del end) would be at the bottom

for(q in 1:length(all_ThreePrime)){ #for the number of all_listA sequences
  for(h in 1:length(all_Gpos)){ #for all the positions in all_listB
    if(all_listA[q,2]==all_listB[h,1]){ #if the sequence name in all_listA matched with that in all_listB
      all_Gid[h]<-c(all_listA[q,4]) #put the all_listA id according to the indices of all_listB
    }
  }
}
all_Gset<-all_listB%>%add_column(all_Gid) #add the id column into all_listB
ID=length(all_ThreePrime)




##################################################################################
all_namesA<-c(all_listA[,2]) #sequence names in the 2nd column of all_listA
all_source<-c(all_listA[,1])
Gid<-c(all_listA[,4])
archGid<-c()
archSource<-c()
ListA<-data.frame(all_source,all_namesA,Gid) #get the sequence names and their corresponding id
ListB<-data.frame(archNames,archGstart,archGend,archLength)
for(e in 1:length(all_namesA)){
  for(f in 1:length(archNames)){
    if(ListA[e,2]==ListB[f,1]){
      archGid[f]<-c(ListA[e,3])
      archSource[f]<-c(ListA[e,1])
    }
  }
}
Archive<-ListB%>%add_column(archGid,.before="archNames")%>%add_column(archSource,.before="archNames")


##################################################################################

AllDel<-ggplot(all_Gset)+geom_line(aes(y=all_Gpos,x=all_Gid))+coord_flip()+labs(title="5'deletions",y="HXB2 position",x="sequence id")+theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank())

AllDel+scale_x_continuous(breaks=seq(0,590,100))+scale_y_continuous(breaks=seq(500,2200,100))


write_xlsx(Archive,"C:\\Users\\jduan\\Archive_hxb2_0405.xlsx")





########################################Color the intact gag

table<-read_xlsx("C:\\Users\\jduan\\Archive_hxb2_0405.xlsx")
SEQ <- read.alignment("5PSD_master_040522_hxb2.fas",format="fasta")

names<-c()
AllSeq<-c()
for(i in 1:SEQ[[1]]){ #for all the sequences
  names<-c(names,SEQ[[2]][i]) #put names into a list
  AllSeq<-c(AllSeq,SEQ[[3]][i]) #put each sequence as an element in a list
}

#####remember to change the sequence titles in the aa file so that they don't contain the "_Gag" in the end
aa <- read.alignment("0405_5psd_GagIntact.aa.fas",format="fasta") #import the file of intact GAG
I_names<-aa[[2]] #get the names of sequences that have intact GAG


GAGrows<-c()

for(e in 1:length(all_namesB)){
  for(f in 2:length(I_names)){ #starting from 2 because we want to exclude HXB2
    if(all_Gset[e,1]==I_names[f]){
      GAGrows<-c(GAGrows,e)
    }
  }
}
GAGset<-all_Gset[GAGrows,]
nonGAGset<-all_Gset[-GAGrows,]

p<-ggplot()+geom_line(data=GAGset,aes(y=all_Gpos,x=all_Gid),color="red")+geom_line(data=nonGAGset,aes(y=all_Gpos,x=all_Gid))+coord_flip()+labs(title="5'deletions with GAG-intact colored",y="HXB2 position",x="sequence id")

p+scale_x_continuous(breaks=seq(0,590,100))+scale_y_continuous(breaks=seq(500,2000,100))+theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),text = element_text(size = 15))   

