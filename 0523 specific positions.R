######looking at specific positions
setwd("C:/Users/jduan")
library(seqinr)
library(tidyverse)
library(ggplot2)
library(dplyr)
library("readxl")
library("writexl")

Archive<-read_excel("5PSD_master_041922_tilda_uniq_allDel_cleared.xlsx")

#1 get all the sequences
SEQ <- read.alignment("5PSD_master_040422_tilda_uniq.fas",format="fasta") #import the file
names<-c()
AllSeq<-c()
for(i in 1:SEQ[[1]]){ #for all the sequences
  names<-c(names,SEQ[[2]][i]) #put names into a list
  AllSeq<-c(AllSeq,SEQ[[3]][i]) #put each sequence as an element in a list
}
HXB2<-c()
for(n in 1:2157){
  HXB2<-c(HXB2,substr(SEQ[[3]][1],n,n)) #get the HXB2 sequences
}

####################Three-way junction 1
#2
status<-c()

for(i in 2:SEQ[[1]]){
  seq<-c(AllSeq[i]) #each round of loop, the individual seq is a new one
  nseq<-c() #extracting the sequence string into individual nucleotides
  for(n in 1:2157){ #here, 2157 can be changed to other numbers according to how much we like
    nseq<-c(nseq,substr(AllSeq[i],n,n))
  }
  str<-paste(nseq[765],nseq[766],nseq[767],nseq[768],nseq[769],nseq[770],nseq[771],sep="",collapse="")
  c=str_count(str,"tttt")
  if(c>0){
    status<-c(status,"present")
  }
  else{
    status<-c(status,"absent")
  }
}
diff_3WJ1<-data.frame(names[2:575],status)




#####################Packaging Signal
#2
SeqName<-c()
NumDiff<-c()
NumGap<-c()
NumNG<-c()
status<-c()
allpos<-c()
alldiff<-c()
for(i in 2:SEQ[[1]]){
  seq<-c(AllSeq[i]) #each round of loop, the individual seq is a new one
  nseq<-c() #extracting the sequence string into individual nucleotides
  for(n in 1:2157){ #here, 2000 can be changed to other numbers according to how much we like
    nseq<-c(nseq,substr(AllSeq[i],n,n))
  }
  
  c=0 #count the number of differences from zero
  Gap=0 #count the number of absences from zero
  position<-c() #for every seq, there's a new start of recording for the positions that are different
  difference<-c()
  for(m in 774:787){ #for the length of the packaging signal
    if(HXB2[m]!=nseq[m]){ #if it's different from the HXB2 position
      c<-c+1 #then we add one
      position<-c(position,m) #archive the position that is different
      difference<-c(difference,nseq[m]) #archive what's different about that position
      if(nseq[m]=="~"){ #if what's making the difference is a gap
        Gap<-Gap+1 #add it as absence
      }
    }
    else{ #if the same
      c=c #then we don't add
    }
  }
  NumGap[i-1]<-Gap
  NumNG[i-1]<-c-Gap
  NumDiff[i-1]<-c #get the number of differences; use i-1 because we start from 2, and we want the table we make later starts from 1
  if(NumDiff[i-1]==0){
    status[i-1]<-'same'
  }
  else if(NumDiff[i-1]==14) { #14 is the length of the sequence for packaging signal
    status[i-1]<-'different'
  }
  else{
    status[i-1]<-'partial'
  }
  allpos<-c(allpos,position) #put the position vector for this sequence as an element into a list called allpos for archive later
  alldiff<-c(alldiff,difference)
  SeqName<-c(SeqName,rep(names[i],length(position))) #repeat the sequence's name for the number of times that equals the number of positions that are different, so that we can archive them together later
}

diff_Psi<-data.frame(names[2:575],status,NumDiff,NumGap,NumNG)
archPsi<-data.frame(SeqName,allpos,alldiff)




####################Three-way junction 2
#2
SeqName<-c()
NumDiff<-c()
NumGap<-c()
NumNG<-c()
status<-c()
allpos<-c()
alldiff<-c()
for(i in 2:SEQ[[1]]){
  seq<-c(AllSeq[i]) #each round of loop, the individual seq is a new one
  nseq<-c() #extracting the sequence string into individual nucleotides
  for(n in 1:2157){ #here, 2000 can be changed to other numbers according to how much we like
    nseq<-c(nseq,substr(AllSeq[i],n,n))
  }
  
  c=0 #count the number of differences from zero
  Gap=0 #count the number of absences from zero
  position<-c() #for every seq, there's a new start of recording for the positions that are different
  difference<-c()
  for(m in 790:793){ #for the length of the 3WJ.2
    if(HXB2[m]!=nseq[m]){ #if it's different from the HXB2 position
      c<-c+1 #then we add one
      position<-c(position,m) #archive the position that is different
      difference<-c(difference,nseq[m]) #archive what's different about that position
      if(nseq[m]=="~"){ #if what's making the difference is a gap
        Gap<-Gap+1 #add it as absence
      }
    }
    else{ #if the same
      c=c #then we don't add
    }
  }
  NumGap[i-1]<-Gap
  NumNG[i-1]<-c-Gap
  NumDiff[i-1]<-c #get the number of differences; use i-1 because we start from 2, and we want the table we make later starts from 1
  if(NumDiff[i-1]==0){
    status[i-1]<-'same'
  }
  else if(NumDiff[i-1]==4) { #4 is the length of the sequence for 3WJ.2
    status[i-1]<-'different'
  }
  else{
    status[i-1]<-'partial'
  }
  allpos<-c(allpos,position) #put the position vector for this sequence as an element into a list called allpos for archive later
  alldiff<-c(alldiff,difference)
  SeqName<-c(SeqName,rep(names[i],length(position))) #repeat the sequence's name for the number of times that equals the number of positions that are different, so that we can archive them together later
}

diff_3WJ2<-data.frame(names[2:575],status,NumDiff,NumGap,NumNG)
arch3WJ2<-data.frame(SeqName,allpos,alldiff)




#####################IPDA PSI probe
#2
SeqName<-c()
NumDiff<-c()
NumGap<-c()
NumNG<-c()
status<-c()
allpos<-c()
alldiff<-c()
for(i in 2:SEQ[[1]]){
  seq<-c(AllSeq[i]) #each round of loop, the individual seq is a new one
  nseq<-c() #extracting the sequence string into individual nucleotides
  for(n in 1:2157){ #here, 2000 can be changed to other numbers according to how much we like
    nseq<-c(nseq,substr(AllSeq[i],n,n))
  }
  
  c=0 #count the number of differences from zero
  Gap=0 #count the number of absences from zero
  position<-c() #for every seq, there's a new start of recording for the positions that are different
  difference<-c()
  for(m in 745:763){ #for the length of the IPDA probe
    if(HXB2[m]!=nseq[m]){ #if it's different from the HXB2 position
      c<-c+1 #then we add one
      position<-c(position,m) #archive the position that is different
      difference<-c(difference,nseq[m]) #archive what's different about that position
      if(nseq[m]=="~"){ #if what's making the difference is a gap
        Gap<-Gap+1 #add it as absence
      }
    }
    else{ #if the same
      c=c #then we don't add
    }
  }
  NumGap[i-1]<-Gap
  NumNG[i-1]<-c-Gap
  NumDiff[i-1]<-c #get the number of differences; use i-1 because we start from 2, and we want the table we make later starts from 1
  if(NumDiff[i-1]==0){
    status[i-1]<-'same'
  }
  else if(NumDiff[i-1]==19) { #19 is the length of the sequence for IPDA probe
    status[i-1]<-'different'
  }
  else{
    status[i-1]<-'partial'
  }
  allpos<-c(allpos,position) #put the position vector for this sequence as an element into a list called allpos for archive later
  alldiff<-c(alldiff,difference)
  SeqName<-c(SeqName,rep(names[i],length(position))) #repeat the sequence's name for the number of times that equals the number of positions that are different, so that we can archive them together later
}

diff_IPDAprobe<-data.frame(names[2:575],status,NumDiff,NumGap,NumNG)
archIPDAprobe<-data.frame(SeqName,allpos,alldiff)



#####################IPDA PSI fwd primer
#2
SeqName<-c()
NumDiff<-c()
NumGap<-c()
NumNG<-c()
status<-c()
allpos<-c()
alldiff<-c()
for(i in 2:SEQ[[1]]){
  seq<-c(AllSeq[i]) #each round of loop, the individual seq is a new one
  nseq<-c() #extracting the sequence string into individual nucleotides
  for(n in 1:2157){ #here, 2000 can be changed to other numbers according to how much we like
    nseq<-c(nseq,substr(AllSeq[i],n,n))
  }
  
  c=0 #count the number of differences from zero
  Gap=0 #count the number of absences from zero
  position<-c() #for every seq, there's a new start of recording for the positions that are different
  difference<-c()
  for(m in 695:716){ #for the length of the IPDA fwd primer
    if(HXB2[m]!=nseq[m]){ #if it's different from the HXB2 position
      c<-c+1 #then we add one
      position<-c(position,m) #archive the position that is different
      difference<-c(difference,nseq[m]) #archive what's different about that position
      if(nseq[m]=="~"){ #if what's making the difference is a gap
        Gap<-Gap+1 #add it as absence
      }
    }
    else{ #if the same
      c=c #then we don't add
    }
  }
  NumGap[i-1]<-Gap
  NumNG[i-1]<-c-Gap
  NumDiff[i-1]<-c #get the number of differences; use i-1 because we start from 2, and we want the table we make later starts from 1
  if(NumDiff[i-1]==0){
    status[i-1]<-'same'
  }
  else if(NumDiff[i-1]==22) { #29 is the length of the IPDA fwd primer
    status[i-1]<-'different'
  }
  else{
    status[i-1]<-'partial'
  }
  allpos<-c(allpos,position) #put the position vector for this sequence as an element into a list called allpos for archive later
  alldiff<-c(alldiff,difference)
  SeqName<-c(SeqName,rep(names[i],length(position))) #repeat the sequence's name for the number of times that equals the number of positions that are different, so that we can archive them together later
}

diff_IPDAfwd<-data.frame(names[2:575],status,NumDiff,NumGap,NumNG)
archIPDAfwd<-data.frame(SeqName,allpos,alldiff)


#####################DIS stem 1
#2
SeqName<-c()
NumDiff<-c()
NumGap<-c()
NumNG<-c()
status<-c()
allpos<-c()
alldiff<-c()
for(i in 2:SEQ[[1]]){
  seq<-c(AllSeq[i]) #each round of loop, the individual seq is a new one
  nseq<-c() #extracting the sequence string into individual nucleotides
  for(n in 1:2157){ #here, 2000 can be changed to other numbers according to how much we like
    nseq<-c(nseq,substr(AllSeq[i],n,n))
  }
  
  c=0 #count the number of differences from zero
  Gap=0 #count the number of absences from zero
  position<-c() #for every seq, there's a new start of recording for the positions that are different
  difference<-c()
  for(m in 699:701){ #for the length of DIS stem 1
    if(HXB2[m]!=nseq[m]){ #if it's different from the HXB2 position
      c<-c+1 #then we add one
      position<-c(position,m) #archive the position that is different
      difference<-c(difference,nseq[m]) #archive what's different about that position
      if(nseq[m]=="~"){ #if what's making the difference is a gap
        Gap<-Gap+1 #add it as absence
      }
    }
    else{ #if the same
      c=c #then we don't add
    }
  }
  NumGap[i-1]<-Gap
  NumNG[i-1]<-c-Gap
  NumDiff[i-1]<-c #get the number of differences; use i-1 because we start from 2, and we want the table we make later starts from 1
  if(NumDiff[i-1]==0){
    status[i-1]<-'same'
  }
  else if(NumDiff[i-1]==3) { #3 is the length of DIS stem 1
    status[i-1]<-'different'
  }
  else{
    status[i-1]<-'partial'
  }
  allpos<-c(allpos,position) #put the position vector for this sequence as an element into a list called allpos for archive later
  alldiff<-c(alldiff,difference)
  SeqName<-c(SeqName,rep(names[i],length(position))) #repeat the sequence's name for the number of times that equals the number of positions that are different, so that we can archive them together later
}

diff_DISstem1<-data.frame(names[2:575],status,NumDiff,NumGap,NumNG)
archDISstem1<-data.frame(SeqName,allpos,alldiff)


#####################DIS Palindrome
#2
SeqName<-c()
NumDiff<-c()
NumGap<-c()
NumNG<-c()
status<-c()
allpos<-c()
alldiff<-c()
for(i in 2:SEQ[[1]]){
  seq<-c(AllSeq[i]) #each round of loop, the individual seq is a new one
  nseq<-c() #extracting the sequence string into individual nucleotides
  for(n in 1:2157){ #here, 2000 can be changed to other numbers according to how much we like
    nseq<-c(nseq,substr(AllSeq[i],n,n))
  }
  
  c=0 #count the number of differences from zero
  Gap=0 #count the number of absences from zero
  position<-c() #for every seq, there's a new start of recording for the positions that are different
  difference<-c()
  for(m in 716:721){ #for the length of the DIS palindrome
    if(HXB2[m]!=nseq[m]){ #if it's different from the HXB2 position
      c<-c+1 #then we add one
      position<-c(position,m) #archive the position that is different
      difference<-c(difference,nseq[m]) #archive what's different about that position
      if(nseq[m]=="~"){ #if what's making the difference is a gap
        Gap<-Gap+1 #add it as absence
      }
    }
    else{ #if the same
      c=c #then we don't add
    }
  }
  NumGap[i-1]<-Gap
  NumNG[i-1]<-c-Gap
  NumDiff[i-1]<-c #get the number of differences; use i-1 because we start from 2, and we want the table we make later starts from 1
  if(NumDiff[i-1]==0){
    status[i-1]<-'same'
  }
  else if(NumDiff[i-1]==6) { #6 is the length of DIS palindrome
    status[i-1]<-'different'
  }
  else{
    status[i-1]<-'partial'
  }
  allpos<-c(allpos,position) #put the position vector for this sequence as an element into a list called allpos for archive later
  alldiff<-c(alldiff,difference)
  SeqName<-c(SeqName,rep(names[i],length(position))) #repeat the sequence's name for the number of times that equals the number of positions that are different, so that we can archive them together later
}

diff_DISpal<-data.frame(names[2:575],status,NumDiff,NumGap,NumNG)
archDISpal<-data.frame(SeqName,allpos,alldiff)


#####################DIS stem 2
#2
SeqName<-c()
NumDiff<-c()
NumGap<-c()
NumNG<-c()
status<-c()
allpos<-c()
alldiff<-c()
for(i in 2:SEQ[[1]]){
  seq<-c(AllSeq[i]) #each round of loop, the individual seq is a new one
  nseq<-c() #extracting the sequence string into individual nucleotides
  for(n in 1:2157){ #here, 2000 can be changed to other numbers according to how much we like
    nseq<-c(nseq,substr(AllSeq[i],n,n))
  }
  
  c=0 #count the number of differences from zero
  Gap=0 #count the number of absences from zero
  position<-c() #for every seq, there's a new start of recording for the positions that are different
  difference<-c()
  for(m in 736:738){ #for the length of the DIS stem 2
    if(HXB2[m]!=nseq[m]){ #if it's different from the HXB2 position
      c<-c+1 #then we add one
      position<-c(position,m) #archive the position that is different
      difference<-c(difference,nseq[m]) #archive what's different about that position
      if(nseq[m]=="~"){ #if what's making the difference is a gap
        Gap<-Gap+1 #add it as absence
      }
    }
    else{ #if the same
      c=c #then we don't add
    }
  }
  NumGap[i-1]<-Gap
  NumNG[i-1]<-c-Gap
  NumDiff[i-1]<-c #get the number of differences; use i-1 because we start from 2, and we want the table we make later starts from 1
  if(NumDiff[i-1]==0){
    status[i-1]<-'same'
  }
  else if(NumDiff[i-1]==3) { #3 is the length of DIS stem 2
    status[i-1]<-'different'
  }
  else{
    status[i-1]<-'partial'
  }
  allpos<-c(allpos,position) #put the position vector for this sequence as an element into a list called allpos for archive later
  alldiff<-c(alldiff,difference)
  SeqName<-c(SeqName,rep(names[i],length(position))) #repeat the sequence's name for the number of times that equals the number of positions that are different, so that we can archive them together later
}

diff_DISstem2<-data.frame(names[2:575],status,NumDiff,NumGap,NumNG)
archDISstem2<-data.frame(SeqName,allpos,alldiff)



#####################DNOVEL1
#2
SeqName<-c()
NumDiff<-c()
NumGap<-c()
NumNG<-c()
status<-c()
allpos<-c()
alldiff<-c()
for(i in 2:SEQ[[1]]){
  seq<-c(AllSeq[i]) #each round of loop, the individual seq is a new one
  nseq<-c() #extracting the sequence string into individual nucleotides
  for(n in 1:2157){ #here, 2000 can be changed to other numbers according to how much we like
    nseq<-c(nseq,substr(AllSeq[i],n,n))
  }
  
  c=0 #count the number of differences from zero
  Gap=0 #count the number of absences from zero
  position<-c() #for every seq, there's a new start of recording for the positions that are different
  difference<-c()
  for(m in 722:725){ #for the length of the DNOVEL1
    if(HXB2[m]!=nseq[m]){ #if it's different from the HXB2 position
      c<-c+1 #then we add one
      position<-c(position,m) #archive the position that is different
      difference<-c(difference,nseq[m]) #archive what's different about that position
      if(nseq[m]=="~"){ #if what's making the difference is a gap
        Gap<-Gap+1 #add it as absence
      }
    }
    else{ #if the same
      c=c #then we don't add
    }
  }
  NumGap[i-1]<-Gap
  NumNG[i-1]<-c-Gap
  NumDiff[i-1]<-c #get the number of differences; use i-1 because we start from 2, and we want the table we make later starts from 1
  if(NumDiff[i-1]==0){
    status[i-1]<-'same'
  }
  else if(NumDiff[i-1]==4) { #4 is the length of DNOVEL1
    status[i-1]<-'different'
  }
  else{
    status[i-1]<-'partial'
  }
  allpos<-c(allpos,position) #put the position vector for this sequence as an element into a list called allpos for archive later
  alldiff<-c(alldiff,difference)
  SeqName<-c(SeqName,rep(names[i],length(position))) #repeat the sequence's name for the number of times that equals the number of positions that are different, so that we can archive them together later
}

diff_DNOVEL1<-data.frame(names[2:575],status,NumDiff,NumGap,NumNG)
archDNOVEL1<-data.frame(SeqName,allpos,alldiff)


#####################D1b
#2
SeqName<-c()
NumDiff<-c()
NumGap<-c()
NumNG<-c()
status<-c()
allpos<-c()
alldiff<-c()
for(i in 2:SEQ[[1]]){
  seq<-c(AllSeq[i]) #each round of loop, the individual seq is a new one
  nseq<-c() #extracting the sequence string into individual nucleotides
  for(n in 1:2157){ #here, 2000 can be changed to other numbers according to how much we like
    nseq<-c(nseq,substr(AllSeq[i],n,n))
  }
  
  c=0 #count the number of differences from zero
  Gap=0 #count the number of absences from zero
  position<-c() #for every seq, there's a new start of recording for the positions that are different
  difference<-c()
  for(m in 730:733){ #for the length of D1b
    if(HXB2[m]!=nseq[m]){ #if it's different from the HXB2 position
      c<-c+1 #then we add one
      position<-c(position,m) #archive the position that is different
      difference<-c(difference,nseq[m]) #archive what's different about that position
      if(nseq[m]=="~"){ #if what's making the difference is a gap
        Gap<-Gap+1 #add it as absence
      }
    }
    else{ #if the same
      c=c #then we don't add
    }
  }
  NumGap[i-1]<-Gap
  NumNG[i-1]<-c-Gap
  NumDiff[i-1]<-c #get the number of differences; use i-1 because we start from 2, and we want the table we make later starts from 1
  if(NumDiff[i-1]==0){
    status[i-1]<-'same'
  }
  else if(NumDiff[i-1]==4) { #4 is the length of D1b
    status[i-1]<-'different'
  }
  else{
    status[i-1]<-'partial'
  }
  allpos<-c(allpos,position) #put the position vector for this sequence as an element into a list called allpos for archive later
  alldiff<-c(alldiff,difference)
  SeqName<-c(SeqName,rep(names[i],length(position))) #repeat the sequence's name for the number of times that equals the number of positions that are different, so that we can archive them together later
}

diff_D1B<-data.frame(names[2:575],status,NumDiff,NumGap,NumNG)
archD1B<-data.frame(SeqName,allpos,alldiff)


#####################D1c
#2
SeqName<-c()
NumDiff<-c()
NumGap<-c()
NumNG<-c()
status<-c()
allpos<-c()
alldiff<-c()
for(i in 2:SEQ[[1]]){
  seq<-c(AllSeq[i]) #each round of loop, the individual seq is a new one
  nseq<-c() #extracting the sequence string into individual nucleotides
  for(n in 1:2157){ #here, 2000 can be changed to other numbers according to how much we like
    nseq<-c(nseq,substr(AllSeq[i],n,n))
  }
  
  c=0 #count the number of differences from zero
  Gap=0 #count the number of absences from zero
  position<-c() #for every seq, there's a new start of recording for the positions that are different
  difference<-c()
  for(m in 751:754){ #for the length of D1c
    if(HXB2[m]!=nseq[m]){ #if it's different from the HXB2 position
      c<-c+1 #then we add one
      position<-c(position,m) #archive the position that is different
      difference<-c(difference,nseq[m]) #archive what's different about that position
      if(nseq[m]=="~"){ #if what's making the difference is a gap
        Gap<-Gap+1 #add it as absence
      }
    }
    else{ #if the same
      c=c #then we don't add
    }
  }
  NumGap[i-1]<-Gap
  NumNG[i-1]<-c-Gap
  NumDiff[i-1]<-c #get the number of differences; use i-1 because we start from 2, and we want the table we make later starts from 1
  if(NumDiff[i-1]==0){
    status[i-1]<-'same'
  }
  else if(NumDiff[i-1]==4) { #4 is the length of D1c
    status[i-1]<-'different'
  }
  else{
    status[i-1]<-'partial'
  }
  allpos<-c(allpos,position) #put the position vector for this sequence as an element into a list called allpos for archive later
  alldiff<-c(alldiff,difference)
  SeqName<-c(SeqName,rep(names[i],length(position))) #repeat the sequence's name for the number of times that equals the number of positions that are different, so that we can archive them together later
}

diff_D1c<-data.frame(names[2:575],status,NumDiff,NumGap,NumNG)
archD1c<-data.frame(SeqName,allpos,alldiff)

#####################IPDA PSI reverse primer
#2
SeqName<-c()
NumDiff<-c()
NumGap<-c()
NumNG<-c()
status<-c()
allpos<-c()
alldiff<-c()
for(i in 2:SEQ[[1]]){
  seq<-c(AllSeq[i]) #each round of loop, the individual seq is a new one
  nseq<-c() #extracting the sequence string into individual nucleotides
  for(n in 1:2157){ #here, 2000 can be changed to other numbers according to how much we like
    nseq<-c(nseq,substr(AllSeq[i],n,n))
  }
  
  c=0 #count the number of differences from zero
  Gap=0 #count the number of absences from zero
  position<-c() #for every seq, there's a new start of recording for the positions that are different
  difference<-c()
  for(m in 783:805){ #for the length of IPDA rev primer
    if(HXB2[m]!=nseq[m]){ #if it's different from the HXB2 position
      c<-c+1 #then we add one
      position<-c(position,m) #archive the position that is different
      difference<-c(difference,nseq[m]) #archive what's different about that position
      if(nseq[m]=="~"){ #if what's making the difference is a gap
        Gap<-Gap+1 #add it as absence
      }
    }
    else{ #if the same
      c=c #then we don't add
    }
  }
  NumGap[i-1]<-Gap
  NumNG[i-1]<-c-Gap
  NumDiff[i-1]<-c #get the number of differences; use i-1 because we start from 2, and we want the table we make later starts from 1
  if(NumDiff[i-1]==0){
    status[i-1]<-'same'
  }
  else if(NumDiff[i-1]==23) { #23 is the length of IPDA rev primer
    status[i-1]<-'different'
  }
  else{
    status[i-1]<-'partial'
  }
  allpos<-c(allpos,position) #put the position vector for this sequence as an element into a list called allpos for archive later
  alldiff<-c(alldiff,difference)
  SeqName<-c(SeqName,rep(names[i],length(position))) #repeat the sequence's name for the number of times that equals the number of positions that are different, so that we can archive them together later
}

diff_IPDArev<-data.frame(names[2:575],status,NumDiff,NumGap,NumNG)
archIPDArev<-data.frame(SeqName,allpos,alldiff)

#####################DNOVEL2
#2
SeqName<-c()
NumDiff<-c()
NumGap<-c()
NumNG<-c()
status<-c()
allpos<-c()
alldiff<-c()
for(i in 2:SEQ[[1]]){
  seq<-c(AllSeq[i]) #each round of loop, the individual seq is a new one
  nseq<-c() #extracting the sequence string into individual nucleotides
  for(n in 1:2157){ #here, 2000 can be changed to other numbers according to how much we like
    nseq<-c(nseq,substr(AllSeq[i],n,n))
  }
  
  c=0 #count the number of differences from zero
  Gap=0 #count the number of absences from zero
  position<-c() #for every seq, there's a new start of recording for the positions that are different
  difference<-c()
  for(m in 798:803){ #for the length of DNOVEL2
    if(HXB2[m]!=nseq[m]){ #if it's different from the HXB2 position
      c<-c+1 #then we add one
      position<-c(position,m) #archive the position that is different
      difference<-c(difference,nseq[m]) #archive what's different about that position
      if(nseq[m]=="~"){ #if what's making the difference is a gap
        Gap<-Gap+1 #add it as absence
      }
    }
    else{ #if the same
      c=c #then we don't add
    }
  }
  NumGap[i-1]<-Gap
  NumNG[i-1]<-c-Gap
  NumDiff[i-1]<-c #get the number of differences; use i-1 because we start from 2, and we want the table we make later starts from 1
  if(NumDiff[i-1]==0){
    status[i-1]<-'same'
  }
  else if(NumDiff[i-1]==6) { #6 is the length of DNOVEL2
    status[i-1]<-'different'
  }
  else{
    status[i-1]<-'partial'
  }
  allpos<-c(allpos,position) #put the position vector for this sequence as an element into a list called allpos for archive later
  alldiff<-c(alldiff,difference)
  SeqName<-c(SeqName,rep(names[i],length(position))) #repeat the sequence's name for the number of times that equals the number of positions that are different, so that we can archive them together later
}

diff_DNOVEL2<-data.frame(names[2:575],status,NumDiff,NumGap,NumNG)
archDNOVEL2<-data.frame(SeqName,allpos,alldiff)

#####################GAGstart
#2
SeqName<-c()
NumDiff<-c()
NumGap<-c()
NumNG<-c()
status<-c()
allpos<-c()
alldiff<-c()
for(i in 2:SEQ[[1]]){
  seq<-c(AllSeq[i]) #each round of loop, the individual seq is a new one
  nseq<-c() #extracting the sequence string into individual nucleotides
  for(n in 1:2157){ #here, 2000 can be changed to other numbers according to how much we like
    nseq<-c(nseq,substr(AllSeq[i],n,n))
  }
  
  c=0 #count the number of differences from zero
  Gap=0 #count the number of absences from zero
  position<-c() #for every seq, there's a new start of recording for the positions that are different
  difference<-c()
  for(m in 798:800){ #for the length of GAGstart
    if(HXB2[m]!=nseq[m]){ #if it's different from the HXB2 position
      c<-c+1 #then we add one
      position<-c(position,m) #archive the position that is different
      difference<-c(difference,nseq[m]) #archive what's different about that position
      if(nseq[m]=="~"){ #if what's making the difference is a gap
        Gap<-Gap+1 #add it as absence
      }
    }
    else{ #if the same
      c=c #then we don't add
    }
  }
  NumGap[i-1]<-Gap
  NumNG[i-1]<-c-Gap
  NumDiff[i-1]<-c #get the number of differences; use i-1 because we start from 2, and we want the table we make later starts from 1
  if(NumDiff[i-1]==0){
    status[i-1]<-'same'
  }
  else if(NumDiff[i-1]==3) { #3 is the length of GAGstart
    status[i-1]<-'different'
  }
  else{
    status[i-1]<-'partial'
  }
  allpos<-c(allpos,position) #put the position vector for this sequence as an element into a list called allpos for archive later
  alldiff<-c(alldiff,difference)
  SeqName<-c(SeqName,rep(names[i],length(position))) #repeat the sequence's name for the number of times that equals the number of positions that are different, so that we can archive them together later
}

diff_GAGstart<-data.frame(names[2:575],status,NumDiff,NumGap,NumNG)
archGAGstart<-data.frame(SeqName,allpos,alldiff)

#####################MSD
#2
SeqName<-c()
NumDiff<-c()
NumGap<-c()
NumNG<-c()
status<-c()
allpos<-c()
alldiff<-c()
for(i in 2:SEQ[[1]]){
  seq<-c(AllSeq[i]) #each round of loop, the individual seq is a new one
  nseq<-c() #extracting the sequence string into individual nucleotides
  for(n in 1:2157){ #here, 2000 can be changed to other numbers according to how much we like
    nseq<-c(nseq,substr(AllSeq[i],n,n))
  }
  
  c=0 #count the number of differences from zero
  Gap=0 #count the number of absences from zero
  position<-c() #for every seq, there's a new start of recording for the positions that are different
  difference<-c()
  for(m in 747:750){ #for the length of MSD
    if(HXB2[m]!=nseq[m]){ #if it's different from the HXB2 position
      c<-c+1 #then we add one
      position<-c(position,m) #archive the position that is different
      difference<-c(difference,nseq[m]) #archive what's different about that position
      if(nseq[m]=="~"){ #if what's making the difference is a gap
        Gap<-Gap+1 #add it as absence
      }
    }
    else{ #if the same
      c=c #then we don't add
    }
  }
  NumGap[i-1]<-Gap
  NumNG[i-1]<-c-Gap
  NumDiff[i-1]<-c #get the number of differences; use i-1 because we start from 2, and we want the table we make later starts from 1
  if(NumDiff[i-1]==0){
    status[i-1]<-'same'
  }
  else if(NumDiff[i-1]==4) { #4 is the length of MSD
    status[i-1]<-'different'
  }
  else{
    status[i-1]<-'partial'
  }
  allpos<-c(allpos,position) #put the position vector for this sequence as an element into a list called allpos for archive later
  alldiff<-c(alldiff,difference)
  SeqName<-c(SeqName,rep(names[i],length(position))) #repeat the sequence's name for the number of times that equals the number of positions that are different, so that we can archive them together later
}

diff_MSD<-data.frame(names[2:575],status,NumDiff,NumGap,NumNG)
archMSD<-data.frame(SeqName,allpos,alldiff)


###############################################

SpecPos<-data.frame(diff_IPDAfwd,diff_IPDAprobe,diff_IPDArev,diff_DISstem1,diff_DISpal,diff_DISstem2,diff_MSD,diff_D1B,diff_D1c,diff_DNOVEL1,diff_DNOVEL2,diff_3WJ1,diff_Psi,diff_3WJ2,diff_GAGstart)

write_xlsx(SpecPos,"C:\\Users\\jduan\\0524SpecPos.xlsx")




