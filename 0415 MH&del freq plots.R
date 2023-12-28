setwd("C:/Users/jduan")
library(seqinr)
library(tidyverse)
library(ggplot2)
library(dplyr)
library("readxl")
library("writexl")

########################microhomology line plot
microhomology<-read_xlsx("5PSD_master_041422_tilda_uniq_consensus_WithSpacer_allMH_updated.xlsx")

MH_hist_info<-hist(microhomology$M_hxb2delend,breaks=seq(590,2100,10))
MH_counts<-c(MH_hist_info$counts)
MH_freq<-MH_counts/574 #divided by number of sequences

Consensus_coordinates=seq(590,2090,10) #should I write HXB2 or consensus coordinates?

data<-data.frame(Consensus_coordinates,MH_freq)
ggplot(data,aes(Consensus_coordinates))+geom_line(aes(y=MH_freq,color="MH freq"),size=1,col="coral2")+scale_x_continuous(breaks=seq(590,2100,by=100))+ylim(0,0.2)+theme_classic()

#######################del freq according to HXB2_Delend (3'end)
archive <- read_excel("5PSD_master_041422_tilda_uniq_allDel_cleared.xlsx") #get the archive of deletion junctions
del_hist_info<-hist(archive$HXB2_Delend,breaks=seq(590,2100,10))
del_3end_counts<-c(del_hist_info$counts)
del_3end_freq<-del_3end_counts/574

Consensus_coordinates=seq(590,2090,10)

data<-data.frame(Consensus_coordinates,del_3end_freq)
ggplot(data,aes(Consensus_coordinates))+geom_line(aes(y=del_3end_freq,color="del freq"),size=1,col="skyblue4")+scale_x_continuous(breaks=seq(590,2100,by=100))+ylim(0,0.2)+theme_classic()



#######################del freq according to HXB2_Delstart (5'end)
del_hist_info<-hist(archive$HXB2_DelStart,breaks=seq(590,2100,10))
del_5end_counts<-c(del_hist_info$counts)
del_5end_freq<-del_5end_counts/574

Consensus_coordinates=seq(590,2090,10)

data<-data.frame(Consensus_coordinates,del_3end_freq)
ggplot(data,aes(Consensus_coordinates))+geom_line(aes(y=del_5end_freq,color="del freq"),size=1,col="grey55")+scale_x_continuous(breaks=seq(590,2100,by=100))+ylim(0,0.2)+theme_classic()



#calculate the frequency bin size 50
#make a histogram
identity<-c(rep(c("del"),length(archive$HXB2_Delend)))
identity<-c(identity,rep(c("MH"),length(microhomology$M_hxb2delend)))
df<-data.frame(identity,count=c(archive$HXB2_Delend,microhomology$M_hxb2delend))
ggplot(df,aes(x=count,fill=identity,color=identity))+geom_histogram(binwidth=50,position="identity",alpha=0.5)+geom_density(col=2)+theme_classic()
hist(archive$HXB2_Delend,main="deletion vs microhomology",xlab="3'consensus coordinates",ylab="number of sequences",xlim=c(600,2100),ylim=c(0,325),breaks=30,col="skyblue3",labels=T,cex.lab=1.3)
hist(microhomology$M_hxb2delend,add=TRUE,col="skyblue1",breaks=30,labels=T,cex.lab=1.3)
#curve(dnorm(x,mean=mean(archive$HXB2_Delend),sd=sd(archive$HXB2_Delend),add=TRUE,col="yellow"))
#legend("topright",c("deletions","microhomology"),fill=c("skyblue3","skyblue1"),bty="n")

