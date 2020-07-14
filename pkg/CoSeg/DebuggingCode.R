# Debugging Code for CoSeg

library(openxlsx)
# library(kinship2)

setwd("C:/Users/kamicz/Desktop/CoSeg Development/")
install.packages("CoSeg.tar.gz", repos = NULL, type="source",INSTALL_opts = c('--no-lock'))
library(CoSeg)

temp=read.xlsx("1011_MSH2_2.xlsx")

temp2=FormatWebToCoSeg2(temp)

# temp2$genotype=AnalyzePedigreeGenotypes(temp2)

PlotPedigree(temp2, affected.vector={temp2$affection==1})

RankMembers(temp2, affected.vector={temp2$affection==1})





# testing package
rm(list = ls())
#can create a tarball with 7zip
setwd("C:/Users/kamicz/Desktop/CoSeg Development/")
install.packages("CoSeg.tar.gz", repos = NULL, type="source",INSTALL_opts = c('--no-lock'))
library(CoSeg)

ped=data.frame(degree=c(3,2,2,3,3,1,1,2,2,3), momid=c(3,NA,7,3,3,NA,NA,7,NA,8), dadid=c(2,NA,6,2,2,NA,NA,6,NA,9), id=1:10, age=c(45,60,50,31,41,68,65,55,62,43), female=c(1,0,1,0,1,0,1,1,0,1), y.born=0, dead=0, geno=2, famid=1, bBRCA1.d=0, oBRCA1.d=0, bBRCA1.aoo=NA, oBRCA1.aoo=NA, proband=0)
ped$y.born=2010-ped$age
ped$geno[c(1,3)]=1
ped$bBRCA1.d[c(1,3)]=1
ped$bBRCA1.aoo[1]=45
ped$bBRCA1.aoo[3]=50
ped$proband[1]=1

ped=ped[c(6,7,2,3,8,9,1,4,5,10),]

#Calculate the likelihood ratio
CalculateLikelihoodRatio(ped=ped, affected.vector={ped$bBRCA1.d|ped$oBRCA1.d}, gene="BRCA1")

#Plot the pedigree
PlotPedigree(ped, affected.vector={ped$bBRCA1.d|ped$oBRCA1.d})

#Pedigree Ranking function
RankMembers(ped, affected.vector={ped$bBRCA1.d|ped$oBRCA1.d})


# Fin
