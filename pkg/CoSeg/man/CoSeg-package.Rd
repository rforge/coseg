\name{CoSeg-package}
\alias{CoSeg-package}
\alias{CoSeg}
\docType{package}
\title{
\packageTitle{CoSeg}
}
\description{
\packageDescription{CoSeg}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{CoSeg}
\packageIndices{CoSeg}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
\packageAuthor{CoSeg}

Maintainer: \packageMaintainer{CoSeg}
}
\references{
~~ Literature or other references for background information ~~
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R ~~
%~~ documentation directory ~~
\keyword{ package }
%%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%%}
\examples{
  # #Load all the data included in the CoSeg package.
  # data(BRCA1Frequencies.df, package="CoSeg")
  # data(BRCA2Frequencies.df, package="CoSeg")
  # data(MLH1Frequencies.df, package="CoSeg")
  # data(USDemographics.df, package="CoSeg")
  # data(ChinaDemographics.df, package="CoSeg")

  #summaries of all the data
  str(BRCA1Frequencies.df)
  str(BRCA2Frequencies.df)
  str(MLH1Frequencies.df)
  str(USDemographics.df)
  str(ChinaDemographics.df)

  #Make a tree with no affection status, g=4 generations above, gdown=2 generations below
  #seed.age=50, and demographics.df=NULL which defaults to USDemographics.df.
  tree1=MakeTree()

  #Make a tree using Chinese demographics instead.
  tree2=MakeTree(demographics.df=ChinaDemographics.df)

  #Add affection statust to tree2 using BRCA1Frequencies.df which gives the BRCA1
  #penetrance function
  tree1a=AddAffectedToTree(tree.f=tree1,frequencies.df=BRCA1Frequencies.df)

  #make a tree with affection status (same as running MakeTree() and then AddAffectedToTree())
  tree3=MakeAffectedTrees(n=1,g=2,gdown=2,frequencies.df=MLH1Frequencies.df)
  #tree4=MakeAffectedTrees(n=1,g=2,gdown=2,frequencies.df=BRCA2Frequencies.df)


  #Depending on the size of the pedigree generated, probands (defined here as members of the
  #pedigree who are carriers of the genotype with the disease) may not always be present in
  #the pedigree.  To alleviate this problem in this example we manually generate a pedigree.
  #Note that this is from the Mohammadi paper where the Likelihood method originates from.
  ped=data.frame(degree=c(3,2,2,3,3,1,1,2,2,3), momid=c(3,NA,7,3,3,NA,NA,7,NA,8),
    dadid=c(2,NA,6,2,2,NA,NA,6,NA,9), id=1:10, age=c(45,60,50,31,41,68,65,55,62,43),
    female=c(1,0,1,0,1,0,1,1,0,1), y.born=0, dead=0, geno=2, famid=1, bBRCA1.d=0, oBRCA1.d=0,
    bBRCA1.aoo=NA, oBRCA1.aoo=NA, proband=0)
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

  #Rank and plot the members of the pedigree with unknown genotypes
  RankMembers(ped=ped, affected.vector={ped$bBRCA1.d|ped$oBRCA1.d}, gene="BRCA1")
}
