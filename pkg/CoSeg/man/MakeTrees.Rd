\name{MakeTrees}
\alias{MakeTrees}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make a tree
}
\description{
This function makes a tree based on historical average age of marriage for male and female, age of death for both sexes, and number of offspring.  For convenience, historical data for the US and China are included in the package.
}
\usage{
MakeTrees(n=1, g = 4, gdown = 2, demographics.df = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
  An integer giving the number of trees to make.
  }
  \item{g}{
    An integer giving the number of ancestral generations of the proband to simulate.
}
  \item{gdown}{
    An integer giving the number of descendent generations of the proband to simulate.
}
  \item{demographics.df}{
    A dataframe giving the demographics to use.  This defaults to \code{\link{USDemographics.df}} which gives historic data on US demographics.  Alternatively, the user can choose \code{\link{ChinaDemographics.df}} or input their own dataframe.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Returns a dataframe for the generated tree.  The dataframe has the columns listed below.
%%  If it is a LIST, use
  \item{degree}{An integer giving the generation number of the individual where ancestors have lower values than descendents.}
  \item{momid}{A numeric giving the ID of the mother of the individual or NA if the individual is a founder.}
  \item{dadid}{A numeric giving the ID of the father of the individual or NA if the individual is a founder.}
  \item{id}{A numeric label for the individual.}
  \item{age}{A numeric giving the current age of the individual.}
  \item{female}{A boolean/integer depicting whether the individual is female(1) or not(0).}
  \item{y.born}{An integer depicting the year the individual was born.}
  \item{geno}{An boolean/integer giving telling whether the individual is a carrier(1) or not(0).}
%% ...
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
John Michael O. Ranola and Brian H. Shirts
}
\note{
  Note that CRAN throws up a note when this package is checked saying "no visible binding for global variable...".  This is normal and does not affect the analysis.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{AddAffectedToTree}}, \code{\link{MakeAffectedTrees}}, and \code{\link{PlotPedigree}}
%\code{\link{BRCA1Frequencies.df}}, \code{\link{BRCA2Frequencies.df}}, and \code{\link{MLH1Frequencies.df}}

}
\examples{

  \dontrun{

    #summaries of all the data
    str(BRCA1Frequencies.df)
    str(BRCA2Frequencies.df)
    str(MLH1Frequencies.df)
    str(USDemographics.df)
    str(ChinaDemographics.df)

    #Make a tree with no affection status, g=4 generations above, gdown=2 generations below,
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
