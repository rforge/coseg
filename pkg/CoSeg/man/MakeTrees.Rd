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
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{AddAffectedToTree}}, \code{\link{MakeAffectedTrees}}, and \code{\link{PlotPedigree}}
%\code{\link{BRCA1Frequencies.df}}, \code{\link{BRCA2Frequencies.df}}, and \code{\link{MLH1Frequencies.df}}

}
\examples{

  \dontrun{
    #Load all the data included in the CoSeg package.
    data(BRCA1Frequencies.df, package="CoSeg")
    data(BRCA2Frequencies.df, package="CoSeg")
    data(MLH1Frequencies.df, package="CoSeg")
    data(USDemographics.df, package="CoSeg")
    data(ChinaDemographics.df, package="CoSeg")

    #summaries of all the data
    str(BRCA1Frequencies.df)
    str(BRCA2Frequencies.df)
    str(MLH1Frequencies.df)
    str(USDemographics.df)
    str(ChinaDemographics.df)

    #Make a tree with no affection status, g=4 generations above, gdown=2 generations below, seed.age=50, and demographics.df=NULL which defaults to USDemographics.df.
    tree1=MakeTree()

    #Make a tree using Chinese demographics instead.
    tree2=MakeTree(demographics.df=ChinaDemographics.df)

    #Add affection statust to tree2 using BRCA1Frequencies.df which gives the BRCA1 penetrance function
    tree2a=AddAffectedToTree(tree.f=tree2,frequencies.df=BRCA1Frequencies.df)

    #make a tree with affection status (same as running MakeTree() and then AddAffectedToTree())
    tree3=MakeAffectedTrees(n=1,g=2,gdown=2,frequencies.df=MLH1Frequencies.df)
    tree4=MakeAffectedTrees(n=1,g=2,gdown=2,frequencies.df=BRCA2Frequencies.df)

    #Plots the pedigree with affection status, proband status, and genotype.  Tree is labeled with either the current age of the inidividual, age at death, or age of first affected.
    PlotPedigree(tree3)
  }

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
