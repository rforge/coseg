\name{AddAffectedToTrees}
\alias{AddAffectedToTrees}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add affection status to a given tree
}
\description{
A function that takes a tree and adds affection status using user given cancer types and penetrances.
}
\usage{
AddAffectedToTrees(tree.f, frequencies.df = NULL, g = 4, benign.boolean=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree.f}{
	Output from \code{\link{MakeTree}} function.
}
  \item{frequencies.df}{
	A dataframe giving the penetrance of the desired affection.  For convenience this can be one of the included penetrance data frames (\code{\link{BRCA1Frequencies.df}}, \code{\link{BRCA2Frequencies.df}}, \code{\link{MLH1Frequencies.df}}), or a user defined penetrance similar to one of those.
}
  \item{g}{
	A parameter giving the generation to select a proband in.
}
  \item{benign.boolean}{
  A boolean variable which tells the program to simulate disease status for a benign variant when benign.boolean is TRUE.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
Adds the affection status the inputted tree.  The added columns are based on the inputted frequencies.df.  For the default BRCA1Frequencies.df the extra columns are listed below.
  \item{famid}{An integer id for the current tree.  This is only useful for multiple trees.}
  \item{bBRCA1.d}{A boolean/integer depicting whether the individual has breast cancer(1) or not(0).}
  \item{oBRCA1.d}{A boolean/integer depicting whether the individual has ovarian cancer(1) or not(0)}
  \item{bBRCA1.aoo}{An integer giving the age of onset of the breast cancer if the individual is affected or NA if not.}
  \item{oBRCA1.aoo}{An integer giving the age of onset of the ovarian cancer if the individual is affected or NA if not.}
  \item{proband}{A boolean/integer labeling the proband(one with the variant and the disease) with a 1 and everyone else with 0.  Note that if no suitable proband is found, all members will receive a -1.}
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
See also \code{\link{MakeTree}}, \code{\link{MakeAffectedTrees}}, and \code{\link{PlotPedigree}}
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
