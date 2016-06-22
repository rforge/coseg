\name{add.affected.to.tree}
\alias{add.affected.to.tree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add affection status to a given tree
}
\description{
A function that takes a tree and adds affection status using user given cancer types and penetrances.
}
\usage{
add.affected.to.tree(tree.f, frequencies.df = NULL, g = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree.f}{
	Output from \code{\link{make.tree}} function.
}
  \item{frequencies.df}{
	A dataframe giving the penetrance of the desired affection.  For convenience this can be one of the included penetrance data frames (\code{\link{BRCA1frequencies.df}}, \code{\link{BRCA2frequencies.df}}, \code{\link{LSfrequencies.df}}), or a user defined penetrance similar to one of those.
}
  \item{g}{
	A parameter giving the generation to select a proband in.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Adds the affection status the inputted tree.  The added columns are based on the inputted frequencies.df.  For the default BRCA1frequencies.df the extra columns are listed below.
  \item{famid}{An integer id for the current tree.  This is only useful for multiple trees.}
  \item{bBRCA1.d}{A boolean/integer depicting whether the individual has breast cancer(1) or not(0).}
  \item{oBRCA1.d}{A boolean/integer depicting whether the individual has ovarian cancer(1) or not(0)}
  \item{bBRCA1.aoo}{An integer giving the age of onset of the breast cancer if the individual is affected or NA if not.}
  \item{oBRCA1.aoo}{An integer giving the age of onset of the ovarian cancer if the individual is affected or NA if not.}
  \item{proband}{A boolean/integer labeling the proband(one with the variant and the disease) with a 1 and everyone else with 0.  Note that if no suitable proband is found, all members will receive a -1.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
John Michael O. Ranola and Brian H. Shirts
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{make.tree}}, \code{\link{make.affected.trees}}, and \code{\link{plot.ped}}
%\code{\link{BRCA1frequencies.df}}, \code{\link{BRCA2frequencies.df}}, and \code{\link{LSfrequencies.df}}
}
\examples{

#Load all the data included in the CoSeg package.
data(BRCA1frequencies.df, package="CoSeg")
data(BRCA2frequencies.df, package="CoSeg")
data(LSfrequencies.df, package="CoSeg")
data(USDemographics.df, package="CoSeg")
data(ChinaDemographics.df, package="CoSeg")

#summaries of all the data
str(BRCA1frequencies.df)
str(BRCA2frequencies.df)
str(LSfrequencies.df)
str(USDemographics.df)
str(ChinaDemographics.df)

#Make a tree with no affection status, g=4 generations above, gdown=2 generations below, seed.age=50, and demographics.df=NULL which defaults to USDemographics.df.
tree1=make.tree()

#Make a tree using Chinese demographics instead.
tree2=make.tree(demographics.df=ChinaDemographics.df)

#Add affection statust to tree2 using BRCA1frequencies.df which gives the BRCA1 penetrance function
tree2a=add.affected.to.tree(tree.f=tree2,frequencies.df=BRCA1frequencies.df)

#make a tree with affection status (same as running make.tree() and then add.affected.to.tree())
tree3=make.affected.trees(n=1,g=2,gdown=2,frequencies.df=LSfrequencies.df)
tree4=make.affected.trees(n=1,g=2,gdown=2,frequencies.df=BRCA2frequencies.df)

#Plots the pedigree with affection status, proband status, and genotype.  Tree is labeled with either the current age of the inidividual, age at death, or age of first affected.
plot.ped(tree3)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
