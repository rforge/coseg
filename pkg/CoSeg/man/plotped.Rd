\name{plotped}
\alias{plotped}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots a tree/pedigree
}
\description{
This function uses the kinship2 package for easily plotting a tree.
}
\usage{
plotped(ped)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ped}{
    A tree generated from \code{\link{make.tree}} or \code{\link{make.affected.trees}}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  No values returned but a plot is displayed.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
  plotped(tree3)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
