\name{PlotPedigree}
\alias{PlotPedigree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots a tree/pedigree
}
\description{
This function uses the kinship2 package for easily plotting a tree.
}
\usage{
PlotPedigree(ped, affected.vector=NULL, legend.location="topleft", legend.radius=0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ped}{
    A tree generated from \code{\link{MakeTree}} or \code{\link{MakeAffectedTrees}}.
    }
  \item{affected.vector}{
    A vector indicating which individuals are afffected(1).  Note that this does not differentiate unaffected(0) and unknown affection status(2).
  }
  \item{legend.location}{
    A string indicating the placement of the legend.  These can be "topleft", "topright", "bottomleft", "bottomright".
  }
  \item{legend.radius}{
    A real number indicating the size of the legend.
  }
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
  No values returned but a plot is displayed.
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

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
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
