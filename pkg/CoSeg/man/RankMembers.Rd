\name{RankMembers}
\alias{RankMembers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to calculate the likelihood ratio
}
\description{
This function calculates the likelihood ratio for an allele causing a disease asssuming that the allele is extremely rare so that all family members who have the allele got it directly from a common ancestor in the pedigree.
}
\usage{
RankMembers(ped, affected.vector, gene="BRCA1")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ped}{
  A dataframe of a pedigree containing id, momid, dadid, age(Current age or age affected), y.born(year born), female(a logical where 1 is female and 0 is male), and genotype(0 is non-carrier, 1 carrier, and 2 unknown) in any order but with these specific column names.
}
  \item{affected.vector}{
  A vector giving the affection status of the individual.  0 is unaffected, 1 is affected, and 2 is unknown.
}
  \item{gene}{
  A text string (e.g. BRCA1, BRCA2, or MLH1) giving the gene name for analysis.
  }
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{likelihood.ratio }{The resultant likelihood ratio}
  \item{separating.meioses }{The number of meioses separating all individuals known to have the genotype including those that are are indirectly known to have the genotype.}
  \item{number.genotypes.found }{The number of permissible genotypes found for the pedigree.}

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
%%\keyword{}
%%\keyword{}% __ONLY ONE__ keyword per line
