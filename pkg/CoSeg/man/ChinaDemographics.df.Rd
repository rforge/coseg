\name{ChinaDemographics.df}
\alias{ChinaDemographics.df}
\docType{data}
\title{
Chinese Demographics Data Frame
}
\description{
A dataframe containing average age of marriage, death, and number offspring by year for individuals in China.
}
\usage{data("ChinaDemographics.df")}
\format{
  A data frame with 12 observations on the following 7 variables.
  \describe{
    \item{\code{in.year}}{A numeric vector giving the input year.}
    \item{\code{out.year}}{A numeric vector giving the output year.}
    \item{\code{female.age.marriage}}{A numeric vector giving the average age of marriage for females by year.}
    \item{\code{male.age.marriage}}{A numeric vector giving the average age of marriage for males by year.}
    \item{\code{female.age.death}}{A numeric vector giving the average age of death for females by year.}
    \item{\code{male.age.death}}{A numeric vector giving the average age of death for males by year.}
    \item{\code{offspring}}{A numeric vector giving the average numbers of offspring for families by year.}
  }
}
%%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%%}
\source{
Data were taken from the following sources.
\url{http://www.digitalhistory.uh.edu/}
\url{http://demog.berkeley.edu/~andrew/1918/figure2.html}
\url{http://www.infoplease.com/ipa/A0005148.html}
\url{http://www.infoplease.com/ipa/A0005061.html}
\url{http://www.gapminder.org}
}
%%\references{
%%  ~~ possibly secondary sources and usages ~~
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
\keyword{datasets}
