\name{MLH1Frequencies.df}
\alias{MLH1Frequencies.df}
\docType{data}
\title{
Lynch Syndrome cancer incidence data frame
}
\description{
A dataframe containing annual incidence of cancers affected by the Lynce Syndrome gene, namely colorectal cancer, endometrial cancer, and other minor cancer.
}
\usage{data("MLH1Frequencies.df")}
\format{
  A data frame with 1200 observations on the following 5 variables.
  \describe{
    \item{\code{age}}{A numeric vector giving the age the data is relevant for.}
    \item{\code{cancer.type}}{A factor with levels \code{bBRCA1} and \code{oBRCA1} which stands for breast BRCA1 cancer and ovarian BRCA1 cancer.}
    \item{\code{female}}{A boolean/numeric vector indicating whether the data is for females(1) or males(0).}
    \item{\code{carrier}}{A boolean/numeric vector depicting whether the data is for a carrier(1) or non-carrier(0).}
    \item{\code{frequencies}}{A numeric vector giving the cancer incidence based on age, sex, carrier status, and gender.}
  }
}
\details{
Data was converted to annual incidence using a spline with 3 df's.  Risk of cancer for people below age 20 was set at 0.
}
\source{
Penetrance data from PMID: 15937084 and	PMID: 19723749.	SEER from \url{http://seer.cancer.gov/}
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
