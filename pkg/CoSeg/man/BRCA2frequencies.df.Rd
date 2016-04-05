\name{BRCA2frequencies.df}
\alias{BRCA2frequencies.df}
\docType{data}
\title{
BRCA2 cancer incidence data frame
}
\description{
A dataframe containing annual incidence of cancers affected by the BRCA2 gene, namely ovarian cancer and breast cancer.
}
\usage{data("BRCA2frequencies.df")}
\format{
  A data frame with 800 observations on the following 5 variables.
  \describe{
    \item{\code{age}}{A numeric vector giving the age the data is relevant for.}
    \item{\code{cancer.type}}{A factor with levels \code{bBRCA2} and \code{oBRCA2} which stands for breast BRCA2 cancer and ovarian BRCA2 cancer.}
    \item{\code{female}}{A boolean/numeric vector indicating whether the data is for females(1) or males(0).}
    \item{\code{carrier}}{A boolean/numeric vector depicting whether the data is for a carrier(1) or non-carrier(0).}
    \item{\code{frequencies}}{A numeric vector giving the cancer incidence based on age, sex, carrier status, and gender.}
  }
}
\details{
Data was converted to annual incidence using a spline with 3 df's.  Risk of cancer for people below age 20 was set at 0.
}
\source{
Penetrance data from PMID: 20204502.
}
\references{
%%  ~~ possibly secondary sources and usages ~~
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
\keyword{datasets}
