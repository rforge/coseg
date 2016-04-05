\name{CalculateLikelihoodRatio}
\alias{CalculateLikelihoodRatio}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to calculate the likelihood ratio
}
\description{
This function calculates the likelihood ratio for an allele causing a disease asssuming that the allele is extremely rare so that all family members who have the allele got it directly from a common ancestor in the pedigree.
}
\usage{
CalculateLikelihoodRatio(ped, affected.boolean)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ped}{
A dataframe of a pedigree containing id, momid, dadid, age(Current age or age affected), y.born(year born), female(a logical where 1 is female and 0 is male), and genotype(0 is non-carrier, 1 carrier, and 2 unknown) in any order but with these specific column names.
}
  \item{affected.boolean}{
A boolean(or logical) vector where 0 is unaffected and 1 is affected.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{likelihood.ratio }{The resultant likelihood ratio}
  \item{separating.meioses }{The number of meioses separating all individuals known to have the genotype including those that are are indirectly known to have the genotype.}
  \item{number.genotypes.found }{The number of permissible genotypes found for the pedigree.}

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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{}
\keyword{}% __ONLY ONE__ keyword per line
