\name{FormatWebToCoSeg}
\alias{FormatWebToCoSeg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots a tree/pedigree
}
\description{
this function takes a pedigree in the format for analyze.myvariant.org which has ordered columns with names that may not be meaningful to one with meaningful names for use in CoSeg
}
\usage{
FormatWebToCoSeg(ped)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ped}{
    A text file that follows the format on \url{www.analyze.myvariant.org}.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
  A pedigree object for use with CoSeg.
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

\dontrun{CoSegPed=FormatWebToCoSeg(WebPed)}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
