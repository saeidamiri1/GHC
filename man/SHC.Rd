\name{SHC}
\alias{SHC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Performs  Stabilized Hybrid Clustering for continuous data.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SHC(x,K,B=200,knmin,knmax)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An nxp data matrix; n is the number of observations and p is the number of dimensions.}
  \item{K}{Size of cluster}
%%     ~~Describe \code{x} here~~
%}
  \item{B}{Number of iterations.
%%     ~~Describe \code{B} here~~
}
 \item{knmin}{Minimum size of cluster to get the stabilized clustering, default knmin=2}
 \item{knmax}{Maximum size of cluster to get the stabilized clustering, default knmax=n/5}
}
\details{
The function 'SHC' generates the Stabilized Hybrid Clustering (SHC) presented in Algorithm 1, see Amiri et al. (2018) for details.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Amiri, S., Clarke, B, & Clarke, J. (2018). A General Hybrid Clustering Technique. To appear in Journal of Computational and Graphical Statistics.
}
\author{Saeid Amiri, Bertrand Clarke and Jennifer Clarke.
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
knmin0<-2
knmax0<-floor(dim(spiral)[1]/5)
CLUS<-SHC(as.matrix(spiral),3,B=20,knmin=knmin0,knmax=knmax0)
# plot the dendrogram
plot(hclust(CLUS[[1]],method="single"),h=-1)
# print the assigned clusters to observation
 print(CLUS[[2]])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
\keyword{clustering}
