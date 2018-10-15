\name{enhcHi}
\alias{enhcHi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Performs ensemble hierarchical clustering for high dimensional categorical data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
enhcHi(x, En, len)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A nxp data matrix of data frame; n is the number of observations and p is the number of features or dimensions.}
%%     ~~Describe \code{x} here~~
%}
  \item{En}{ Number of ensembles.
%%     ~~Describe \code{B} here~~
}
  \item{len}{ Range of sizes of clusterings to run and ensemble.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Amiri, S., Clarke, B. & Clarke J. (2015). Clustering categorical data via ensembling dissimilarity matrices.  arXiv:1506.07930. 
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
data("alphadata")
### The following code generates the dissimilary matrix of sequence data stored in  alphadata
### The ensembling repeated for 100 times and the clustering range is generated randomly 
### from distrete uniform from (2,10)
#ens<-enhcHi(alphadata,En=100,len=c(2,10))
### Calculate the hamming distance
#dis0<-hammingD(ens)
### Save as distance format
#REDIST<-as.dist(dis0)
#hc0 <- hclust(REDIST,method = "average")
#plot(hc0,label=xlab1,hang =-1) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
\keyword{highdimension}