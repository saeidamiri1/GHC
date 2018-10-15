\name{enhc}
\alias{enhc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Performs ensemble hierarchical clustering for categorical data.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
enhc(x, En)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An nxp data matrix; n is the number of observations and p is the number of dimensions.}
%%     ~~Describe \code{x} here~~
%}
  \item{En}{ Number of ensembles.
%%     ~~Describe \code{B} here~~
}
}
\details{
The function 'enhc' generates a dissimilarity via the ensemble, the dissimilarity can be used to obtain the clustering of data, see Amiri et al. (2015) for details.
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
\author{ Saeid Amiri, Bertrand Clarke and Jennifer Clarke.
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
### zoo includes the zoo data downloaded from UCI
### MAchine Learning Repository
data('zoo')
### Calculate ensemble dissimilarities with 150 ensembles
disten<-enhc(zoo$obs,En=150)
###  hclust is R function to perform a hierarchical cluster analysis using 
### dissimilarities obtained by ensembling
en<-hclust(disten,method='average')
### Dendrogram can be obtained using
plot(en,label=zoo$lab,hang=-1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
\keyword{lowdimension}