\name{hammingD}
\alias{hammingD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate the hamming distance between data points.  
%%  ~~function to do ... ~~
}
\description{Hamming distance is defined on categorical vectors.  It counts the number of times the cordinates in two data vectors differ. 
}
\usage{
hammingD(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ data should be a matrix or data frame.}
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
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%See Also as \link{alphadata},
}
\examples{
### The running is time consuming
### Run hamming distance
#dis0<-hammingD(alphadata)
### Save as distance format
#REDIST<-as.dist(dis0)
### Run a hierarchical clustering using average linkage
#hc0 <- hclust(REDIST,method = "average")
### plot the dendrogram
#plot(hc0,label=xlab1,hang =-1)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
