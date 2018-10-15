\name{mixmodCluster}
\alias{mixmodCluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Model-based clustering of categorical data
%%  ~~function to do ... ~~
}
\description{
This function runs model-based clustering as implemented in library(Rmixmod).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
 mixmodCluster(data, nbCluster,model,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ frame containing quantitative, qualitative or composite data. Rows correspond to observations and columns correspond to variables. Data values should be characters.
%%     ~~Describe \code{data} here~~
}
  \item{nbCluster}{ number of clusters
%%     ~~Describe \code{k} here~~
}
  \item{model}{ Model of data, for the categorical data mixmodMultinomialModel().
%%     ~~Describe \code{k2} here~~
}
}
\details{
See library(Rmixmod) for details.  
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(Rmixmod)
data(zoo)
### the following commands change the data to character values
zoo2<-matrix(,nrow=nrow(zoo$obs),ncol=ncol(zoo$obs))
for(ii in 1:ncol(zoo$obs)){
  zoo2[,ii]<-as.character(zoo$obs[,ii])
}
### clustering data into seven clusters
xem2<-mixmodCluster(data.frame(zoo2), 7,model = mixmodMultinomialModel())
### output of clustering can be obtained by the following 
xem2[12]
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
\keyword{lowdimension}