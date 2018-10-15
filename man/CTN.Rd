\name{CTN}
\alias{CTN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{convert genetic data (nucleotides) to numerical values
%%  ~~function to do ... ~~
}
\description{This function converts genetic data (nucleotides) to numeric data. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CTN(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{data  should be in fasta format
%%     ~~Describe \code{data} here~~
}
}

\details{
R is faster on the numerical data and storage of data via numerical values takes less capacity. Genetic data includes nucleotide data {A,T,C,G} and are usually saved in Fasta format. After downloading the data from one of the bioinformatics Repository, e.g., \url{http://www.viprbrc.org/}, and importing it to R, this function converts the data to numerical values. 
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
### load libarary 'seqnir' that includes several functions to handle genetics data 
library('seqinr')
### import fasta data to R.
#x.dna0 <- read.fasta("dna.fasta")
###  convert data to numerical values
#x.dna<-CTN(x.dna0)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
