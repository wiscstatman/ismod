\name{Gcandi}
\alias{Gcandi}
\title{Generate candidate sets}
\description{
  Make a candidat set L= L1 or L2
}
\usage{
Gcandi(data,pval1,pval2,error)
}
\arguments{
  \item{data}{cCGH data profile}
  \item{pval1}{probability of upper tail in binomial distribution,
    ex:.001 or .005  (to select L1) }
  \item{pval2}{percentage of tumors considered as a candidate set
  ex: .1 or .2}
  \item{error}{error term for positive internal & negative external}
}
\examples{
data(rcc)

rcc2<-rcc$data

pval1 <- 10^(-15); pval2<- .1; error <- .005

candid <- Gcandi(rcc2,pval1,pval2,error)
}
\value{
  matrix (number of candidate x number of aberration) A. A[i,j]==1
  indicates ith candidate set containing jth aberration. 
}
