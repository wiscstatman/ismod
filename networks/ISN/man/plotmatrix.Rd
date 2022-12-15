\name{plot.matrix}
\alias{plot.matrix}
\title{plot cCGH profile}
\description{plot cCGH profile}
\usage{plot.matrix(data}}
\arguments{
  \item{data}{cCGH data}
}
\value{
  plot the number of interactions
}
\author{David Dahl}
\examples{
  data(rcc) # renal cell carcinoma data]

  rcc2<-rcc$data

  plot.matrix(rcc2,xlab="Aberrations", ylab="RCC Tumor ID") 
}
\keyword{plots}
