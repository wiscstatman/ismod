\name{Trisn}
\alias{Trisn}
\title{generate random binary vectors}
\description{
  Sample random binary vectors from the tree-like network using
  rejection sampling method
}
\usage{
Trisn(rates,partition,parent,neutral,numtumor)
}
\arguments{
  \item{rates}{(alpha, p.a, p.b), here
    
      alpha = Prob[ X=1 ] : estimated from not significant aberrations

      p.a = Prob[ Aberration | X=1 ] : error terms; .99 or .98...

      p.b = Prob[ Aberration | X=0 ] : error terms; .1 or .2...}
    \item{partition}{n component vector holding distinct edge labels
      0 is coded for neutral aberration.}
    \item{parent}{like partition, but holds labels for parent edges.
       0 is coded for root node}
    \item{neutral}{logical vector indicating neutral components}
    \item{numtumor}{number of tumors to be sampled}
}
\examples{
  para <- c(0.05, 0.01, 0.01)

  partition<-c(1,2,3,4,5,6, rep(1,14)) 

  parent<-c(0,0,1,1,2,2,rep(0,14))

  neutral<-c(rep(0,6),rep(1,14)) 

  numtumor<-100

  data1 <- Trisn(para,partition, parent,neutral, numtumor)
}
\value{
  binary vectors having 1 (indicating aberration) and 0.
}
