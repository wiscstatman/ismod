\name{Grisn}
\alias{Grisn}
\title{generate random binary vectors}
\description{
  Sample random binary vectors from general netowrk, C,
  using rejection sampling method
}
\usage{
Grisn(params,C,numtumor)
}
\arguments{
  \item{params}{(alpha, p.a, p.b), here
    
      alpha = Prob[ X=1 ] : estimated from not significant aberrations

      p.a = Prob[ Aberration | X=1 ] : error terms; .99 or .98...

      p.b = Prob[ Aberration | X=0 ] : error terms; .1 or .2...}
  \item{C}{General network: kxn matrix; k ensembles, n aberrations
       C[i,j]=1 if jth aberration resides on ith ensemble}
  \item{numtumor}{number of tumors to be sampled}
}
\value{
  binary vectors having 1 (indicating aberration) and 0.
}
