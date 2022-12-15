\name{Glik}
\alias{Glik}
\title{Calculate the log joint probability density}
\description{
  Calculate the log joint probability density of aberrations in the ISGN
  model. Implementation was suggested by Jay Kadane.
}
\usage{
Glik(C,data,params)
}
\arguments{
  \item{C}{(General network): kxn matrix; k ensembles, n aberrations.
    C[i,j]=1 if jth aberration resides on ith ensemble}
  \item{data}{cCGH data profile}
  \item{params}{(alpha, p.a, p.b), here
    
      alpha = Prob[ X=1 ] : estimated from not significant aberrations

      p.a = Prob[ Aberration | X=1 ] : error terms; .99 or .98...

      p.b = Prob[ Aberration | X=0 ] : error terms; .1 or .2...
    }
  }
\value{
    log joint probability density 
}
