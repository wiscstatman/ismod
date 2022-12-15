\name{ISGN}
\alias{ISGN}
\title{main MCMC function}
\description{
  fit the ISGN model using 3 MCMC move type
}
\usage{
ISGN(data, mcmc, para, pri, candid)
}
\arguments{
  \item{data}{cCGH data profile}
  \item{mcmc}{list with driving instructions
    
    padI   = roportion of switches to change in last move type
    
    nskip  = size of stretch between saved states
    
    nsave  = number of saved states
    
    tmax   = number of maximum candidble in a network
  }
  \item{para}{(alpha, p.a, p.b), here

    alpha = Prob[ X=1 ] : estimated from not significant aberrations

    p.a = Prob[ Aberration | X=1 ] : error terms; .99 or .98...

    p.b = Prob[ Aberration | X=0 ] : error terms; .1 or .2...
  }
  \item{pri}{prior, uniform or cadidional uniform}
  \item{candid}{candidate sets from Gcandi}
}
\value{return 'nsave' general networks with their likelihood and
  posterior probability, and best network.
}
