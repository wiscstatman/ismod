\name{ISTN}
\alias{ISTN}
\title{main MCMC function}
\description{
  fit the ISTN model using 8 MCMC move types
}
\usage{
ISTN(data,mcmc,tau,para,res)
}
\arguments{
  \item{data}{cCGH data profile}
  \item{mcmc}{list with driving instructions

    nskip  = size of stretch between saved states

    nsave  = number of saved states

    nperm  = number of switches permuted in 2nd movetype update

    padI =  proportion of switches to change in first neutral vector
    update (3rd movetype)

    padII (<1/2) = chance of full activation / deactivation in second
    neutral vector update (4th movetype)
  }
  \item{tau}{Positive number, cluster prior hyperparameter}
  \item{para}{(alpha, p.a, p.b), here

    alpha = Prob[ X=1 ] : estimated from not significant aberrations

    p.a = Prob[ Aberration | X=1 ] : error terms; .99 or .98...

    p.b = Prob[ Aberration | X=0 ] : error terms; .1 or .2...
  }
  \item{res}{Probability matrix used for making a sub-tree. Refer
    function res. Recommended to save res<-res(20).} 
}
\examples{
  res<-dget('res') # it is highly recommended to save res<-res(20)

  mcmc <- list(nsave=10, nskip=100, nperm=3, padI=.15, padII=.15)

  data(rcc)

  data<-rcc$data

  para<-c(0.05, 0.01, 0.01);tau<-4

  mh <- ISTN(data,mcmc,tau,para,res)
}
\value{return 'nsave' tree-like networks with their likelihood and
  posterior probability, and best network.
}
