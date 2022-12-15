
# Following is an example of how the different functions work.

# Necessary functions

source("risn.R")         # Function to sample a random profile, ISN model
source("disn.R")         # Loglikelihood function, ISN model
source("dpolya1.R")      # Logprior function, for augmented-state partition
source("dpolya2.R")      # Logprior function, for actual partition
source("rpolya.R")       # Function to sample a Uniform Polya sequence
source("mhisn.R")        # Function to run Metropolis-Hastings sampler
source("mixsample.R")    # random sample of T/F's for MH sampler
source("postpairs.R")    # Marginal posterior probability summaries 

# Simulate a small data set with 12 potential aberrations (switches)
# arranged on three oncogenic pathways, and with 6 inactive switches,
# and with 100 tumors.

parta <- c( rep(.1,2), .2, rep(.3,3), rep(0,6)  )
yy <- risn( c(.2,.05), partition=parta,  ntumor=100 )

# Evaluate the log probability of these data under, say, the true 
# parameter values, and another value.

loglika <- disn( c(.2,.05), data=yy, partition=parta, null=(parta==0) )

partb <- parta
partb[3] <- .1      # This partition has two paths instead of three
loglikb <- disn( c(.2,.05), data=yy, partition=partb, null=(partb==0) )

# Evaluate the log prior probability of the two partitions.

lpriora <- dpolya2( parta, null=(parta==0), tau=5 )
lpriorb <- dpolya2( partb, null=(partb==0), tau=5 )

# Run the Metroplis Hastings algorithm for a few scans 
# (With on the order of 100 tumors and 50-100 aberrations,
#  I have used nsave=2000, nskip=1000; this can take a lot of CPUs!)

mcmc <- list(nsave=1000, nskip=500, eps=c(0.01,0.01), pup=0.10,
                nperm=5, nmove=4, padI=(3/100), padII=(1/4) )

mh <- mhisn(dat=yy,mcmc=mcmc,tau=1,theta0=c(.2,.05) )

# One posterior summary includes marginal posterior probabilities
# that pairs of aberrations are on the network and on the same path
# The diagonal of this matrix approximates the posterior relevance probability.

# MAKE SURE THAT THE MCMC SAMPLER IS MIXING WELL BEFORE BELIEVING THESE
# APPROXIMATIONS!

ppairs <- ( postpairs( mh$partitions ) )$smat
