
# Following is an example of how the different functions work.
## copied from example, and bringing in data from Rich's exome
## study/27-dat.RData



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

## data matrix yy should have rows for polyps and columns for aberrations
load("27-dat.RData")

###ysub <- t( zz2[,p2$NEW.GROUP.FULL=="Large_And_Growing"] )
ysub <- t( zz2[,p2$NEW.GROUP.FULL=="Small_And_NotGrowing"] )
tmp  <- colSums(ysub) > 0
ysub <- ysub[,tmp]  ## 45 polyps x 181 genes

# Run the Metroplis Hastings algorithm for a few scans 

# I'll use a big nskip here to favor 
mcmc <- list(nsave=1000, nskip=2000, eps=c(0.01,0.01), pup=0.10,
                nperm=5, nmove=4, padI=(3/100), padII=(1/4) )

mh <- mhisn(dat=ysub,mcmc=mcmc,tau=1,theta0=c(.2,.05) )

# One posterior summary includes marginal posterior probabilities
# that pairs of aberrations are on the network and on the same path
# The diagonal of this matrix approximates the posterior relevance probability.

# MAKE SURE THAT THE MCMC SAMPLER IS MIXING WELL BEFORE BELIEVING THESE
# APPROXIMATIONS!

ppairs <- ( postpairs( mh$partitions[-(1:10),] ) )$smat  ## a bit of burn in


save( mh, ppairs, file="02-LG.RData" )

