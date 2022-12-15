###--- >>> `dpolya2' <<<----- ~~function to do ... ~~

	## alias	 help(dpolya2)

##___ Examples ___:

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--         or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(partition,null,tau)
 {
  ####################################################################
  #
  # Evaluate the log of the probability mass of a double-Polya prior.
  # This is used as the prior distribution for the network structure
  # in the  Instability-Selection-Network Model.
  #
  # partition = nswitch vector holding distinct path labels
  # null = logical vector indicating active (i.e. relevant) switches
  # tau = clustering parameter: bigger <-> more paths
  #
  ####################################################################

  # first the null part of the logprior  (beta-binomial piece)
  nswitch <- length(partition)
  nnull <- sum( null )
  logp <- lgamma(nnull+1)+lgamma(nswitch-nnull+1)-lgamma(nswitch+2)

  # Now the cluster prior on the other switches
  if( !all(null) )
   {
    mm <- nswitch-nnull
    apart <- partition[!null] 
    incidence <- outer( apart, unique(apart), "==" )
    sizes <- t(incidence) 
    npath <- length(sizes)
    logp <- logp + npath*log(tau) + sum( lgamma(sizes) ) 
    logp <- logp + lgamma(tau)-lgamma(mm+tau)
   }
  return( logp )
 }

## Keywords: '~kwd1', '~kwd2'.


