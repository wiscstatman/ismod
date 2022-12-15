###--- >>> `dpolya1' <<<----- ~~function to do ... ~~

	## alias	 help(dpolya1)

##___ Examples ___:

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--         or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(partition,null,tau)
 {
  #####################################################################
  #
  # Evaluate the log prior mass of a Double Polya Prior-augmented state.
  #
  # This is the prior for an augmented-state partition in ISN model
  # (Used in the data-augmentation part of the MCMC algorithm.)
  #
  #
  # partition = nswitch vector holding distinct path labels
  # null = logical vector indicating active (i.e. relevant) switches
  # tau = clustering parameter: bigger <-> more paths
  #
  #####################################################################

  nswitch <- length(partition)

  # The partition componenent
  incidence <- outer( partition, unique(partition), "==" )
  sizes <- t(incidence) 
  npath <- length(sizes)
  logp <-  npath*log(tau) + sum( lgamma(sizes) ) 
  logp <- logp + lgamma(tau)-lgamma(nswitch+tau)

  # The activation component
  nnull <- sum( null )
  logp <- logp+lgamma(nnull+1)+lgamma(nswitch-nnull+1)-lgamma(nswitch+2)
  return( logp )
 }

## Keywords: '~kwd1', '~kwd2'.


