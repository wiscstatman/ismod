"disn" <-
function(theta,data,partition,null)
{
 #######################################################################
 #
 # Calculate the log joint probability density of 
 # aberrations in the Instability-Selection-Network model
 #
 # theta = (alpha, beta)  in (0,1)^2
 # data = ntumor x nswitch matrix of binaries 0/1 indicating abnormality
 # partition = nswitch vector holding distinct path labels
 # null = logical vector indicating null switches
 #        (note this is slightly different than the null label in risn.R)
 #
 #######################################################################

  alpha <- theta[1]; beta <- theta[2]
  gamma <- 0

  # First get the prior part of the log density
  s0 <- sum( data==0 );  s1 <- sum( data==1 )
  logprob <- s0*log( 1-alpha ) + s1*log( alpha ) 
  
  # In case all switches are in the null path we're done; otherwise
  if( !( all( null ) ) )
   {
    #  Now compute the remaining sufficient statistics. 
    #  i.e., for each path j, in each tumor i, we need S_ij, the
    #  number of non-defects.
    #  To do so create an incidence matrix nswitch x npath 
    #  indicating the path label of each switch.
    apart <- partition[!null]    # active partition
    adat <-  data[,!null]        # active data
    incidence <- outer( apart, unique(apart), "==" )
    npath <- ncol(incidence)
    stats <- (1-adat) %*% incidence  # this is ntumor x npath
    # Now compute the numerator of the log-density
    one <- matrix(1,npath,1)
    tmp <- c( log( 1-beta^stats ) %*% one )
    numer <- 1 - (1-gamma)*exp(tmp)
    # Compute the denominator (same for all tumors)
    # First get the path lengths
    plen <- rep(1,sum(!null)) %*% incidence
    denom <- 1-(1-gamma)*prod( 1-(1-(1-alpha)*(1-beta) )^plen  )
    # Put it all together
    logprob <- logprob + sum( log(numer) - log( denom ) )
   }
  return(  logprob )
 }
