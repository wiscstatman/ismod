"rpolya" <-
function(n,tau)
	{
        #############################################################
        #
        # Simulate a Standard Uniform Polya partition, slowly
        #   (used to start the MCMC sampler)
        #
	# n = number of objects to be partitioned
	# tau > 0 = hyperparameter
        #
        #############################################################
	v <- runif(n)
	x <- v
	if( n >= 2 )
 	 {
  	  for( i in 2:n )
	   {
	    pwild <- tau/( tau + i - 1 )
	    if(  runif(1) > pwild ) # select an old one
	     {
	       ind <- sample( 1:(i-1), size=1 )
	       x[i] <- x[ind]
	     }
	   }
	 }
	return(x)
	}
