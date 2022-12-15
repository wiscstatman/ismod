"mixsample" <-
function(n)
	{
	# n >= 2
	# Returns a random vector of T's and F's uniformly among those
	# that are not all T or all F
	# Method: rejection sampling
	notdone <- T
	while( notdone )
	 {
	  x <- sample( c(T,F), size=n, replace=T )
	  notdone <- (sum(x)==0) | (sum(x) == n)
	 }
	return(x)
	}
