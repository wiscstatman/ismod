"postpairs" <-
function( partitions, labels=NULL )
	{
        # Function to process a sample of partitions
        # and record three symmetric matrices:
        #
        #   smat[i,j] = P[ i and j in network and on same path | data  ]
        #   amat[i,j] = P[ i and j in network | data ]
        #   dmat[i,j] = P[ i and j in network and on diff. paths | data ]
 	#
	# Partitions is a matrix each of whose rows denotes a 
	# partition.  Zero entries are inactive; common values on same path.
	#
	nswitch <- ncol(partitions)
	nsave <- nrow(partitions)
	smat <- matrix(0,nswitch,nswitch)
	amat <- smat
	dimnames( smat ) <- list( labels, labels )
	dimnames( amat ) <- list( labels, labels )
	for( i in 1:nsave )
 	{
	  active <- (partitions[i,] > 0)
	  aa <- active %o% active    # active?
	  cl <- outer( partitions[i,], partitions[i,], "==" )  # same path?
	  smat <- smat + aa*cl
	  amat <- amat + aa    # just askes if active
	 }
	dmat <- amat-smat  # active and on different paths
	amat <- amat/nsave
	smat <- smat/nsave
	dmat <- dmat/nsave
	return( list( smat=smat, amat=amat, dmat=dmat ) )
	}
