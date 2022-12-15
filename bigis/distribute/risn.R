"risn" <-
function(theta,partition,null=0,ntumor)
  {
   #########################################################
   #  Sample random binary vectors from the
   #  Instability-Selection-Network Model  (Newton, 2001)
   #
   # theta       (alpha, beta)
   #              alpha     rate of observable aberations
   #              beta      rate of unobservable aberations
   # partition   an nswitch-length vector with different labels
   #             for the different clusters
   # null        label of the null cluster
   # ntumor      the number of tumors to be simulated
   #
   # Method: Rejection sampling.
   #
   #########################################################

   alpha <- theta[1]; beta <- theta[2]; gamma <- 0
   nswitch <- length(partition)

   # Incidence matrix connecting switches and paths
   incidence <- outer( partition, unique(partition), "==" )
   okswitches <- !( unique(partition) == null )  # non-null switches

   xsave <- matrix(NA,ntumor,nswitch)
   for( i in 1:ntumor )
    {
     notdone <- T
     while( notdone )
      {
       # sample from prior
       xi <- ifelse( runif(nswitch) <= alpha, 1, 0 )
       # evaluate the envelope function
       stats <- (1-xi) %*% incidence
       stats <- stats[ okswitches ] 
       env <- 1 - (1-gamma)*prod(1 - beta^stats)
       # rejection step
       notdone <- ifelse( runif(1) <= env, F, T )
      }
     xsave[i,] <- xi
    }
   return(xsave)
  }
