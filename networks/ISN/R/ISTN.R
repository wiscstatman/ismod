plot.matrix <- function(A, rsize=0.25, main=NULL, xlab=NULL, ylab=NULL,
                        axes=TRUE, col=c("wheat","black","red","blue"), ...) {
#  plot.matrix, an R function to display matrices contains zeros or ones.
#  Written by David Dahl, 8 February 2000
    n <- nrow(A)
    m <- ncol(A)

    max.nm <- max(n,m)
    min.dim <- max.nm*(1+rsize)

    plot.new()
    min.pin <- min(par()$pin)

    new <- par()$pin/min.pin * min.dim
    par(usr=c(0,new[1],0,new[2]))

    data <- matrix(NA,nrow=n*m,ncol=3)
    rect(0,0,m,n,col=col[1],border=col[2])
    for(i in 1:n) {
      for(j in 1:m) {
        if(A[i,j])
          rect((j-1),(i-1),j,i, col=col[2], border=col[1], lty="blank", lwd=1)
      }
    }
    rect(0,0,m,n,border=col[2])
    buffer <- max.nm/50
    abs.bp.height <- (min.dim - max.nm) - buffer

    rtotals <- (A %*% matrix(1,nrow=m))[,]
    rmultiplier <- abs.bp.height / max(rtotals)
    for(i in 1:n) {
      rect((m+buffer),(i-1),(m+buffer+rmultiplier*rtotals[i]),i,col=col[3])
    }

    ctotals <- (matrix(1,ncol=n) %*% A)[,]
    cmultiplier <- abs.bp.height / max(ctotals)
    for(j in 1:m) {
      rect(j-1,n+buffer,j,n+buffer+cmultiplier*ctotals[j],col=col[4])
    }

    if (!is.null(main)) {
      mtext(main, side=3, at=m/2, line=1, cex=1.5)
    }
    if (!is.null(xlab)) {
      mtext(xlab, side=1, at=m/2, line=3)
    }
    if (!is.null(ylab)) {
      mtext(ylab, side=2, at=n/2, line=3)
    }
    if (axes) {
        axis(1, pos=0, at=seq(0,m,by=10), ...)
        axis(2, pos=0, at=seq(0,n,by=10), las=2, ...)
    }    
    invisible()
}

"res"<-function(nswitch){
  #########################################################################  
  # Probability matrix used for making a sub-tree.
  # For M edges, res[M,n] indicates the probability that n (3<=n<=M) is being
  # together to make a sub-tree. It is highly recommended to save 20 by 20
  # 'res' matrix using this function and use it.
  # nswitch: number of edges
  #########################################################################  
  res<-matrix(0, nswitch,nswitch)
  diag(res)<-1
  for(i in 4:nswitch){
    gg<-p1(i, null)$set
    res[i,c(3:(i-1))]<-gg/(gg[i-3]+1)
  }
  return(res)
}

"subp" <- function(g,j){
  ########################################################################
  # Used inside function "pri" and p1
  k1<-min(g,j)
  tmp<-0
  if(g>=3){
    for(i in 3:k1){
      m<-p1(i-1,set)
      tmp<-tmp+choose(g,i)*i*m$ans*subp(g-i, i)
    }
    return(1+tmp)
  }
  else
    return(1)
}

"p1"<-function(n,set){
  ########################################################################
  # Used inside function "pri"
  tmp<-0; set<-NULL
  if(n>3){
    for(i in 3:(n-1)){
      m<-p1(i-1,set)
      tmp<-tmp+choose(n,i)*i*m$ans*subp(n-i, i)
      set<-c(set, tmp)
    }
    ans<-1+tmp; 
    result <- list(ans=ans, set=set)
    return(result)
  }
  if(n==1) ans<-0
  if(n==2||n==3)  ans<-1;
  result <- list(ans=ans, set=set)
  return(result)
}

"pri"<-function(n){
  ######################################################################
  # Return the number of possible trees using n edges
  # n: number of edges. 
  ######################################################################
  if(n>2) return(p1(n, NULL)$ans+n*p1(n-1, NULL)$ans)
  else return(1)
}

"ISTN" <- function(data,mcmc,tau,para,res){
  #########################################################################  
  # data = binary matrix, ntumor x nswitch, holding CGH profiles
  # mcmc = list with driving instructions
  #        nskip  = size of stretch between saved states
  #        nsave  = number of saved states
  #        nperm  = number of switches permuted in 2nd movetype update
  #        padI =  proportion of switches to change in first neutral vector
  #                update (3rd movetype)
  #        padII (<1/2) = chance of full activation / deactivation in second
  #                       neutral vector update (4th movetype)
  # tau > 0.  The cluster prior hyperparameter
  # para      theta, gamma, delta
  # res       Probability matrix used for making a sub-tree.
  #           Refer function res.
  # 
  # ex.
  # source('tree.R')
  # res<-dget('res') # it is highly recommended to save res<-res(20)
  # mcmc <- list(nsave=10, nskip=100, nperm=3, padI=.15, padII=.15)
  # dd<-dget('trim'); data<-dd$gradeI; 
  # para<-c(0.05, 0.01, 0.01);tau<-4
  # mh <- ISTN(data,mcmc,tau,para,res)
  ##########################################################################
  nsave<-mcmc$nsave; nskip<-mcmc$nskip; nperm<-mcmc$nperm;
  padI<-mcmc$padI; padII<-mcmc$padII; 
  nswitch<-ncol(data)
  # given n edges, the number of possible trees.
  # one can get these numbers using function "pri" but it is highly recommended
  # to save the result and use it.
  # prior<-rep(0,25); for(i in 1:25){prior[i]<-pri(i)}
  prior<-c(1,1,4,17,116,997,10270,115257,1539928, 23414021,3.911842e+08, 
     7.173987e+09,1.440993e+11, 3.170871e+12, 7.486385e+13, 1.895879e+15,
     5.117341e+16, 1.475958e+18,4.497267e+19, 1.450948e+21, 4.923465e+22, 
     1.760436e+24, 6.581871e+25, 2.579816e+27,1.053847e+29)
  # given n edges, the number of possible typeII trees; a tree composed of n 
  # edges and whose root has more than two child edges. 
  # rres<-p1(n, NULL)$ans
  rres <- c(0,1,1,13,51,691,5433,71793,893791,14476111,2.319470e+08, 
     4.390623e+09,8.702118e+10, 1.952574e+12, 4.557523e+13,1.166675e+15, 
     3.133992e+16, 9.118389e+17,2.764773e+19, 8.979930e+20,3.037680e+22, 
     1.092146e+24, 4.069935e+25,1.603032e+27, 6.530885e+28)
  # a probability being a typeII tree. 
  zerop <- c(rres/prior)
  prior<- -log(prior)
  index<-c(1:nswitch)

  # Create objects to save MCMC output.
  nscan <- nsave*nskip
  pasave <- matrix(NA,nsave,nswitch);  # saved parents
  psave <- matrix(NA,nsave,nswitch);   # saved partitions
  pneu <- matrix(NA,nsave,nswitch);    # saved neutral
  liksave <- numeric(nsave);      # saved log-likelihood 
  acrate <- rep(0,8);           # acceptance rates
  
  # make an inital tree
  init <- Tnet(nswitch, tau, zerop, res)
  partition <- init$x
  parent <- init$tree
  neutral <- init$neutral
 
  # Evaluate the loglikelihood and logprior of starting state
  loglik <- Tlik(para, partition, parent, neutral,data)
  logpri <- dpolya1(partition,neutral,tau, prior)

  ####################################################################
  # Run MCMC.
  # Propose a new network configuration (8 move types)
  ####################################################################

  skipcount <- 0
  isave <- 1
  for(iscan in 1:nscan){
    skipcount <- skipcount+1; 

    #######################################################################
    # Update I: change the tree structure; choose one parent edge and among
    # its child edges, choose some and change the whole structure. 
    ########################################################################
    logmh <- NA
    cindex <- sample(parent, size=1);  
    cchild <- unique(partition[parent==cindex])
    lcchild <- length(cchild)
    lpa <- sum(is.element(parent, cindex)) 
    csize<-sample(c(1:lcchild), size=1) 
    # choose some of child edges
    ccandi<- sample(cchild, size=csize) 
    leaf<-setdiff(partition, parent)
    ok<-T; tmp<-NULL; lc<-1
    while(ok){
      tmp<-partition[is.element(parent,ccandi)]
      ccandi<-unique(c(ccandi,tmp))
      nlc<-length(ccandi)
      if((nlc-lc)>0) lc<-nlc
      else ok<-F
    }
    tindex<-index[is.element(partition, ccandi)]
    lden<-0; lnom<-0; 
    ltindex <- length(tindex)
    if(ltindex>1){
      npar <- partition; 
      rok<-T
      while(rok){
        npar[tindex]<-runif(ltindex); 
        if(sum(is.element(npar[tindex], partition))==0) rok<-F
      }
      epartition <- npar[tindex[1]]
      np<-parent; np[tindex]<-cindex
      for( i in 2:ltindex){
        pwild <- tau/(tau+i-1)
        if( runif(1) > pwild ){ # don't take the new label
          if(length(epartition)==1) npar[tindex[i]]<-epartition
          else{
            sam <- sample( epartition, size=1 ) 
            npar[tindex[i]] <- sam 
          }
        }
        epartition <- c(epartition, npar[tindex[i]] )
      }#for end
      tp<-np[tindex]
      tpar<-npar[tindex]
      ltpar <- length(unique(tpar))
      
      # if it can make a new tree
      if(!(ltpar==1 & lcchild==csize & cindex!=0)){
        lnom <- dpolya1(tpar, neutral[tindex], tau, prior)
        lden <- dpolya1(partition[tindex], neutral[tindex], tau, prior)
        if(ltpar > 2  ){
          bzero<-runif(1)
          if(bzero > zerop[ltpar] & ((lcchild !=csize)||(cindex==0))){      
            epar <- sample(unique(tpar), size=1)
            tindex1 <- tindex[tpar!=epar]
            tpar <- tpar[tpar!=epar]
            tp <- rep(epar, length(tindex1))
            np[tindex1] <- struc(tindex1, unique(tpar), tp, tpar,res,zerop) 
          }
          else np[tindex]<-struc(tindex, unique(tpar), tp, tpar,res,zerop)
        }
        nlpa <- sum(is.element(np, cindex)) 
        ncnum<-length(unique(npar[np==cindex]))
        plogpri <- dpolya1(npar,neutral,tau, prior) #p(old|new)-p(new|old)
        ploglik <- Tlik(para, npar, np, neutral,data )
        logmh <- plogpri-logpri +lden -lnom + ploglik - loglik
        logmh<-logmh-log(lpa*(1/lcchild)*(1/choose(lcchild, csize)))
        logmh<-logmh+log(nlpa*(1/ncnum)*(1/choose(ncnum,(lcchild-csize))))
       }
     }

    ok1 <- !is.na(logmh)
    if( ok1 ) {
      accept <- ( log(runif(1)) <= logmh )
      if( accept ) {
         parent <- np; partition<-npar
         logpri <- plogpri; loglik<-ploglik
         acrate[1] <- acrate[1]+1
      }
    }

   ##################################################################
   # Update II: shuffle the partition; choose nprem aberrations and 
   # shuffle them.
   ###################################################################
    logmh <- NA
    ii <- sample(c(1:nswitch), size=nperm, replace=F)
    jj <- sample(ii)  # permute these
    npar <- partition; np<-parent
    npar[ii] <- partition[jj]; np[ii]<-parent[jj]
    # Check whether or not to accept if proposal != partition
    if( any( npar != partition ) ){
      ploglik <-Tlik(para, npar, np, neutral, data )
      logmh <- ploglik-loglik
      ok <- !is.na(logmh)
      if( ok ){
        accept <- ( log(runif(1)) <= logmh )
        if( accept ){
           parent<-np; partition<-npar
           loglik <-  ploglik
           acrate[2] <- acrate[2]+1
        }
      }
    }
    ###################################################################
    # Update method III:
    # Change the activation status (active/null) of a random set
    ####################################################################
    change <- sample(c(1:nswitch), size=1) ;logmh<-NA
    pneutral <- neutral; ok<-T
    while(ok){
      pneutral[change] <- !neutral[change]
      if(!all(pneutral)==T) ok<-F 
      else change <- sample(c(1:nswitch), size=1) 
    }
    ploglik <- Tlik(para, partition, parent, pneutral, data)
    plogpri <- dpolya1(partition, pneutral, tau, prior)
    logmh <- plogpri-logpri+ploglik-loglik
    ok <- !is.na(logmh)
    if( ok ){
      accept <- ( log(runif(1)) <= logmh )
      if( accept ){
         neutral <- pneutral
         logpri <- plogpri
         loglik <- ploglik
         acrate[3] <- acrate[3]+1
       }
     }
    ##################################################################
    # Update method IV:
    # Change the activation status of a whole path (of length >= 2)
    ##################################################################
    logmh<-NA
    path <- sample( partition, size=1 )
    ind <- (partition == path) 
    mm <- sum(ind)
    if( mm > 1 ){
      ok<-T
      while(ok){ 
        nullpart <- neutral[ind]
        mm0 <- sum( nullpart )
        lqnum <- log( padII )   # log( q(state* , state ) )
        if( mm0 > 0 & mm0 < mm )
          lqnum <- log( 1 - 2*padII ) - log( 2^mm - 2 ) 
        xi <- runif(1)
        if( xi <= padII ){  # fully activate path (make sure padII < 1/2)
          nullnew <- rep(F,mm)
          lqden <- log( padII )
        }
        else if( xi <= 2*padII ){  # fully deactivate path
          nullnew <- rep(T,mm)
          lqden <- log( padII )
        }
        else{  # sample one of the 2^mm - 2 mixed arrangements
          nullnew <- mixsample( mm )
          lqden <- log( 1 - 2*padII ) - log( 2^mm - 2 ) 
        } 
        nullstar <-  neutral
        nullstar[ind] <- nullnew
        if(!all(nullstar)==T) ok<-F
      }
      ploglik <- Tlik(para, partition, parent, nullstar,data )
      plogpri <- dpolya1(partition,nullstar,tau, prior)
      logmh <- plogpri-logpri+lqnum-lqden +ploglik-loglik
      ok <- !is.na(logmh)
      if( ok ){
         accept <- ( log(runif(1)) <= logmh )
         if( accept ){
           neutral <- nullstar
           logpri <- plogpri
           loglik <- ploglik
           acrate[4] <- acrate[4]+1
         }
       }
    }   

    ##################################################################
    # Update V: change the tree structure; change the position of sub tree.
    # choose one parent and its whole child edges. Put the whole sub-tree
    # to a new position but reserve the original structure of sub-tree. 
    ########################################################################
    logmh <- NA
    if(!all(parent==0)){
      cindex <- sample(setdiff(unique(parent),0), size=1);  
      cchild <- unique(partition[parent==cindex])
      ccandi <- cchild
      leaf<-setdiff(partition, parent)
      ok<-T; tmp<-NULL; 
      # choose whole child edges. 
      while(ok){
        tmp<-partition[is.element(parent,cchild)]
        ttmp<-setdiff(tmp, cchild)
        cchild<-unique(c(cchild,tmp))
        if(sum(!is.element(ttmp, leaf))<=0) ok<-F;
      }
      tindex<-index[is.element(partition, ccandi)]
      poss<-setdiff(leaf, cchild)
      if(length(poss)!=0){
        rep<-sample(poss,size=1)
        np<-parent
        np[tindex]<-rep
        ploglik <- Tlik(para, partition, np, neutral, data)
        logmh<-ploglik-loglik
        ok <- !is.na(logmh)
        if( ok ){
          accept <- ( log(runif(1)) <= logmh )
          if( accept ){
            parent <- np
            loglik <- ploglik
            acrate[5] <- acrate[5]+1
          }
        }
      }
    }
    ##################################################################
    # Update VI: change the tree structure; 
    # choose two edges and change their position. 
    ###################################################################
    logmh <- NA
    lup<-length(unique(partition))
    if(lup>2){
      ii <- sample(unique(partition), size=2, replace=F)
      fpa <- unique(parent[partition==ii[1]])
      spa <- unique(parent[partition==ii[2]])
      leaf <- setdiff(partition, parent)
      if(!(all(is.element(ii,leaf)) & fpa==spa)){      
        fg<-index[partition==ii[1]]
        sg<-index[partition==ii[2]]
        # change the position of partition and parent
        npar <- partition; np<-parent
        npar[fg] <- ii[2]; np[fg]<-spa
        npar[sg] <- ii[1]; np[sg]<-fpa
        ploglik <- Tlik(para, npar, np, neutral, data )
        logmh <- ploglik-loglik
        ok <- !is.na(logmh)
        if( ok ){
          accept <- ( log(runif(1)) <= logmh )
          if( accept ){
             parent<-np; partition<-npar
             loglik <-  ploglik;
             acrate[6] <- acrate[6]+1
          }
        }
      }
    }

    #################################################################
    # Update VII: change the tree structure; choose two partition edges, 
    # and make new two partition edges.  
    ##############################################################
    upar<-unique(partition)
    if(length(upar)>=2){
      # choose two partition edges; e1, e2
      ii <- sample(upar, size=2, replace=F)
      pindex<-index[is.element(partition,ii)]
      lpindex<-length(pindex)
      if(lpindex>2){
        np<-parent; npar<-partition
        # at least one aberration should assigned in e1 and e2.
        k<-sample(ii, size=lpindex-2, replace=T)
        tpar<-sample(c(ii,k))
        npar[pindex]<-tpar
        tp<-parent[pindex]
        upar<-partition[pindex]
        np[pindex]<-tp[match(tpar, upar)]
        plogpri <- dpolya1(npar,neutral,tau, prior)
        ploglik <-Tlik(para, npar, np, neutral, data )
        logmh <- plogpri-logpri +ploglik-loglik
        ok <- !is.na(logmh)
        if( ok ){
          accept <- ( log(runif(1)) <= logmh )
          if( accept ){
            parent<-np; partition<-npar
            logpri<-plogpri; loglik <-  ploglik
            acrate[7] <- acrate[7]+1
          }
        }
      }
    }
    ######################################################################
    # Update VIII: change the tree structure; choose one avaiable leaf
    # edge whose parent edge has more then 3 child edges and change its
    # position.  
    #######################################################################
    if(!all(parent==0)){   
      kkk<-match(unique(partition), partition)
      leaf<-setdiff(partition, parent)
      eindex<-is.element(partition,leaf)
      up<-unique(parent)
      k<-outer(up, parent, "==")
      t<-as.matrix(k[,kkk]); tt<-eindex[kkk]
      nchild <- apply(t,1,sum)
      elength <- apply(t[,tt],1,sum)
      target<-up[nchild>=3 & elength>=1]
      if(nchild[up==0]>=2 & elength[up==0]>=1) target<-c(target,0)
      if(length(target)>0){
        pindex<-index[is.element(parent,target)] 
        peindex<-eindex[is.element(parent,target)]
        poss<-intersect(pindex[peindex],kkk)
        ldom<-length(poss)
        if(ldom==1) pick<-poss
        else pick<-sample(poss, size=1)
        pick<-index[partition==partition[pick]]
        np<-parent
        np[pick]<-sample(setdiff(up, parent[pick]), size=1)
        nleaf<-setdiff(partition, np)
        eindex<-is.element(partition,nleaf)
        up<-unique(np)
        k<-outer(up, np, "==")
        t<-as.matrix(k[,kkk]); tt<-eindex[kkk]
        nchild<-apply(t,1,sum)
        elength<-apply(t[,tt],1,sum)
        target<-up[nchild>=3 & elength>=1]
        if(nchild[up==0]>=2 & elength[up==0]>=1) target<-c(target,0)
        npindex<-index[is.element(np,target)]
        npeindex<-eindex[is.element(np,target)]
        nposs<-intersect(npindex[npeindex],kkk)
        lnom<-length(nposs)
        ploglik <-Tlik(para, partition, np, neutral, data )
        logmh<-log(ldom/lnom) +ploglik-loglik #p(old|new)-p(new|old)
        ok <- !is.na(logmh)
        if( ok ){
          accept <- ( log(runif(1)) <= logmh )
          if( accept ){
            parent<-np; 
            loglik <-  ploglik
            acrate[8] <- acrate[8]+1
          }
        }
      }
    }       
    #################################################################
    #
    # Store summary statistics periodically..
    if(skipcount == nskip ){
      skipcount <- 0
      psave[isave,] <- partition; 
      pasave[isave,] <- parent
      pneu[isave,]<-neutral
      liksave[isave] <- loglik
      isave <- isave+1
    }
  }
  post <- liksave; index<-c(1:nsave)
  netsave <- index[max(post)==post][1]
  mtree<-ematrix(psave[netsave,],pneu[netsave,],pasave[netsave,],data)
  result<-list(partition=psave, parent=pasave, neutral=pneu,
     loglik=liksave, acrate=acrate, data=data, mtree=mtree$tree)
  return( result )
}

"SISTN" <- function(data,mcmc,tau,para,res){
  #########################################################################  
  # data = binary matrix, ntumor x nswitch, holding CGH profiles
  # mcmc = list with driving instructions
  #        nskip  = size of stretch between saved states
  #        nsave  = number of saved states
  #        nperm  = number of switches permuted in 2nd movetype update
  #        padI =  proportion of switches to change in first neutral vector
  #                update (3rd movetype)
  #        padII (<1/2) = chance of full activation / deactivation in second
  #                       neutral vector update (4th movetype)
  # tau > 0.  The cluster prior hyperparameter
  # para      Starting values for theta, gamma, delta
  # res       Probability matrix used for making a sub-tree.
  #           Refer function res.
  # 
  # ex.
  # source('tree.R')
  # res<-dget('res') # it is highly recommended to save res<-res(20)
  # mcmc <- list(nsave=10, nskip=100, nperm=3, padI=.15, padII=.15)
  # dd<-dget('trim'); data<-dd$gradeI; 
  # para<-c(0.05, 0.01, 0.01);tau<-4
  # mh <- ISTN(data,mcmc,tau,para,res)
  ##########################################################################
  nsave<-mcmc$nsave; nskip<-mcmc$nskip; nperm<-mcmc$nperm;
  padI<-mcmc$padI; padII<-mcmc$padII; 
  nswitch<-ncol(data)
  # given n edges, the number of possible trees.
  # one can get these numbers using function "pri" but it is highly recommended
  # to save the result and use it.
  # prior<-rep(0,25); for(i in 1:25){prior[i]<-pri(i)}
  prior<-c(1,1,4,17,116,997,10270,115257,1539928, 23414021,3.911842e+08, 
     7.173987e+09,1.440993e+11, 3.170871e+12, 7.486385e+13, 1.895879e+15,
     5.117341e+16, 1.475958e+18,4.497267e+19, 1.450948e+21, 4.923465e+22, 
     1.760436e+24, 6.581871e+25, 2.579816e+27,1.053847e+29)
  # given n edges, the number of possible typeII trees; a tree composed of n 
  # edges and whose root has more than two child edges. 
  # rres<-p1(n, NULL)$ans
  rres <- c(0,1,1,13,51,691,5433,71793,893791,14476111,2.319470e+08, 
     4.390623e+09,8.702118e+10, 1.952574e+12, 4.557523e+13,1.166675e+15, 
     3.133992e+16, 9.118389e+17,2.764773e+19, 8.979930e+20,3.037680e+22, 
     1.092146e+24, 4.069935e+25,1.603032e+27, 6.530885e+28)
  # a probability being a typeII tree. 
  zerop <- c(rres/prior)
  prior<- -log(prior)
  index<-c(1:nswitch)

  # Create objects to save MCMC output.
  nscan <- nsave*nskip
  pasave <- matrix(NA,nsave,nswitch);  # saved parents
  psave <- matrix(NA,nsave,nswitch);   # saved partitions
  pneu <- matrix(NA,nsave,nswitch);    # saved neutral
  liksave <- numeric(nsave);      # saved log-likelihood 
  acrate <- rep(0,4);           # acceptance rates
  
  # make an inital tree
  init <- Tnet(nswitch, tau, zerop, res)
  partition <- init$x
  parent <- init$tree
  neutral <- init$neutral
 
  # Evaluate the loglikelihood and logprior of starting state
  loglik <- Tlik(para, partition, parent, neutral,data)
  logpri <- dpolya1(partition,neutral,tau, prior)

  ####################################################################
  # Run MCMC.
  # Propose a new network configuration (8 move types)
  ####################################################################

  skipcount <- 0
  isave <- 1
  for(iscan in 1:nscan){
    skipcount <- skipcount+1; 

    #######################################################################
    # Update I: change the tree structure; choose one parent edge and among
    # its child edges, choose some and change the whole structure. 
    ########################################################################
    logmh <- NA
    cindex <- sample(parent, size=1);  
    cchild <- unique(partition[parent==cindex])
    lcchild <- length(cchild)
    lpa <- sum(is.element(parent, cindex)) 
    csize<-sample(c(1:lcchild), size=1) 
    # choose some of child edges
    ccandi<- sample(cchild, size=csize) 
    leaf<-setdiff(partition, parent)
    ok<-T; tmp<-NULL; lc<-1
    while(ok){
      tmp<-partition[is.element(parent,ccandi)]
      ccandi<-unique(c(ccandi,tmp))
      nlc<-length(ccandi)
      if((nlc-lc)>0) lc<-nlc
      else ok<-F
    }
    tindex<-index[is.element(partition, ccandi)]
    lden<-0; lnom<-0; 
    ltindex <- length(tindex)
    if(ltindex>1){
      npar <- partition; 
      rok<-T
      while(rok){
        npar[tindex]<-runif(ltindex); 
        if(sum(is.element(npar[tindex], partition))==0) rok<-F
      }
      epartition <- npar[tindex[1]]
      np<-parent; np[tindex]<-cindex
      for( i in 2:ltindex){
        pwild <- tau/(tau+i-1)
        if( runif(1) > pwild ){ # don't take the new label
          if(length(epartition)==1) npar[tindex[i]]<-epartition
          else{
            sam <- sample( epartition, size=1 ) 
            npar[tindex[i]] <- sam 
          }
        }
        epartition <- c(epartition, npar[tindex[i]] )
      }#for end
      tp<-np[tindex]
      tpar<-npar[tindex]
      ltpar <- length(unique(tpar))
      
      # if it can make a new tree
      if(!(ltpar==1 & lcchild==csize & cindex!=0)){
        lnom <- dpolya1(tpar, neutral[tindex], tau, prior)
        lden <- dpolya1(partition[tindex], neutral[tindex], tau, prior)
        if(ltpar > 2  ){
          bzero<-runif(1)
          if(bzero > zerop[ltpar] & ((lcchild !=csize)||(cindex==0))){      
            epar <- sample(unique(tpar), size=1)
            tindex1 <- tindex[tpar!=epar]
            tpar <- tpar[tpar!=epar]
            tp <- rep(epar, length(tindex1))
            np[tindex1] <- struc(tindex1, unique(tpar), tp, tpar,res,zerop) 
          }
          else np[tindex]<-struc(tindex, unique(tpar), tp, tpar,res,zerop)
        }
        nlpa <- sum(is.element(np, cindex)) 
        ncnum<-length(unique(npar[np==cindex]))
        plogpri <- dpolya1(npar,neutral,tau, prior) #p(old|new)-p(new|old)
        ploglik <- Tlik(para, npar, np, neutral,data )
        logmh <- plogpri-logpri +lden -lnom + ploglik - loglik
        logmh<-logmh-log(lpa*(1/lcchild)*(1/choose(lcchild, csize)))
        logmh<-logmh+log(nlpa*(1/ncnum)*(1/choose(ncnum,(lcchild-csize))))
       }
     }

    ok1 <- !is.na(logmh)
    if( ok1 ) {
      accept <- ( log(runif(1)) <= logmh )
      if( accept ) {
         parent <- np; partition<-npar
         logpri <- plogpri; loglik<-ploglik
         acrate[1] <- acrate[1]+1
      }
    }

   ##################################################################
   # Update II: shuffle the partition; choose nprem aberrations and 
   # shuffle them.
   ###################################################################
    logmh <- NA
    ii <- sample(c(1:nswitch), size=nperm, replace=F)
    jj <- sample(ii)  # permute these
    npar <- partition; np<-parent
    npar[ii] <- partition[jj]; np[ii]<-parent[jj]
    # Check whether or not to accept if proposal != partition
    if( any( npar != partition ) ){
      ploglik <-Tlik(para, npar, np, neutral, data )
      logmh <- ploglik-loglik
      ok <- !is.na(logmh)
      if( ok ){
        accept <- ( log(runif(1)) <= logmh )
        if( accept ){
           parent<-np; partition<-npar
           loglik <-  ploglik
           acrate[2] <- acrate[2]+1
        }
      }
    }
    ###################################################################
    # Update method III:
    # Change the activation status (active/null) of a random set
    ####################################################################
    change <- sample(c(1:nswitch), size=1) ;logmh<-NA
    pneutral <- neutral; ok<-T
    while(ok){
      pneutral[change] <- !neutral[change]
      if(!all(pneutral)==T) ok<-F 
      else change <- sample(c(1:nswitch), size=1) 
    }
    ploglik <- Tlik(para, partition, parent, pneutral, data)
    plogpri <- dpolya1(partition, pneutral, tau, prior)
    logmh <- plogpri-logpri+ploglik-loglik
    ok <- !is.na(logmh)
    if( ok ){
      accept <- ( log(runif(1)) <= logmh )
      if( accept ){
         neutral <- pneutral
         logpri <- plogpri
         loglik <- ploglik
         acrate[3] <- acrate[3]+1
       }
     }
    ##################################################################
    # Update method IV:
    # Change the activation status of a whole path (of length >= 2)
    ##################################################################
    logmh<-NA
    path <- sample( partition, size=1 )
    ind <- (partition == path) 
    mm <- sum(ind)
    if( mm > 1 ){
      ok<-T
      while(ok){ 
        nullpart <- neutral[ind]
        mm0 <- sum( nullpart )
        lqnum <- log( padII )   # log( q(state* , state ) )
        if( mm0 > 0 & mm0 < mm )
          lqnum <- log( 1 - 2*padII ) - log( 2^mm - 2 ) 
        xi <- runif(1)
        if( xi <= padII ){  # fully activate path (make sure padII < 1/2)
          nullnew <- rep(F,mm)
          lqden <- log( padII )
        }
        else if( xi <= 2*padII ){  # fully deactivate path
          nullnew <- rep(T,mm)
          lqden <- log( padII )
        }
        else{  # sample one of the 2^mm - 2 mixed arrangements
          nullnew <- mixsample( mm )
          lqden <- log( 1 - 2*padII ) - log( 2^mm - 2 ) 
        } 
        nullstar <-  neutral
        nullstar[ind] <- nullnew
        if(!all(nullstar)==T) ok<-F
      }
      ploglik <- Tlik(para, partition, parent, nullstar,data )
      plogpri <- dpolya1(partition,nullstar,tau, prior)
      logmh <- plogpri-logpri+lqnum-lqden +ploglik-loglik
      ok <- !is.na(logmh)
      if( ok ){
         accept <- ( log(runif(1)) <= logmh )
         if( accept ){
           neutral <- nullstar
           logpri <- plogpri
           loglik <- ploglik
           acrate[4] <- acrate[4]+1
         }
       }
    }   
    #################################################################
    #
    # Store summary statistics periodically..
    if(skipcount == nskip ){
      skipcount <- 0
      psave[isave,] <- partition; 
      pasave[isave,] <- parent
      pneu[isave,]<-neutral
      liksave[isave] <- loglik
      isave <- isave+1
    }
  }
  post <- liksave; index<-c(1:nsave)
  netsave <- index[max(post)==post][1]
  mtree<-ematrix(psave[netsave,],pneu[netsave,],pasave[netsave,],data)
  result<-list(partition=psave, parent=pasave, neutral=pneu,
     loglik=liksave, acrate=acrate, data=data, mtree=mtree$tree)
  #result<-list(partition=psave, parent=pasave, neutral=pneu,
  #   loglik=liksave, acrate=acrate)
  return( result )
}

"Tlik" <- function(rates,partition,parent,neutral,data) {  
 #######################################################################
 # Calculate the log joint probability density of aberrations in the 
 # Instability Selection Tree-like Network model: 
 #
 # rates = (theta,delta, gamma)  in (0,1)^3 
 #           theta = Prob[ Aberration ]
 #           gamma = Prob[ X = 1 | no Aberration ]
 #           delta = Prob[ X = 0 |  Aberration ]
 # data = ntumor x nswitch matrix of binaries 0/1 indicating abnormality
 #            length should be positive.
 # partition = ncomponent vector holding distinct edge labels
 # parent    = like partition, but holds labels for parent edges
 #             0 is code for root node
 # neutral = logical vector indicating neutral components
 #######################################################################

  #theta <- rates[1]; delta <- rates[2]; gamma <- rates[3]
  #alpha <- (theta)*(1-delta) + (1-theta)*gamma  
  # marginal frequency; neutrals
  #beta1 <- theta*delta/(1-alpha)
  #beta2 <- theta*(1-delta)/alpha
  alpha <- rates[1]; beta2<- rates[2]; beta1<-rates[1]
  theta <- beta1*(1-alpha)+beta2*alpha
 
  # First get the prior part of the log density
  s0 <- sum( data==0 );  s1 <- sum( data==1 )
  logprob <- s0*log( 1-alpha[1] ) + s1*log( alpha[1] )

  # In case all switches are in the neutral path we're done; otherwise
  if( !( all( neutral ) ) ) {
    apart <- partition[!neutral]; anse <- apart
    tmp <- NULL; ok <- T; lc <- 1
    while(ok){
      tmp<-parent[is.element(partition,apart)]
      apart<-unique(c(apart,tmp))
      nlc<-length(apart)
      if((nlc-lc)>0) lc<-nlc
      else ok<-F
    }
    k<-ifelse(is.element(partition, setdiff(apart,anse)), F,T)
    tneutral<-k*neutral # for non-neutral branch: make F, 
                        # eleminate it(make T) when all leaf ==T(none) 
    npart <- partition[!tneutral]
    nparent <- parent[!tneutral]
    tmpin <- outer( npart, unique(npart), "==" )
    nedge <- ncol(tmpin)       # number of edges
   
    edgelengths <- c( rep(1,sum(!tneutral)) %*% tmpin )
    edgelabs <- unique( npart )
    tmp <- match( edgelabs, npart )
    eparent <- nparent[tmp]     # parent information at edge-level

    # degrees of nodes in the tree
    children <- outer( eparent, edgelabs, "==" )  # columns hold children
                                                 # of each edge
    nchild <- c( rep(T,nedge) %*% children )
    adat <-  data[, !tneutral]        # active data
    stats1 <- (1-adat) %*% tmpin  # this is ntumor x nedge
    stats2 <- (adat) %*% tmpin  # this is ntumor x nedge

    tmps <- ifelse(is.element(unique(npart), anse), F,T) 
    ntumor <- nrow(stats1)
   
    teprob.close <- log1p( -(beta1^stats1)*(beta2^stats2) ) 
                # log Prob[ edge closed |x ]
    teprob.close <- rbind( teprob.close, log1p( -theta^edgelengths ) )
                  # This deals with unconditional case also
    eprob.close<-matrix(-Inf, (ntumor+1), nedge)
    eprob.close[,!tmps]<-teprob.close[,!tmps]
    bprob.close <- eprob.close # log Prob[not SEL on branch|x] ; or not cond.
   
    # Any edge which is not a parent of another edge should be
    # dealt with first. (selection from that edge = prob that edge open )
   
    finished <- (nchild==0)
    while( !all(finished) ) {
      # Which edges are being updated on this iteration?
      tmp <- c(finished %*% children)
      allcdone <- (tmp==nchild)  # all children finished?
      beingdone <- ( allcdone & !finished )
                                                 
      zz <- ( bprob.close %*% children )[,beingdone]
      yy <- ( eprob.close )[,beingdone]
      mm <- pmax(zz,yy)
      kkk<- expm1(yy-mm) + expm1(zz-mm) - expm1(yy+zz-mm)
      tmp<- mm+log1p(kkk) 
      tmp[is.nan(tmp)] <- -Inf    # correct when (-inf + inf)

      bprob.close[,beingdone] <- tmp
      finished[beingdone] <- T
    }
                  
    # Now put it all together.
    fromroot <- as.matrix( eparent == 0 )
    tmp <- bprob.close[1:ntumor,] %*% fromroot
    lnumer <- sum( log1p( -1-expm1(tmp)) )              #log P{SEL|data}
    tmp <- sum( bprob.close[ntumor+1,(eparent==0)]  )
    ldenom <- ntumor*log1p( -1-expm1(tmp) )            # log P{SEL}
    logprob <- logprob + lnumer-ldenom
  }
  return(as.numeric(logprob))
}


"dpolya1" <- function(partition,neutral,tau, prior) {

  #####################################################################
  #
  # Evaluate the log prior mass of a Double Polya Prior-augmented state.
  #
  # This is the prior for an augmented-state partition in ISTN model
  # (Used in the data-augmentation part of the MCMC algorithm.)
  #
  # partition = n vector holding distinct edge labels
  # neutral = logical vector indicating relevant(1)/neutral(0) aberration
  # tau = clustering parameter: bigger <-> more edges
  # prior= probability of choosing one tree structure given number of edge.
  #####################################################################

  nswitch <- length(partition)
  upar<-unique(partition)
  logp <- prior[length(upar)]

  # The partition componenent
  incidence <- outer( partition, upar, "==" )
  sizes <- t(incidence) %*% rep(1,nswitch)
  npath <- length(sizes)
  logp <-  logp + npath*log(tau) + sum( lgamma(sizes) ) 
  logp <- logp + lgamma(tau)-lgamma(nswitch+tau)

  # The activation component
  nnull <- sum( neutral )
  logp <- logp+lgamma(nnull+1)+lgamma(nswitch-nnull+1)-lgamma(nswitch+2)
  return( logp )
 }

"mixsample" <- function(n) {
  ######################################################################
  # n >= 2
  # Returns a random vector of T's and F's uniformly among those
  # that are not all T or all F
  # Method: rejection sampling
  ######################################################################
  notdone <- T
  while( notdone ) {
    x <- sample( c(T,F), size=n, replace=T )
    notdone <- (sum(x)==0) | (sum(x) == n)
  }
  return(x)
}

"Trisn" <- function(rates,partition,parent,neutral,numtumor){ 
  ####################################################################
  #  Sample random binary vectors from the
  #  Instability-Selection-Tree like Network Model.
  #  rates = (theta,delta, gamma)  in (0,1)^3 
  #           theta = Prob[ Aberration ]
  #           gamma = Prob[ X = 1 | no Aberration ]
  #           delta = Prob[ X = 0 |  Aberration ]
  #  partition = n component vector holding distinct edge labels
  #              0 is coded for neutral aberration.
  #  parent = like partition, but holds labels for parent edges
  #           0 is coded for root node and -1 is coded for neutral aberration.
  #  neutral = logical vector indicating neutral components
  #  numtumor = the number of tumors to be simulated
  #
  #  Method: Rejection sampling
  #####################################################################
  #theta <- rates[1]; delta <- rates[2]; gamma <- rates[3] 
  #alpha <- theta*(1-delta) + (1-theta)*gamma 
  #beta1 <- theta*delta/(1-alpha) 
  #beta2 <- theta*(1-delta)/alpha 
  alpha <- rates[1]; beta2 <- rates[2]; beta1 <- rates[3] 
  nswitch <- length(partition) 
  apart <- partition[!neutral]    # active partition 
  aparent <- parent[!neutral] 
  incidence <- outer( apart, unique(apart), "==" ) 
  nedge <- ncol(incidence)       # number of edges 
  edgelengths <- c( rep(1,sum(!neutral)) %*% incidence ) 
  edgelabs <- unique( apart ) 
  tmp <- match( edgelabs, apart ) 
  eparent <- aparent[tmp]     # parent information at edge-level 
  # degrees of nodes in the tree 
  children <- outer( eparent, edgelabs, "==" ) 
  nchild <- c( rep(T,nedge) %*% children ) 
  xsave <- matrix(NA,numtumor,nswitch) 
  for( i in 1:numtumor ){ 
    notdone <- T 
    while( notdone ){
      xi <- ifelse( runif(nswitch) <= alpha, 1, 0 ) 
      adat <-  xi[!neutral]        # active data 
      stats1 <- (1-adat) %*% incidence  # this is ntumor x nedge 
      stats2 <- (adat) %*% incidence  # this is ntumor x nedge 
      ntumor <- nrow(stats1) 
      eprob.close <- log1p( -(beta1^(stats1))*(beta2^stats2) ) 
      bprob.close <- eprob.close 
      finished <- (nchild==0) 
      while( !all(finished) ){
        # Which edges are being updated on this iteration? 
        tmp <- c(finished %*% children) 
        allcdone <- (tmp==nchild)  # all children finished? 
        beingdone <- ( allcdone & !finished ) 
        zz <- ( bprob.close %*% children )[,beingdone] 
        yy <- ( eprob.close )[,beingdone] 
        mm <- pmax(zz,yy) 
        kkk<- expm1(yy-mm) + expm1(zz-mm) - expm1(yy+zz-mm)  
        tmp<- mm+log1p(kkk)
        tmp[is.nan(tmp)] <- -Inf    # correct when (-inf + inf) 
        bprob.close[,beingdone] <- tmp 
        finished[beingdone] <- T 
      }
      # Now put it all together. 
      fromroot <- as.matrix( eparent == 0 ) 
      tmp <- bprob.close[1:ntumor,] %*% fromroot 
      env <- sum( log1p(-1 -expm1(tmp)) ) 
      notdone <- ifelse( log(runif(1)) <= env, F, T ) 
    }
    xsave[i,] <- xi 
  } 
  return(xsave) 
}

"Tnet" <- function(n,tau, zerop, res) {
  ##############################################################
  # Simulate a tree-like network.
  # zerop = a probability being a typeII tree. 
  # n = number of objects to be partitioned
  # tau > 0 = hyperparameter
  # zerop = # of typeI tree / #of possible trees
  #         it is recommended to save zerop and use it. 
  #         prior<-rep(0,25); rres<-rep(0,25)
  #         for(i in 1:25){
  #           prior[i]<-pri(i)}
  #           rres[i]<-p1(n, NULL)$ans
  #         }
  #         zerop <- c(rres/prior) 
  # 
  # res = Probability matrix used for making a sub-tree.
  #       Refer to function res
  #############################################################
  theta <- runif(1); ok<-T
  while(ok){
    neutral <- sample(c(1,0), n, replace = TRUE, prob = c(theta, 1-theta))
    if(sum(!neutral)>0) ok<-F
  }
  v <- runif(n)
  x <- v
  if( n >= 2 ){
     for( i in 2:n ){
       pwild <- tau/( tau + i - 1 )
       if( runif(1) > pwild ){ # select an old one
         ind <- sample( 1:(i-1), size=1 )
         x[i] <- x[ind]
       }
     }
   }
   
   pa<-rep(0,n); tindex<-c(1:n)
   ltpar <- length(unique(x))
   nx<-x; 
   if(ltpar > 2){ 
     bzero<-runif(1)
     if(bzero > zerop[ltpar]){      
       epar <- sample(unique(nx), size=1)
       tindex <- tindex[epar!=nx]
       nx <- nx[epar!=nx]
       npa <- rep(epar, length(tindex))
       pa[tindex] <- struc(c(1:length(tindex)), unique(nx), npa, nx, res,zerop) 
       tree<-pa
     }
     else tree<-struc(c(1:length(tindex)), unique(x), pa, x, res,zerop)
   }
   else tree<-pa
   result<-list(neutral=neutral, x=x, tree=tree)
   return(result)
 }

"Tnet1" <- function(n,tau, zerop, res) {
  ############################################################# 
  # Simulate a tree-like network; Same with function maketree except
  #  it only generate partition and parent vector. 
  #############################################################  
  v <- runif(n)
  x <- v
  if( n >= 2 ){
     for( i in 2:n ){
       pwild <- tau/( tau + i - 1 )
       if( runif(1) > pwild ){ # select an old one
         ind <- sample( 1:(i-1), size=1 )
         x[i] <- x[ind]
       }
     }
   }
   pa<-rep(0,n); tindex<-c(1:n)
   ltpar <- length(unique(x))
   nx<-x; 
   if(ltpar > 2){ 
     bzero<-runif(1)
     if(bzero > zerop[ltpar]){      
       epar <- sample(unique(nx), size=1)
       tindex <- tindex[epar!=nx]
       nx <- nx[epar!=nx]
       npa <- rep(epar, length(tindex))
       pa[tindex] <- struc(c(1:length(tindex)), unique(nx), npa, nx, res,zerop) 
       tree<-pa
     }
     else tree<-struc(c(1:length(tindex)), unique(x), pa, x, res,zerop)
   }
   else tree<-pa
   result<-list(x=x, tree=tree)
   return(result)
 }

struc<-function(tindex, tempar, partp, partpa, res,zerop){
  ########################################################################
  # Among N edges (N>3), choose m edges using 'res' matrix and make a 
  # sub-tree. Recursively call struc function and make sub-tree again
  # using N-m edges.
  # 
  # tindex = index of aberrations being chosen to make a new sub-tree 
  # structure.
  # tempar = unique(partpa)
  # partp = sub vector of parent. 
  # partpa =  sub vector of partition. 
  # res = pre-calculated  probability matrix from function res.
  ########################################################################

  nbranch<-length(tempar)
  r<-runif(1); r<-rep(r, nbranch); rindex<-c(1:nbranch)
  re<-res[nbranch,c(1:nbranch)]
  ngroup<-min(rindex[r<re]);
  ftempar<-sample(tempar, size=ngroup)
  findex<-tindex[is.element(partpa, ftempar)]; fuindex<-findex; 

  if(nbranch>3){      
    fp<-partp[is.element(tindex, findex)]
    if((ngroup<nbranch) && ngroup>=3){
      fh<-sample(ftempar, size=1)
      ftg <-setdiff(ftempar,fh) 
      ftindex<-tindex[is.element(partpa, ftg)]
      ftp<-rep(fh, length(ftindex))
      fpa<-partpa[is.element(partpa, ftg)]
      fp[is.element(findex, ftindex)]<-struc(ftindex, ftg, ftp, fpa, res,zerop);
    }
    tempar<-setdiff(tempar, ftempar)
    tm<-fp
    nbranch<-length(tempar)
      
    while(nbranch>2){
      r<-runif(1); r<-rep(r, nbranch); rindex<-c(1:nbranch)
      re<-res[nbranch,]
      ngroup<-min(rindex[r<re]); 
      ftempar<-sample(tempar, size=ngroup)
      findex<-tindex[is.element(partpa, ftempar)]
      fuindex<-c(fuindex, findex)
      fp<-partp[is.element(tindex, findex)]

      if((ngroup<=nbranch) && ngroup>=3){
        fh<-sample(ftempar, size=1)
        ftg <-setdiff(ftempar,fh) 
        ftindex<-tindex[is.element(partpa, ftg)]
        ftp<-rep(fh, length(ftindex))
        fpa<-partpa[is.element(partpa, ftg)]
        fp[is.element(findex, ftindex)]<-struc(ftindex, ftg, ftp, fpa, res,zerop);
      }
      else{
        fpa<-partpa[is.element(partpa, ftempar)]
        fp<-struc(findex, ftempar, fp, fpa, res,zerop);
      }

      tm<-c(tm, fp)
      tempar<-setdiff(tempar, ftempar)
      nbranch<-length(tempar) 
    }
    rest<- setdiff(tindex, fuindex)
    fp<-partp[is.element(tindex, rest)]; 
    fuindex<-c(fuindex, rest)
    tm<-c(tm, fp)
    rp<-tm[order(fuindex)]
    return(rp)
  }
  return(partp);
}

Tcluster <- function(mcmc,name){ 
  cov1 <- hcluster (mcmc, name)
  kkk<-dist(cov1)/max(dist(cov1)) 
  kkk<-hclust(kkk);
  plclust(kkk,main="Posterior probability clustering",xlab="",
          sub="", ylab="Height")
}
hcluster <-function(dd, name){
  m<-ncol(dd$parent); index<-c(1:m)
  cov<-matrix(0,m,m)
  for(i in 1:nrow(dd$partition)){
    partition<-dd$partition[i,]; parent<-dd$parent[i,]
    for(j in 1:m){
      apart<-partition[j]
      tmp<-NULL; ok<-T; lc<-1
      while(ok){
        tmp<-parent[is.element(partition,apart)]
        apart<-unique(c(apart,tmp))
        nlc<-length(apart)
        if((nlc-lc)>0) lc<-nlc
        else ok<-F
      }
      apart<-partition[j]
      tmp1<-NULL; ok<-T; lc<-1
      while(ok){
        tmp1<-partition[is.element(parent,apart)]
        apart<-unique(c(apart,tmp1))
        nlc<-length(apart)
        if((nlc-lc)>0) lc<-nlc
        else ok<-F
      }
      tmp2<-unique(c(tmp1, partition[j],tmp))
      ind<-index[is.element(partition,tmp2)]
      cov[j,ind]<-cov[j,ind]+1
    }
  }
  dimnames(cov)<-list(name, name)
  return(cov)
}

#######################################################################
# make a small box at each node. Used in function Tdraw.
#######################################################################
sw <- function(xcv,ylv,eps=.05, col="grey", txt=" "){
  for(i in 1:length(xcv)){
    xc<-xcv[i]; yl<-ylv[i]
    polygon( c(xc-eps/2, xc-eps/2,xc+eps/2,xc+eps/2), 
             c(yl-eps/2,yl+eps/2,yl+eps/2,yl-eps/2), col=col )
    lines( c(xc-eps/2, xc-eps/2 ), c(yl-eps/2,yl+eps/2), lwd=3 )
    lines( c(xc+eps/2, xc+eps/2 ), c(yl-eps/2,yl+eps/2), lwd=3 )
    lines( c(xc-eps/2, xc+eps/2 ), c(yl+eps/2,yl+eps/2), lwd=3 )
    lines( c(xc-eps/2, xc+eps/2 ), c(yl-eps/2,yl-eps/2), lwd=3 )
  }        
}

######################################################################
# Draw a tree-like network (top) and the number of tumors which 
# experiance (all aberrations reside on that edge are active)
# each edge (bottom). 
#
# partition = ncomponent vector holding distinct edge labels
# neutral = logical vector indicating neutral components
# parent    = like partition, but holds labels for parent edges
#             0 is code for root node
# data = ntumor x nswitch matrix of binaries 0/1 indicating abnormality
#            length should be positive.
#
# It eliminates the neutral aberrations and shrink the tree,
# if there is neutral edge (all aberrations on an edge are netral). Then
# it draw network using plot.tree and plot.stree functions.
#
# The number on each edge in the second figure indicates the number of 
# tumors whose all aberrations reside on that edge are active (not neutral)
# The figure at the bottom of each ensemble indicates the number of tumors
# whose all aberrations reside on that ensembes are active.
########################################################################

Tdraw <- function(partition,neutral,parent,data,line1,num1,letter1,
                  main1=NULL, main2=NULL){
  names<-dimnames(data)[[2]] 
  np<-parent; npa<-partition 
  if(any(neutral)){
    leaf <- setdiff(partition, parent); tl<-length(partition)
    apart <- partition[!neutral]
    incidence <- outer( apart, unique(partition), "==" )
    edgelengths <- c( rep(1,sum(!neutral)) %*% incidence )
    upartition <- unique(partition)
    lc <- upartition[edgelengths==0]
    lc1 <- setdiff(lc, leaf)
    if(length(lc1)>0){
      for(ii in 1:length(lc1)){
        np[is.element(np, lc1[ii])]<-unique(np[partition==lc1[ii]])
      }
    }
    lc2<-intersect(lc, leaf)
    if(length(lc2)>0){
      for(ii in 1:length(lc2)){
        kk<-npa[(np==np[npa==lc2[ii]])*(!neutral)==1]
        kk<-setdiff(kk,lc2[ii])
        if(length(kk)==1){
          if(unique(np[partition==kk])!=0){
            upp<-unique(np[npa==kk])
            if(all(neutral[np==upp])){
              ok<-T
              while(ok){
                upp<-unique(parent[partition==upp])
                ok<-all(neutral[partition==upp])
              }
            }
            npa[partition==kk]<- upp
            np[partition==kk]<-unique(np[partition==upp])
            if(!is.element(kk, leaf)) np[np==kk]<-upp
          }
        }
      }
    }
    np<-np[!neutral]; npa<-npa[!neutral]
    names<-names[!neutral]; data<-data[,!neutral]
  }
  par(mfrow=c(2,1), mar=c(2, 2, 2, 2))
  a<-plot.tree(npa,np,names,main=main1,line1,num1,letter1)
  iplot.tree(a$csave, a$inter, data, npa,np, names, main=main2,
             line1,num1)
}

###########################################################
# Function plot.tree draw the edges directly from the root. 
# It is used in function Tdraw. 
###########################################################
plot.tree<-function(partition,parent, names, main, line1,num1,letter1){
  upar<-unique(partition)
  index<-letters[1:length(upar)]
  n<-length(partition)
  leaf<-setdiff(partition, parent)
  nleaf<-length(leaf)
  layer<-1;
  ok <- !(nleaf==length(upar))
  xpos<-NULL; ypos<-NULL; nchild<-NULL
  ccandi <- unique(partition[parent==0])
  cindex<- index[is.element(upar, ccandi)]
  ctmp<-ccandi; csave<-NULL; inter<-NULL 
  nchild<-length(ctmp); csave[1]<-nchild
  tmp1<-outer(parent, ccandi, "==")
  a<-match(unique(partition), partition)
  csave[2]<-max(apply(as.matrix(tmp1[a,]),2,sum))
  while(ok){
    layer<-layer+1
    tmp<-partition[is.element(parent,ccandi)]
    ttmp<-setdiff(tmp, ccandi)
    ccandi<-unique(c(ccandi,tmp))
    tmp1<-outer(parent, ttmp, "==")
    csave[layer+1]<-max(apply(tmp1[a,],2,sum))
    if(sum(!is.element(ttmp, leaf))<=0) ok<-F;
  }
  inter<-rep(1, layer)
  if(layer>1){
    for(i in 1:(layer-1)){
      tmp<-inter[layer+1-i]
      inter[layer-i]<- tmp*csave[layer+1-i]
    }
  }
  xmax<-prod(csave[c(1:layer)]); ymax<-max(10, length(upar)+1)
  plot(1,1,type="n",xlab="", ylab="",axes=F,xlim=c(-1,xmax),ylim=c(-.6,ymax))
  if(!is.null(main)) mtext(main, side=3, line=1, cex=1.2)
  yint <- ymax/layer;  #xint <- xmax/nleaf; 
  if(nchild>1){
    xpos <- rep(0,nchild)+c(1:nchild)*inter[1]
    xpos<- xpos-mean(xpos)+xmax/2
  }
  else xpos<-xmax/2
  xx <- rep(xmax/2,2*nchild)
  xx[seq(2,2*nchild,2)]<-xpos
  yy<-rep(c(ymax,ymax-yint),nchild)
  lines(xx,yy, lty=2,col=4, lw=line1)
  sw(xx,yy)
  eps<- .2; cin1<-NULL; ctm1<-NULL; cin<-NULL; ctm<-NULL
  for(i in 1:nchild){
    xp <-xpos[i]+ (xmax/2-xpos[i])*.6 
    yp <- ymax-yint+yint*.6
    text( xp, yp, cindex[i], cex=num1,col=2)    
    cin1<-c(cin1, cindex[i]); ctm1<-c(ctm1, ctmp[i])
    ccandi <- unique(partition[parent==ctmp[i]])    
    ok<-T; tmp<-NULL; lc<-1
    while(ok){
      tmp<-partition[is.element(parent,ccandi)]
      ccandi<-unique(c(ccandi,tmp))
      nlc<-length(ccandi)
      if((nlc-lc)>0) lc<-nlc
      else ok<-F
    }
    kk <- is.element(partition, ccandi)
    if(sum(kk)>0){
      a<-plot.stree(NULL, NULL, upar,index,ctmp[i], partition[kk],parent[kk],
                   names[kk],xpos[i],yy[2],yint,inter[-1],line1,num1,letter1)
      cin<-c(cin, a$cin)
      ctm<-c(ctm, a$ctm)
    }
  }  
  cin<-c(unique(cin1), unique(cin)); ctm<-c(unique(ctm1),unique(ctm))
  for(ii in 1:length(cin)){
    text(-1,ymax-ii, paste(cin[ii],":",paste(names[partition==ctm[ii]],
                  collapse=","), collapse=","), adj=0, cex=letter1)
  }
  result<-list(csave=csave, inter=inter)
  return(result)
}   

######################################################################
# Function plot.stree draws the edges not directly from the root. 
# It recursively call itself. Used inside plot.tree, and tanal function.
########################################################################
plot.stree<-function(cin, ctm, upar, index,inip, partition, parent, names,
                     sx, sy, yint, inter, line1,num1,letter1){
  leaf<-setdiff(partition, parent)
  xpos<-NULL; ypos<-NULL; nchild<-NULL
  ctmp <- unique(partition[parent==inip])
  cindex<- index[is.element(upar, ctmp)]
  nchild<-length(ctmp)
  eps<-.15
  xpos <- rep(0,nchild)+c(1:nchild)*inter[1]
  xpos<- xpos-mean(xpos)+sx
  xx <- rep(sx,2*nchild)
  xx[seq(2,2*nchild,2)]<-xpos
  yy<-rep(c(sy,sy-yint),nchild)
  lines(xx,yy, col=4, lty=2, lw=line1)
  sw(xx,yy); cin1<-NULL; ctm1<-NULL   
  for(i in 1:nchild){
    xp <-xpos[i]+ (sx-xpos[i])*.2
    yp <- sy-yint+yint*.5
    text( xp, yp, cindex[i], cex=num1, col=2)
    cin<- c(cin, cindex[i])
    ctm<- c(ctm, ctmp[i])
    ccandi <- unique(partition[parent==ctmp[i]]) 
    ok<-T; tmp<-NULL; lc<-1
    while(ok){
      tmp<-partition[is.element(parent,ccandi)]
      ccandi<-unique(c(ccandi,tmp))
      nlc<-length(ccandi)
      if((nlc-lc)>0) lc<-nlc
      else ok<-F
    }
    kk <- is.element(partition, ccandi)
    if(sum(kk)>0){
      a<-plot.stree(cin, ctm,upar,index, ctmp[i], partition[kk],parent[kk],
                    names[kk], xpos[i],yy[2],yint,inter[-1],line1,num1,letter1)
      cin<-c(cin, a$cin)
      ctm<-c(ctm, a$ctm)
    }
  }
  result<-list(cin=cin, ctm=ctm)
  return(result)
}   

#######################################################################
# Function iplot.tree used inside the function tanal and draws the 
# second figure.
#######################################################################
iplot.tree<-function(csave, inter,data, partition,parent, names, main,line1,num1){
  upar<-unique(partition)
  n<-length(partition)
  leaf<-setdiff(partition, parent)
  nleaf<-length(leaf)
  layer<-1;
  ok <- !(nleaf==length(upar))
  xpos<-NULL; ypos<-NULL; nchild<-NULL
  ccandi <- unique(partition[parent==0])
  ctmp<-ccandi; 
  nchild<-length(ctmp);
  while(ok){
    layer<-layer+1
    tmp<-partition[is.element(parent,ccandi)]
    ttmp<-setdiff(tmp, ccandi)
    ccandi<-unique(c(ccandi,tmp))
    if(sum(!is.element(ttmp, leaf))<=0) ok<-F;
  }
  xmax<-prod(csave[c(1:layer)]); ymax<-max(10, length(upar)+1)
  plot(1,1, type="n",xlab="", ylab="", axes=F,xlim=c(0,xmax),ylim=c(-.6,ymax))
  if(!is.null(main)) mtext(main, side=3, line=1, cex=1.2)

  yint <- ymax/layer;  
  if(nchild>1){ 
    xpos <- rep(0,nchild)+c(1:nchild)*inter[1]
    xpos<- xpos-mean(xpos)+xmax/2
  }
  else xpos<-xmax/2
  xx <- rep(xmax/2,2*nchild)
  xx[seq(2,2*nchild,2)]<-xpos
  yy<-rep(c(ymax,ymax-yint),nchild)
  lines(xx,yy, lty=2,lwd=line1,col=4)
  sw(xx,yy)
  eps<- .2; 
  for(i in 1:nchild){ 
    pdata<-as.matrix(data[,partition==ctmp[i]])
    lpcandi<-sum(apply(pdata,1,prod))
    xp <-xpos[i]+ (xmax/2-xpos[i])*.6 
    yp <- ymax-yint+yint*.6 
    text( xp, yp, lpcandi, cex=num1, col=3)
    ccandi <- unique(partition[parent==ctmp[i]])    
    ok<-T; tmp<-NULL; lc<-1
    while(ok){
      tmp<-partition[is.element(parent,ccandi)]
      ccandi<-unique(c(ccandi,tmp))
      nlc<-length(ccandi)
      if((nlc-lc)>0) lc<-nlc
      else ok<-F
    }
    kk <- is.element(partition, ccandi)
    if(sum(kk)>0)
    iplot.stree(pdata, data[,kk],ctmp[i], partition[kk],parent[kk],names[kk],
                xpos[i],yy[2], yint,inter[-1],line1,num1)
    else text(xpos[i],ymax-yint-.5,lpcandi,cex=num1,col=2)
  }   
}   
 
  
########################################################################
# Function iplot.stree used inside the function tanal and iplot.tree
########################################################################
iplot.stree<-function(pdata, data,inip, partition, parent, names, sx, sy,
                      yint, inter,line1,num1){
  leaf<-setdiff(partition, parent)
  xpos<-NULL; ypos<-NULL; nchild<-NULL
  ctmp <- unique(partition[parent==inip])
  #cindex<- index[is.element(upar, ctmp)]
  nchild<-length(ctmp)
  eps<-.15
  xpos <- rep(0,nchild)+c(1:nchild)*inter[1]
  xpos<- xpos-mean(xpos)+sx
  xx <- rep(sx,2*nchild)
  xx[seq(2,2*nchild,2)]<-xpos
  yy<-rep(c(sy,sy-yint),nchild)
  lines(xx,yy, col=4, lty=2, lw=line1)
  sw(xx,yy); 
  for(i in 1:nchild){
    ppdata<- as.matrix(data[,partition==ctmp[i]])
    tdata<- as.matrix(apply(cbind(ppdata, pdata), 1, prod))
    lpcandi<-sum(apply(ppdata,1,prod))
    xp <-xpos[i]+ (sx-xpos[i])*.2 
    yp <- sy-yint+yint*.5
    text( xp, yp, lpcandi, cex=num1,col=3)
    ccandi <- unique(partition[parent==ctmp[i]]) 
    ok<-T; tmp<-NULL; lc<-1
    while(ok){
      tmp<-partition[is.element(parent,ccandi)]
      ccandi<-unique(c(ccandi,tmp))
      nlc<-length(ccandi)
      if((nlc-lc)>0) lc<-nlc
      else ok<-F
    }
    kk <- is.element(partition, ccandi)
    if(sum(kk)>0)
       iplot.stree(tdata,data[,kk],ctmp[i], partition[kk],parent[kk],names[kk],
                   xpos[i],yy[2], yint, inter[-1],line1,num1)
    else text( xpos[i], yy[2]-.5, sum(tdata), cex=num1, col=2)
    }
}   


######################################################################
# Function ematrix returns the ensemble of a tree and the number of 
# tumors whose all aberrations on an ensemble are non-neutral.
####################################################################### 
ematrix<-function(partition, neutral, parent, data){
  index1<-c(1:length(partition))
  index<-dimnames(data)[[2]]
  np<-parent; npa<-partition 
  if(any(neutral)){
    leaf<-setdiff(partition, parent); tl<-length(partition)
    apart<-partition[!neutral]
    incidence <- outer( apart, unique(partition), "==" )
    edgelengths <- c( rep(1,sum(!neutral)) %*% incidence )
    upartition<-unique(partition)
    lc<-upartition[edgelengths==0]
    lc1<-setdiff(lc, leaf)
    if(length(lc1)>0){
      for(ii in 1:length(lc1)){
        np[is.element(np, lc1[ii])]<-unique(np[partition==lc1[ii]])
      }
    }
    lc2<-intersect(lc, leaf)
    if(length(lc2)>0){
      for(ii in 1:length(lc2)){
        kk<-npa[(np==np[npa==lc2[ii]])*(!neutral)==1]
        kk<-setdiff(kk,lc2[ii])
        if(length(kk)==1){
          if(unique(np[partition==kk])!=0){
            upp<-unique(np[npa==kk])
            if(all(neutral[np==upp])){
              ok<-T
              while(ok){
                upp<-unique(parent[partition==upp])
                ok<-all(neutral[partition==upp])
              }
            }
            npa[partition==kk]<- upp
            np[partition==kk]<-unique(np[partition==upp])
            if(!is.element(kk, leaf)) np[np==kk]<-upp
         }
       }
     }
    }
    parent<-np[!neutral]; partition<-npa[!neutral]; index<-index[!neutral]
    index1<-index1[!neutral]
  }
  ################# ensemble
  if(!all(parent==0)){
    leaf<-setdiff(partition, parent); lleaf<-length(leaf)
    ensm<-list(); ensm1<-list()
    for(i in 1:lleaf){
      a<-leaf[i]; b<-a; ok1<-b!=0 
      while(ok1){
        b<-unique(parent[partition==b])
        ok1<- b!=0
        a<-c(a,b)
      }
      ens <- a[a!=0]
      ensm[[i]]<-index[is.element(partition, ens)]
      ensm1[[i]]<-index1[is.element(partition, ens)]
    }
    lensm<-length(ensm)
    nensm<-matrix(0,nrow(data), lensm)
    for(ii in 1:nrow(data)){
      for(j in 1:lensm){
        nensm[ii,j] <- all(data[ii, as.numeric(ensm1[[j]])])
      }
    }
    dimnames(nensm)<-list(dimnames(data)[[1]], ensm)
  }
  else{
    edge <- list(); upar <- unique(partition); edge1<-list()
    lupar <- length(upar); npa<-nrow(data)
    nensm<-matrix(0, npa, lupar)
    for(i in 1:lupar){
      edge[[i]]<-index[is.element(partition, upar[i])]
      edge1[[i]]<-index1[is.element(partition, upar[i])]
      for(ii in 1:npa){
        nensm[ii,i] <- all(data[ii, as.numeric(edge1[[i]])])
      }
    }
    dimnames(nensm)<-list(dimnames(data)[[1]], edge)
  }
  result <- list(nensm=nensm, tree=apply(nensm,2,sum))
  return(result)
}

########################################################################
# Function edmatrix returns the edges of a tree and the number of
# tumors whose all aberrations on a edge are non-neutral.
########################################################################
edmatrix<-function(partition, neutral, data){
  index1<-c(1:length(partition))
  index<-dimnames(data)[[2]]
  index<-index[!neutral]; partition<-partition[!neutral]
  index1<-index1[!neutral]
  edge <- list(); upar <- unique(partition); edge1<-list()
  lupar <- length(upar); npa<-nrow(data)
  nedge<-matrix(0, npa, lupar)
  for(i in 1:lupar){
    edge[[i]]<-index[is.element(partition, upar[i])]
    edge1[[i]]<-index1[is.element(partition, upar[i])]
    for(ii in 1:npa){
      nedge[ii,i] <- all(data[ii, as.numeric(edge1[[i]])])
    }
  }
  dimnames(nedge)<-list(dimnames(data)[[1]], edge)
  tedge<-apply(nedge, 2, sum)
  return(tedge)
}

