try2<-function(n){
  # List all possible result of n Bernoulli try. Rerutn is n by 2^n-1 matrix.
  re<-NULL
  for(i in 1:n){
    re<-rbind(re,try1(n,i,NULL))
  }
  return(t(re))
}

try1<-function(n,l, re){
  # Recursive function used inside of try2. 
  if(n==l){
    re<-matrix(rep(1,n),1,n)
  }
  else{
    if(l==1){
      re<-matrix(0,n,n)
      diag(re)<-1
    }
    else{
      for(i in 1:(n-l+1)){
        aa<-try1(n-i,l-1, NULL)
        bb<-matrix(c(rep(0,(i-1)*nrow(aa)),rep(1,nrow(aa))),nrow(aa),i)
        re<-rbind(re,cbind(bb,aa))
      }
    }
  }
  return(re)
}

Ginter <- function(data,k1,k2,k3,k4){
  #########################################################################
  # Generate single, pair, triple and quadruple interactions experienced by
  # more than k1 (single), k2(pair), k3(triple) and k4(quadruple) tumors.
  #########################################################################
  ntum<-nrow(data)
  nab<-ncol(data)
  ensem<-NULL; 
  for(j in 1:nab){
    if(sum(data[,j]==1)>=k1){
      tmp<-matrix(0,1,nab); tmp[,j]<-1; ensem<-rbind(ensem, tmp)
    }
  }
  for(i in 1:(nab-1)){
    for(j in (i+1):nab){
      if(sum(data[,i]*data[,j])>=k2){
        tmp<-matrix(0,1,nab); tmp[,c(i,j)]<-1; ensem<-rbind(ensem, tmp)
      } 
    }
  }
  for(i in 1:(nab-2)){
    for(j in (i+1):(nab-1)){
      for(g in (j+1):nab){
        if(sum(data[,i]*data[,j]*data[,g])>=k3){
          tmp<-matrix(0,1,nab); tmp[,c(i,j,g)]<-1; ensem<-rbind(ensem, tmp)
        }         
      }
    }
  }
  for(i in 1:(nab-3)){
    for(j in (i+1):(nab-2)){
      for(g in (j+1):(nab-1)){
        for(h in (g+1):nab){
          if(sum(data[,i]*data[,j]*data[,g]*data[,h])>=k4){
            tmp<-matrix(0,1,nab); tmp[,c(i,j,g,h)]<-1; ensem<-rbind(ensem, tmp)
          }
        }         
      }
    }
  }  
  return(ensem)
}

# Plot the number of interactions vs. number of aberrations.
Ginternum <- function(data){
  nab<-ncol(data);   ntum<-nrow(data); 
  re1<-rep(0,ntum); re2<-rep(0,ntum); re3<-rep(0,ntum);re4<-rep(0,ntum); 
  for(i in 1:nab){
    t1 <- sum(data[,i])
    re1[t1+1] <- re1[t1+1]+1 
  }
  re1 <- nab - cumsum(re1)

  for(i in 1:(nab-1)){
    for(j in (i+1):nab){
      t2 <- sum(data[,i]*data[,j])
      re2[t2+1]<-re2[t2+1]+1 
    }
  }
  re2 <-  nab*(nab-1)/2 - cumsum(re2)
   
  for(i in 1:(nab-2)){
    for(j in (i+1):(nab-1)){
      for(g in (j+1):nab){
        t3<- sum(data[,i]*data[,j]*data[,g])
        re3[t3+1]<-re3[t3+1]+1 
      }
    }
  }
  re3 <-  nab*(nab-1)*(nab-2)/6 - cumsum(re3)

  for(i in 1:(nab-3)){
    for(j in (i+1):(nab-2)){
      for(g in (j+1):(nab-1)){
        for(h in (g+1):nab){
          t4 <- sum(data[,i]*data[,j]*data[,g]*data[,h])
          re4[t4+1]<-re4[t4+1]+1        
        }
      }
    }
  }
  re4 <-  nab*(nab-1)*(nab-2)*(nab-3)/24 - cumsum(re4)
  
  plot(c(1:ntum),re1,log="y",ylim=c(1,re4[1]), type="s",
       xlab="Number of tumors", ylab="Number of interactions")
  legend(re1[1],10, c("Single", "Pair","Triple", "Quadruple"), bty="n",cex=.8,
         lty = c(1:4))
  lines(c(1:ntum), re2, lty=2, type="s")
  lines(c(1:ntum), re3, lty=3, type="s")
  lines(c(1:ntum), re4, lty=4, type="s")  
}

Grisn <- function(params,C,numtumor){
  ########################################################################
  # Sample random binary vectors (numtumor tumors) from general netowrk, C,
  #  using rejection sampling.
  #
  # params = (alpha, p.a, p.b)  
  #       alpha = Prob[ X=1 ] : estimated from not significant aberrations
  #       p.a = Prob[ Aberration | X=1 ] : error terms; .99 or .98...
  #       p.b = Prob[ Aberration | X=0 ] : error terms; .1 or .2...
  #   C(General network): kxn matrix; k ensembles, n aberrations
  #   C[i,j]=1 if jth aberration resides on ith ensemble
  # numtumor : number of tumors to sample. 
  ########################################################################
  p.a <- params[2]
  p.b <- params[3]
  alpha <- params[1]
  theta <- p.b*(1-alpha)+p.a*alpha

  C <- t(C)
  k <- ncol(C)
  n <- nrow(C)
  # B is lists all the ways of combining columns for JK's alg.
  # strongly recommend to store the results for 1 <= k < 20 and use them. 
  B <- try2(k) 
  if(k==1) B <- t(as.matrix(B))
  n.one <- c( rep(1,k) %*% B ) 
  Bstar <- 1*( C%*%B > 0 )  
  lengths <- c( rep(1,n)  %*% Bstar )
  p.sel <- sum(p.a^lengths[n.one==1] )   # selection probability

  xsave <- matrix(NA,numtumor,n) # save the random vectors
  for( i in 1:numtumor ){ 
    notdone <- T 
    while( notdone ){ 
      xi <- ifelse( runif(n) <= alpha, 1, 0 )
      xi <- t(as.matrix(xi)) 
      # for each ensemble combination indicate which aberrations are included
      stat.a <- xi %*% Bstar
      stat.b <- (1-xi) %*% Bstar   
      in.ex <- (p.a^stat.a)*(p.b^stat.b)
      # inclusion exclusion probabilities, ncol(B) x ntumor
      p.sel.x <- in.ex %*% ((-1)^(n.one+1))  # added up
      # put together p(sel|x) p(x)/p(sel)
      env <- sum( log(p.sel.x)- log(p.sel) ) 
      notdone <- ifelse( log(runif(1)) <= env, F, T ) 
    } 
    xsave[i,] <- xi 
  } 
  return( xsave )
 }

Gcandi<-function(data,pval1,pval2,error){
  #########################################################################
  # Make a candidat set L = L1 or L2
  # the higest order of a candidate is 10.
  # L1 : 
  # pval1 = probability of upper tail in binomial distribution (to select L1)
  # pval2 = percentage of tumor number not included in L2 
  # error = error term for positive internal & negative external covariance
  #########################################################################
  dat<-data; nab<-ncol(dat); ntum<-nrow(dat)
  i<-c(1:10) # the higest order of a candidate is 10. 
  p1<-sum(dat)/length(dat)
  c1 <- qbinom(1-pval1,ntum,p1^i) # threshold for L1
  c2 <- ntum*pval2 # threshold for L2
  c1[c1<c2]<-c2    # threshold L1 should bigger than that of L2
  
  set1 <- NULL; set2 <- NULL; tem2 <- rep(T, nab)
  ok <- T; l<-1; nabt<-nab
  while(ok){
    tset1 <- NULL; tset2 <- NULL;  
    a <- try1(nab, l, NULL);
    tmp <- apply(dat %*% t(a)==l, 2, sum)
    fset <- tmp >= c1[l]; sset <- (tmp >= c2 & tmp < c1[l])
    sfset <- sum(fset); ssset <- sum(sset)
    # select L1
    if( sfset > 0 ){
      if( l==1 ) tset1 <- a[fset,]
      if( l > 1 & sfset==1){
        if(  min( cov(dat[,a[fset,]==1]) ) >= -error ){
          tset1 <- a[fset,]
        }
      }
      if( l > 1 & sfset > 1){
        aa <- a[fset,]
        for(j in 1:nrow(aa)){
          if(  min( cov(dat[,aa[j,]==1]) ) >= -error ){
            tset1 <- rbind(tset1,aa[j,])
          }
        }
      }
      if(!is.null(tset1) & any(!tem2) ){
        tiset1 <- matrix(0,nrow(tset1),nabt)
        tiset1[,tem2]<- tset1
      }
      else tiset1 <- tset1
      set1 <- rbind(set1, tiset1)
    }
    # select L2
    if( ssset > 0){
      if( l==1 ) tset2 <- a[sset,]
      if( l > 1 & ssset==1){
        if(  min( cov(dat[,a[sset,]==1]) ) >= -error ){
          tset2 <- a[sset,]
        }
      }
      if( l > 1 & ssset > 1){
        aa <- a[sset,]
        for(j in 1:nrow(aa)){
          if(  min( cov(dat[,aa[j,]==1]) ) >= -error ){
            tset2 <- rbind(tset2,aa[j,])
          }
        }
      }
      if(!is.null(tset2) & any(!tem2) ){
        tiset2 <- matrix(0,nrow(tset2),nabt)
        tiset2[,tem2]<- tset2
      }
      else tiset2 <- tset2
      set2 <- rbind(set2, tiset2)
    }
    # reduce data if possible to save computational time
    if((!is.null(set1)||!is.null(set2))){
      if(l>3||is.null(tset1)||is.null(tset2)){
        tem1 <- rbind(set1,set2) 
        tem2 <- apply(tem1, 2, sum)>0
        nab<-sum(tem2); dat<-data[,tem2]
        if(nab==1) dat <- as.matrix(dat)
      }
      l<-l+1; ok <- !is.null(tset1) || !is.null(tset2)
    }
  }
  # L2 has externally negative covariance with at least one candidate in L1  
  set <- NULL
  if((!is.null(set1)&!is.null(set2))){
    fnum<-nrow(set1); snum<-nrow(set2); tmp<-rep(F,snum)
    for(sj in 1:snum){
      ok<-T; si<-1
      while(ok & si <= fnum){
        if(max(cov(data[,set1[si,]==1],data[,set2[sj,]==1])) <= error){
          ok<-F; tmp[sj]<-T; 
        }
        si<- si+1; 
      }
    } 
    set2<-set2[tmp==1,]; 
    set <-rbind(set1,set2)
  }
  else if((is.null(set1)||!is.null(set2))) set <- set2
  else if((!is.null(set1)||is.null(set2))) set <- set1
  return(set)
}

mpara<-function(alpha){
  ########################################################################
  # Given alpha (Prob[ X=1 ]), return possible parameter settings
  # (theta,delta, gamm), when P(z=0|x=1)=1-p.a, P(z=1|x=0)=p.b changes
  # .1 to .4
  ########################################################################
  res<-NULL
  pa<-seq(.99 ,.96, -.01)
  pb<-seq(.01,.04,.01)
  for(i in 1:length(pa)){
    for(j in 1:length(pb)){
      pa1<-pa[i]; pb1<-pb[j]
      theta<- pb1*(1-alpha)+pa1*alpha
      delta <- pb1*(1-alpha)/theta
      gamma <- (alpha-theta*(1-delta))/(1-theta)
      res<-rbind(res, c(pa1,pb1,theta,delta, gamma))
    }
  }
  colnames(res) <- c("p.a","p.b","theta","delta", "gamma")
  return(res)
}

Glik <- function(C,dat,params){
  ########################################################################
  # Calculate the log joint probability density of 
  # aberrations in the ISGN model. Implementation was suggested by Jay Kadane.
  #
  # C (General network): kxn matrix; k ensembles, n aberrations.
  #   C[i,j]=1 if jth aberration resides on ith ensemble
  # dat : m x n binary matrix, data from m patients and observe n aberrations
  # from each tumor.
  # params = (p.a, p.b, alpha, theta) IS parameter values
  #       alpha = Prob[ X=1 ] : estimated from not significant aberrations
  #       p.a = Prob[ Aberration | X=1 ] : error terms; .99 or .98...
  #       p.b = Prob[ Aberration | X=0 ] : error terms; .1 or .2...
  #       theta = Prob[ Z=1 ] 
  ########################################################################
  #alpha <- params[1]; delta <- params[2]; gamma <- params[3]
  #theta<- (alpha - gamma)/(1-delta-gamma)  
  #p.b <- theta*delta/(1-alpha) # pb
  #p.a <- theta*(1-delta)/alpha # pa

  p.a <- params[1]
  p.b <- params[2]
  alpha <- params[3]
  theta <- p.b*(1-alpha)+p.a*alpha
  #theta <- params[4] 
  C<-t(C)
  k <- ncol(C)
  n <- nrow(C)
  # B is lists all the ways of combining columns for JK's alg.
  # strongly recommend to store the results for 1 <= k < 20 and use them. 
  B<-try2(k) 
  if(k==1) B<-t(as.matrix(B))
  n.one <- c( rep(1,k) %*% B ) 

  # for each ensemble combination indicate which aberrations are included
  Bstar <- 1*( C%*%B > 0 )  
  lengths <- c( rep(1,n)  %*% Bstar )
  p.sel <- sum(  (-1)^(n.one+1) *theta^lengths )   # selection probability
  # statistics for every ensemble combination for every tumor
  stat.a <- dat %*% Bstar
  stat.b <- (1-dat) %*% Bstar   
  in.ex <- (p.a^stat.a)*(p.b^stat.b)
  # inclusion exclusion probabilities, ncol(B) x ntumor
  p.sel.x <- in.ex %*% ((-1)^(n.one+1))  # added up
  # put together p(sel|x) p(x)/p(sel)
  loglik <- sum( log(p.sel.x) - log(p.sel) ) +
      sum(dat)*log(alpha) + sum(1-dat)*log(1-alpha)
  #tloglik<-list(top=sum(log(p.sel.x)), bot=116*log(p.sel), loglik=loglik)
  return( loglik )
 }


check<-function(network, candi){
  #########################################################################
  # check if a candidate is a subset of any current ensemble
  # if not (if return is T), the candidate can be a part of network.
  #
  # network : binary matrix, (number of ensemble) by (number of tumor)
  #           each row represents an ensemble
  # candi : binary vector encodes a candidate
  #########################################################################
  tmp <- candi %*% t(network)
  anetwork <- apply(network,1, sum); scandi<-sum(candi)
  tmp1<-(anetwork - scandi)>=0
  okg <- all(tmp[tmp1]<scandi) & all(tmp[!tmp1]-anetwork[!tmp1]<0)
  return(okg)
}


checkm <- function(network){
  #########################################################################
  # check if network is legal; if each ensemble (row) is not subset of
  # the other row then it returns T.
  #
  # network : binary matrix, (number of ensemble) by (number of tumor)
  #           each row represents an ensemble
  #########################################################################
  tmp <- network %*% t(network)
  diag(tmp)<- 0
  okg <- all(apply(network, 1, sum)-apply(tmp,1,max)>0) 
  return(okg)
}

Gnet <- function(candid, tmaxe){
  #######################################################################
  # Simulate a network (used to start the ISGN).
  # From candidate sets construct a network having less or equal to tmaxe
  # ensembles.
  #
  # candid : candidate set
  # tmax : number of maximum ensemble in a network
  #######################################################################
  ncandi<- nrow(candid); gind <- rep(F, ncandi); ok<-T; index<-c(1:ncandi)
  while(ok){ # construct a network till it satisfies no subset condition
    size <- sample(c(1:tmaxe),1)
    tmp<-sample(index, size)
    network <- candid[tmp,]
    if(size==1) network<-t(as.matrix(network))
    if(checkm(network)) ok<-F
  }
  gind<-rep(F,ncandi); gind[tmp]<-T
  re<-list(network = network, gindex = gind)  
  return(re)
}


ISGN <- function(data, mcmc, para, pri, candid){
  ###############################################################
  # main MCMC function
  #
  # data =  m x n binary matrix, data from m patients and observe n aberrations
  # mcmc = list with driving instructions
  #        padI   = roportion of switches to change in last move type
  #        nskip  = size of stretch between saved states
  #        nsave  = number of saved states
  #        tmax   = number of maximum candidble in a network
  # para = (alpha, p.a, p.b) 
  #     alpha = Prob[ X=1 ] : estimated from not significant aberrations
  #       p.a = Prob[ Aberration | X=1 ] : error terms; .99 or .98...
  #       p.b = Prob[ Aberration | X=0 ] : error terms; .1 or .2...
  #  Note that the likelihood function is defined using (theta, delta, gamma)
  #  but the input parameters is alpha, p.a and p.b. where, 
  #     theta = Prob[ Z=1 ] : estimated from not significant aberrations
  #     delta = Prob[ X=1 | no Aberration ] 
  #     gamma = Prob[ X=0 | Aberration ]
  # pri =  prior, either uniform or uniform over the number of ensemble.
  # candid = candidate set
  ###############################################################

  padI  <- mcmc$padI
  nscan <- mcmc$nsave*mcmc$nskip
  nsave <- mcmc$nsave
  nskip <- mcmc$nskip
  tmax  <- mcmc$tmax

  alpha <- para[1]; p.a <- para[2]; p.b <- para[3]
  theta <- p.b*(1-alpha)+p.a*alpha 
  param<-c(p.a, p.b, alpha, theta)

  ntum <- nrow(candid); nab <- ncol(candid);
  index<-c(1:ntum); 
  abnum <- matrix(1, nab,1); acandid <- apply(candid,1,sum)

  gsave <- matrix(F, nsave, ntum);
  # saved network; each row represents one network; 
  liksave <- numeric(nsave);      # saved log-likelihood
  prisave <- numeric(nsave);      # saved prior
  acrate <- rep(0,3); 
  a <- Gnet(candid,tmax);  # generate one network

  network <- a$network;
  gindex<-rep(0,ntum); gindex[a$gindex]<-T  
  # Evaluate the loglikelihood and logprior of starting state
  loglik <- Glik(network, data, param)
  logpri <- pri[sum(gindex)]


  ####################################################################
  # Run MCMC.
  # Propose a new network configuration (3 move types)
  ####################################################################
  skipcount <- 0
  isave <- 1

  for(iscan in 1:nscan){
    skipcount <- skipcount+1;
    ###################################################################
    # add one or remove one candidate
    lgindex<-sum(gindex); logmh<-NA; lgindexs<-0;gindexs<-gindex
    candi<-sample(index, 1)
    if(gindex[candi]){ # REMOVE
      gindexs[candi]<-F
      lgindexs<- lgindex -1
    }
    else{
      gcandi<-candid[candi,]
      if(check(network, gcandi)){ # ADD 
        gindexs[candi]<-T       
        lgindexs<- lgindex +1
      }
    }
    if(lgindexs<=tmax & lgindexs>0){
      networks<- candid[gindexs==1,]
      if(lgindexs==1) networks<-t(as.matrix(networks))
      ploglik <- Glik(networks, data, param)
      plogpri <-pri[lgindexs]
      logmh <- ploglik-loglik + plogpri-logpri
      ok <- !is.na(logmh)
      if( ok ){
        accept <- ( log(runif(1)) <= logmh )
        if( accept ){
          gindex <- gindexs
          network <- networks
          logpri <- plogpri
          loglik <- ploglik
          acrate[1] <- acrate[1]+1
        }
      }
    }

    ###################################################################
    # swap two candidates
    if( any(gindex==1) & any(gindex==0)){
      lgindex <- sum(gindex); logmh <- NA;
      lgindexs <- 0; gindexs<-gindex
      # choose one from the network and remove it
      if(lgindex==1) candiT <- index[gindex==1]
      else candiT<-sample(index[gindex==1],1)
      # choose one from the candidate set and add it
      if(lgindex+1==ntum) candiF <- index[gindex==0]
      else candiF<-sample(index[gindex==0],1)

      gindexs[candiT]<- F; lnetworks<- candid[candiF,]
      networks<-candid[gindexs==1,]
      if(lgindex==2) networks<-t(as.matrix(networks))
      if(check(networks, lnetworks)){
        gindexs[candiF]<-T
        networks<- candid[gindexs==1,]
        if(lgindex==1) networks<-t(as.matrix(networks))
        ploglik <- Glik(networks, data, param)
        logmh <- ploglik-loglik 
        ok <- !is.na(logmh)
        if( ok ){
          accept <- ( log(runif(1)) <= logmh )
          if( accept ){
            gindex<-gindexs
            network <- networks
            loglik <- ploglik
            acrate[2] <- acrate[2]+1
          }
        }
      }
    }
    #################################################################
    # swap sets of candidate
    lgindex <- sum(gindex); logmh <- NA; lgindexs <- 0; gindexs <- gindex
    change <- ( runif(ntum) <= padI ); okc<-F
    if( any(change) ){
      gindexs[change] <- !gindex[change]
      networks <- candid[gindexs==1,]
      lgindexs <- sum(gindexs)
      if(lgindexs==1 & lgindexs>0){
        okc<-T; networks <- t(as.matrix(networks))
      }
      else
        if(lgindexs>0 & lgindexs<=tmax) okc <- checkm(networks)
      if(okc){
        plogpri <-pri[lgindexs]
        ploglik <- Glik(networks, data, param)
        logmh <- ploglik-loglik + plogpri-logpri
        ok <- !is.na(logmh)
        if( ok ){
          accept <- ( log(runif(1)) <= logmh )
          if( accept ){
            gindex<-gindexs
            network <- networks
            logpri <- plogpri
            loglik <- ploglik
            acrate[3] <- acrate[3]+1
          }
        }
      }
    }

    #######################################################################
    # Store summary statistics periodically..
    if(skipcount == nskip ){
      #if(is.element(isave, testp)) print(gindex)
      skipcount <- 0
      gsave[isave,]<-gindex
      liksave[isave] <- loglik
      prisave[isave]<-logpri
      isave <- isave+1
    }
  }
  post <- liksave+prisave; index<-c(1:nsave)
  netsave <- index[max(post)==post][1]
  netsave <- gsave[netsave,]
  snet <- sum(netsave==1)
  bestnet <-candid[netsave==1,]
  if(snet==1) bestnet <- t(as.matrix(bestnet))  

  mnetwork <- gnetwork(colnames(data), data, candid, bestnet, snet)
  result<-list(data=data, candid=candid, gsave=gsave, logpri=prisave,
    acrate=acrate,loglik=liksave, bestnet=bestnet,mnetwork=mnetwork)
  return( result )
}

gnetwork <- function(cname, data, candid, bestnet,snet){
  k<- list(); num <- NULL
  for(i in 1:snet){
    k[[i]] <- cname[bestnet[i,]==1]
    saa <- sum(bestnet[i,]==1)
    if(saa > 1) num <- c(num, sum(apply(data[,bestnet[i,]==1],1,min)))
    else num <- c(num, sum(data[,bestnet[i,]==1]))
  }
  names(num) <- k
  return(num)
}

Gscore <- function(data, network){
  ########################################################################
  # Scoring a network based on how well it can explane the observed data
  # It returns 4 different scores :
  #   i) score based on match/mismatch the data
  #  ii) adjuste score using marginal aberration rate
  # iii) score based on covariance
  #  iv) score based on distance - absolute distance
  #
  # data =  m x n binary matrix, data from m patients and observe n aberrations
  # network = general network
  ########################################################################
  
  te<-matrix(0,1,ncol(data)) # empty network
  tf<-matrix(1,1,ncol(data)) # full network
  tdata <- t(data); 
  # i) score based on match and mismatch 
  k <- network%*%tdata+(1-network)%*%(1-tdata)- # match
    (1-network)%*%tdata- network%*%(1-tdata)    # mismatch
  # ii) score adjust the marginal aberration rate
  be<-(1-te)%*%(1-tdata); bbe<-mean(be); bf<-tf%*%tdata; bbf<-mean(bf)
  k1 <- network%*%tdata/bbf+(1-network)%*%(1-tdata)/bbe-
       (1-network)%*%tdata/bbf-network%*%(1-tdata)/bbe
  # maximum score among scores from each ensemble
  if(nrow(network)==1){kk<-k; kk1<-k1}
  else{ kk<-apply(k,2, max); kk1<-apply(k1,2, max)}
  sc <- matrix(0, nrow(data), nrow(network)) # score based on covariance
  sd <- matrix(0, nrow(data), nrow(network)) # score based on distance
  for(i in 1:nrow(data)){
    for(j in 1:nrow(network)){
      sc[i,j]<-cov(data[i,], network[j,])
      sd[i,j]<-dist(rbind(data[i,], network[j,]), method="manhattan")
     }
   } 
  cc<-apply(sc,1,max); dd<-apply(sc,1,min)
  a<-c(mean(kk)/ncol(data), mean(kk1), mean(sc), mean(sd))
  names(a)<-c("score","adj.s","cov.s","dist.s")
  return(round(a,3))
}


