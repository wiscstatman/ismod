library(ISN)

## test of ISGN functions.
## generate 50 data from general network C
para <- c( .08,.99,.01); numtumor <- 50
C <- matrix(0,4,20); 
C[1,]<-c(1,1,0,0,1,0,1,rep(0,13)); C[2,]<-c(0,0,0,0,1,1,1,rep(0,13)); 
C[3,]<-c(0,0,0,1,0,0,0,rep(0,13)); C[4,]<-c(1,1,1,0,0,0,0,rep(0,13))
data <- Grisn(para,C,numtumor)
colnames(data)<- letters[1:ncol(C)]
## generate a candidate set
pval1 <- 10^(-10); pval2<- .1; error <- .005
candid <- Gcandi(data,pval1,pval2,error)

## run MCMC
mcmc<-list(padI=.01,nsave=1000,nskip=1, tmax=4)
pri<-rep(0,25)
a<-ISGN(data, mcmc, para, pri, candid)

## check best network
a$mnetwork
## calculate the score
Gscore(data, a$bestnet)

#source('ISN/R/ISTN.R')

## test of ISTN functions
## generate 100 data from a tree-like network
partition <- c(1:3,0,5,6)
parent<-c(0,1,1,-1,2,2)
para<-c(0.05, 0.01, 0.01); numtumor<-100
data <- Trisn(para, partition, parent, numtumor)
colnames(data)<-letters[1:length(partition)]

## run MCMC
mcmc <- list(nsave=10, nskip=100, nperm=3, padI=.15, padII=.15)
res<-res(6)
mh <- ISTN(data,mcmc=mcmc,tau=4,para,res)

# Hierachical clustering
Tcluster(mh, colnames(data))

# check the best tree-like network and draw it. 
mh$mtree
# draw a tree-like network
post <- mh$loglik; index<-c(1:mcmc$nsave)
netsave <- index[max(post)==post][1]
Tdraw(mh$partition[netsave,],mh$neutral[netsave,],mh$parent[netsave,],data)

