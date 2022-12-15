\name{Tdraw}
\alias{Tdraw}
\title{plot the tree-like network}
\description{
  Draw a tree-like network (top) and the number of tumors which 
  experiance (all aberrations reside on that edge are active)
  each edge (bottom).
}
\usage{
 Tdraw(partition,neutral,parent,data,line1,num1,letter1,
                  main1=NULL, main2=NULL)
}
\arguments{
  \item{partition}{ncomponent vector holding distinct edge labels}
  \item{neutral}{logical vector indicating neutral components}
  \item{parent}{like partition, but holds labels for parent edges
    0 is code for root node}
  \item{data}{ntumor x nswitch matrix of binaries 0/1 indicating
    abnormality length should be positive}
  \item{line1}{line width}
  \item{num1}{number of tumor size}
  \item{letter1}{aberration text size}
  \item{main}{title for the top and bottom figure}
}
\examples{
index<-c(1:mcmc$nsave) # mcmc is result from ISTN.R

i<-index[max(mh$loglik)==mh$loglik][1]

Tdraw(mh$partition[i,],mh$neutral[i,],mh$parent[i,],sdata,
line1=1,num1=1,letter1=1,main1="MLE tree-like network",
main2="Number of tumors experiencing each edge/ensemble")  
}
\value{
  plot tree-like network.
}
\keyword{plots}
