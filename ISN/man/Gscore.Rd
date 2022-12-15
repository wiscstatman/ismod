\name{Gscore}
\alias{Gscore}
\title{scoring a network}
\description{
  Scoring a network based on how well it can explane the observed data 
}
\usage{
Gscore(data, network)
}
\arguments{
  \item{data}{cCGH data profile}
  \item{newtork}{number of ensemble by number of aberration matrix}
}
\value{return 4 different scores :

  i) score based on match/mismatch the data
  
  ii) adjuste score using marginal aberration rate
  
  iii) score based on covariance
  
  iv) score based on distance - absolute distance
}
