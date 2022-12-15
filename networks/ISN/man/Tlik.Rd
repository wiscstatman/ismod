\name{Tlik}
\alias{Tlik}
\title{Calculate the log joint probability density}
\description{
  Calculate the log joint probability density of aberrations in the 
  Instability Selection Tree-like Network model.
}
\usage{
Glik(rates,partition,parent,neutral,data)
}
\arguments{
  \item{rates}{(theta,delta, gamma), here theta = Prob[ Aberration ],
    gamma = Prob[ X = 1 | no Aberration ] and
    delta = Prob[ X = 0 |  Aberration ]}
  \item{partition}{ncomponent vector holding distinct edge labels}
  \item{parent}{like partition, but holds labels for parent edges
    0 is code for root node}
  \item{neutral}{logical vector indicating neutral components}
  \item{data}{ntumor x nswitch matrix of binaries 0/1 indicating abnormality}
}
\value{
    log joint probability density 
}
